# Interactive dashboard for the pipeline results (isolates or controls), runs locally in the browser.
#   install.packages(c("shiny", "DT", "plotly", "bslib", "readr", "dplyr", "tidyr"))   # once
#   Rscript -e 'shiny::runApp("dashboard", launch.browser = TRUE)'                     # from the repo folder
# Point it at a summary folder (results/_combined/summary or results_controls/_combined/summary,
# or a copy downloaded from BioCloud) in the box on the left.

library(shiny); library(bslib); library(plotly); library(DT)
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})

candidates <- c(Sys.getenv("BGC_SUMMARY"), "../results_controls/_combined/summary", "../results/_combined/summary",
                "../controls_summary", "../summary")
start_dir <- c(candidates[candidates != "" & file.exists(file.path(candidates, "master_summary.tsv"))], "")[1]

novelty_col <- c(known = "#2271B2", related = "#E69F00", novel_candidate = "#D55E00")
tool_lab <- c(antismash = "antiSMASH", gecco = "GECCO", deepbgc = "DeepBGC", bgc_prophet = "BGC-Prophet")
read_tab <- function(dir, f) {
  p <- file.path(dir, f)
  if (file.exists(p)) read_tsv(p, col_types = cols(.default = "c"), progress = FALSE, na = c("", "NA")) else NULL
}
num <- function(x) suppressWarnings(as.numeric(x))
loci_cols <- c("priority_score", "novelty", "isolate", "locus_id", "contig", "start", "end", "length_kb",
               "n_tools", "tools", "antismash_type", "deepbgc_activity", "closest_MIBiG", "closest_MIBiG_compound",
               "MIBiG_similarity_pct", "resistance_gene_in_locus", "bigscape_family", "contig_edge",
               "gtdb_classification", "possible_new_species")

ui <- page_sidebar(
  title = "BGC pipeline results",
  theme = bs_theme(version = 5, primary = "#2271B2"),
  sidebar = sidebar(
    width = 320,
    textInput("dir", "Summary folder", value = start_dir, placeholder = "…/_combined/summary"),
    actionButton("load", "Load", class = "btn-primary"),
    textOutput("status"),
    hr(),
    selectizeInput("iso", "Isolates", choices = NULL, multiple = TRUE, options = list(placeholder = "all")),
    checkboxGroupInput("nov", "Novelty", choices = names(novelty_col), selected = names(novelty_col)),
    sliderInput("score", "Minimum priority score", min = -1, max = 12, value = -1, step = 1),
    sliderInput("ntools", "Found by at least … tools", min = 1, max = 4, value = 1, step = 1),
    checkboxInput("resist", "Only loci with a resistance gene", FALSE),
    textInput("search", "Search (compound, class, MIBiG…)", "")
  ),
  navset_card_tab(
    nav_panel("BGC loci",
      plotlyOutput("scatter", height = "380px"),
      DTOutput("loci")),
    nav_panel("Per isolate",
      plotlyOutput("per_isolate", height = "420px"),
      plotlyOutput("per_tool", height = "380px")),
    nav_panel("Genomes (QC + taxonomy)", DTOutput("genomes")),
    nav_panel("Controls",
      uiOutput("ctrl_note"),
      plotlyOutput("recovery", height = "380px"),
      DTOutput("ctrl_loci"),
      DTOutput("ctrl_expected"))
  )
)

server <- function(input, output, session) {
  data <- eventReactive(input$load, {
    d <- trimws(path.expand(input$dir))
    validate(need(file.exists(file.path(d, "master_summary.tsv")),
                  paste("No master_summary.tsv in", d)))
    master <- read_tab(d, "master_summary.tsv") |>
      mutate(across(c(priority_score, n_tools, length_kb, MIBiG_similarity_pct), num))
    list(dir = d, master = master, genomes = read_tab(d, "genome_summary.tsv"),
         recovery = read_tab(d, "controls_recovery.tsv"), ctrl_loci = read_tab(d, "controls_loci.tsv"),
         expected = read_tab(d, "controls_expected_bgcs.tsv"))
  }, ignoreNULL = FALSE)

  observeEvent(data(), {
    m <- data()$master
    updateSelectizeInput(session, "iso", choices = sort(unique(m$isolate)), server = TRUE)
    updateSliderInput(session, "score", min = min(m$priority_score, na.rm = TRUE),
                      max = max(m$priority_score, na.rm = TRUE), value = min(m$priority_score, na.rm = TRUE))
  })
  output$status <- renderText(sprintf("%d loci, %d isolates", nrow(data()$master), n_distinct(data()$master$isolate)))

  filtered <- reactive({
    m <- data()$master
    if (length(input$iso)) m <- filter(m, isolate %in% input$iso)
    m <- filter(m, novelty %in% input$nov, priority_score >= input$score, n_tools >= input$ntools)
    if (input$resist) m <- filter(m, !resistance_gene_in_locus %in% c("none", NA))
    if (nzchar(input$search)) {
      hit <- apply(m, 1, function(r) any(grepl(input$search, r, ignore.case = TRUE, fixed = TRUE)))
      m <- m[hit, ]
    }
    m
  })

  output$scatter <- renderPlotly({
    m <- filtered()
    validate(need(nrow(m) > 0, "No loci match the filters"))
    m <- m |> mutate(sim = coalesce(MIBiG_similarity_pct, 0),
                     hover = paste0("<b>", isolate, "</b> ", locus_id, "<br>", contig, ":", start, "-", end,
                                    " (", length_kb, " kb)<br>tools: ", tools, "<br>type: ", antismash_type,
                                    "<br>MIBiG: ", closest_MIBiG_compound, " (", closest_MIBiG, ", ",
                                    MIBiG_similarity_pct, " %)<br>resistance: ", resistance_gene_in_locus))
    plot_ly(m, x = ~sim, y = ~jitter(priority_score, 0.4), color = ~novelty, colors = novelty_col,
            size = ~n_tools, sizes = c(8, 22), text = ~hover, hoverinfo = "text",
            type = "scatter", mode = "markers", marker = list(opacity = 0.75)) |>
      layout(xaxis = list(title = "Similarity to closest MIBiG BGC (%)"),
             yaxis = list(title = "Priority score"), legend = list(title = list(text = "Novelty")),
             shapes = list(list(type = "line", x0 = 0, x1 = 100, y0 = 5, y1 = 5,
                                line = list(dash = "dash", color = "grey"))))
  })

  output$loci <- renderDT({
    m <- filtered() |> select(any_of(loci_cols)) |> arrange(desc(priority_score))
    datatable(m, filter = "top", rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", buttons = c("csv", "excel"))) |>
      formatStyle("novelty", color = styleEqual(names(novelty_col), novelty_col), fontWeight = "bold")
  })

  output$per_isolate <- renderPlotly({
    m <- filtered() |> count(isolate, novelty)
    validate(need(nrow(m) > 0, "No loci match the filters"))
    plot_ly(m, y = ~isolate, x = ~n, color = ~novelty, colors = novelty_col, type = "bar", orientation = "h") |>
      layout(barmode = "stack", xaxis = list(title = "BGC loci"), yaxis = list(title = "", categoryorder = "total ascending"),
             title = "Loci per isolate (filtered)")
  })

  output$per_tool <- renderPlotly({
    m <- filtered() |> separate_rows(tools, sep = ";") |> filter(tools %in% names(tool_lab)) |>
      count(isolate, tool = tool_lab[tools])
    validate(need(nrow(m) > 0, ""))
    plot_ly(m, x = ~isolate, y = ~n, color = ~tool, type = "bar") |>
      layout(barmode = "group", yaxis = list(title = "Loci the tool contributes to"), xaxis = list(title = ""),
             title = "Which tools found the loci")
  })

  output$genomes <- renderDT({
    g <- data()$genomes
    validate(need(!is.null(g), "No genome_summary.tsv in this folder"))
    datatable(g, filter = "top", rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("csv", "excel")))
  })

  output$ctrl_note <- renderUI({
    if (is.null(data()$recovery))
      p("No control tables here. For a control run: bash setup/bgc_pipeline.sh check-controls, then load ",
        code("results_controls/_combined/summary"), ".")
  })
  output$recovery <- renderPlotly({
    r <- data()$recovery
    validate(need(!is.null(r), ""))
    long <- r |> pivot_longer(c(names(tool_lab), "any_tool"), names_to = "tool", values_to = "k") |>
      mutate(k = num(k), n = num(expected), pct = 100 * k / n,
             tool = factor(c(tool_lab, any_tool = "Any tool")[tool], levels = c(tool_lab, "Any tool")))
    plot_ly(long, x = ~tool, y = ~isolate, z = ~pct, type = "heatmap", colors = c("#F2F5F9", "#2271B2"),
            zmin = 0, zmax = 100, text = ~paste0(k, "/", n), texttemplate = "%{text}", hoverinfo = "text") |>
      layout(title = "Known MIBiG BGCs recovered (positive controls)", xaxis = list(title = ""), yaxis = list(title = ""))
  })
  output$ctrl_loci <- renderDT({
    validate(need(!is.null(data()$ctrl_loci), ""))
    datatable(data()$ctrl_loci, rownames = FALSE, options = list(pageLength = 15, dom = "t")) |>
      formatStyle("verdict", color = styleEqual(c("PASS", "CHECK", "NOT_IN_RUN"), c("#009E73", "#D55E00", "grey")),
                  fontWeight = "bold")
  })
  output$ctrl_expected <- renderDT({
    validate(need(!is.null(data()$expected), ""))
    datatable(data()$expected, filter = "top", rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })
}

shinyApp(ui, server)
