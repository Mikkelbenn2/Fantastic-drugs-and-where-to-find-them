# Interactive dashboard for the pipeline results: BGC loci, QC, metadata and controls.
# On BioCloud (reads the results in place):   bash setup/bgc_pipeline.sh dashboard
# On your own computer (repo folder):         Rscript -e 'shiny::runApp("dashboard", launch.browser = TRUE)'
#   needs: install.packages(c("shiny", "DT", "plotly", "bslib", "readr", "dplyr", "tidyr"))
# Data: every <folder>/_combined/summary with a master_summary.tsv inside the project folder, plus
#   metadata.tsv            your own sample information (first column "isolate", any other columns)
#   controls/control_panel.tsv  strain, accession and papers of the controls (added automatically)

library(shiny); library(bslib); library(plotly); library(DT)
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})

project <- normalizePath(Sys.getenv("BGC_PROJECT", ".."), mustWork = FALSE)
find_summaries <- function() {
  d <- c(Sys.getenv("BGC_SUMMARY"), Sys.glob(file.path(project, c("*/_combined/summary", "*summary"))))
  d <- unique(d[d != "" & file.exists(file.path(d, "master_summary.tsv"))])
  setNames(d, sub(paste0("^", project, "/?"), "", d))
}

novelty_col <- c(known = "#2271B2", related = "#E69F00", novel_candidate = "#D55E00", low_support = "#BBBBBB")
tool_lab <- c(antismash = "antiSMASH", gecco = "GECCO", deepbgc = "DeepBGC", bgc_prophet = "BGC-Prophet")
num <- function(x) suppressWarnings(as.numeric(x))
read_tab <- function(p) {
  if (!file.exists(p)) return(NULL)
  read_tsv(p, col_types = cols(.default = "c"), progress = FALSE, na = c("", "NA"))
}
qc_num <- c("raw_reads", "raw_bases", "raw_read_N50", "raw_meanQ", "filtered_reads", "filtered_bases",
            "filtered_read_N50", "filtered_meanQ", "genome_length", "contigs", "N50", "GC_percent",
            "circular_contigs", "mean_depth", "completeness", "contamination", "CDS", "closest_ANI",
            "BGCs_antismash", "BGCs_gecco", "BGCs_deepbgc", "BGCs_bgc_prophet", "resistance_genes",
            "BGC_loci_all_tools")
loci_cols <- c("priority_score", "novelty", "isolate", "locus_id", "contig", "start", "end", "length_kb",
               "n_tools", "tools", "antismash_type", "deepbgc_activity", "closest_MIBiG", "closest_MIBiG_compound",
               "MIBiG_similarity_pct", "resistance_gene_in_locus", "bigscape_family", "contig_edge")

# Sample metadata: your metadata.tsv + the control panel, one row per isolate
read_metadata <- function() {
  own  <- read_tab(file.path(project, "metadata.tsv"))
  ctrl <- read_tab(file.path(project, "controls", "control_panel.tsv"))
  if (!is.null(ctrl)) ctrl <- select(ctrl, any_of(c("isolate", "control", "strain", "data_accession", "expected_antibiotics")))
  if (!is.null(own)) names(own)[1] <- "isolate"
  m <- bind_rows(own, ctrl)
  if (nrow(m) == 0) return(tibble(isolate = character()))
  distinct(m, isolate, .keep_all = TRUE)
}

ui <- page_sidebar(
  title = "BGC pipeline results",
  fillable = FALSE,                      # let the page scroll; otherwise stacked tables are squeezed and overlap
  theme = bs_theme(version = 5, primary = "#2271B2"),
  sidebar = sidebar(
    width = 330,
    selectInput("dir", "Results", choices = NULL),
    textInput("dir_other", NULL, placeholder = "or type another summary folder"),
    actionButton("load", "Load / refresh", class = "btn-primary"),
    textOutput("status"),
    hr(),
    selectizeInput("iso", "Isolates", choices = NULL, multiple = TRUE, options = list(placeholder = "all")),
    selectInput("meta_col", "Colour / filter by metadata", choices = c("(none)")),
    selectizeInput("meta_val", NULL, choices = NULL, multiple = TRUE, options = list(placeholder = "all values")),
    checkboxGroupInput("nov", "Novelty", choices = names(novelty_col), selected = setdiff(names(novelty_col), "low_support")),
    sliderInput("score", "Minimum priority score", min = -1, max = 12, value = -1, step = 1),
    sliderInput("ntools", "Found by at least … tools", min = 1, max = 4, value = 1, step = 1),
    checkboxInput("resist", "Only loci with a resistance gene", FALSE),
    textInput("search", "Search (compound, class, MIBiG…)", "")
  ),
  navset_card_tab(
    nav_panel("BGC loci", plotlyOutput("scatter", height = "380px"), DTOutput("loci")),
    nav_panel("Per isolate", plotlyOutput("per_isolate", height = "420px"), plotlyOutput("per_tool", height = "360px")),
    nav_panel("QC",
      layout_columns(col_widths = c(6, 6, 6, 6),
        plotlyOutput("qc_checkm", height = "320px"), plotlyOutput("qc_assembly", height = "320px"),
        plotlyOutput("qc_reads", height = "320px"), plotlyOutput("qc_depth", height = "320px")),
      DTOutput("genomes")),
    nav_panel("Metadata", uiOutput("meta_note"), DTOutput("metadata")),
    nav_panel("Controls", uiOutput("ctrl_note"), plotlyOutput("recovery", height = "380px"),
              h5("Verdict per control"), DTOutput("ctrl_loci"),
              h5("Expected MIBiG BGCs and the tools that found them"), DTOutput("ctrl_expected"))
  )
)

server <- function(input, output, session) {
  observe({
    s <- find_summaries()
    updateSelectInput(session, "dir", choices = s, selected = s[1])
  })

  data <- eventReactive(list(input$load, input$dir), {
    d <- if (nzchar(trimws(input$dir_other))) trimws(path.expand(input$dir_other)) else input$dir
    if (length(d) == 1 && !startsWith(d, "/")) d <- file.path(project, d)   # typed paths are relative to the project
    if (length(d) == 1 && file.exists(file.path(d, "_combined/summary/master_summary.tsv"))) d <- file.path(d, "_combined/summary")
    validate(need(length(d) == 1 && file.exists(file.path(d, "master_summary.tsv")),
                  paste0("No master_summary.tsv in ", d, ". The run is not finished, or its last (combined) job failed: ",
                         "check squeue --me and <results folder>/slurm_logs/combined_*.out")))
    meta <- read_metadata()
    genomes <- read_tab(file.path(d, "genome_summary.tsv"))
    if (!is.null(genomes)) genomes <- genomes |> mutate(across(any_of(qc_num), num)) |>
      left_join(meta, by = "isolate", suffix = c("", ".meta"))
    master <- read_tab(file.path(d, "master_summary.tsv")) |>
      mutate(across(c(priority_score, n_tools, length_kb, MIBiG_similarity_pct), num)) |>
      left_join(meta, by = "isolate", suffix = c("", ".meta"))
    list(dir = d, master = master, genomes = genomes, meta = meta,
         recovery = read_tab(file.path(d, "controls_recovery.tsv")),
         ctrl_loci = read_tab(file.path(d, "controls_loci.tsv")),
         expected = read_tab(file.path(d, "controls_expected_bgcs.tsv")))
  })

  observeEvent(data(), {
    m <- data()$master
    updateSelectizeInput(session, "iso", choices = sort(unique(m$isolate)), server = TRUE)
    rng <- range(m$priority_score, na.rm = TRUE)
    updateSliderInput(session, "score", min = rng[1], max = rng[2], value = rng[1])
    updateSelectInput(session, "meta_col", choices = c("(none)", setdiff(names(data()$meta), "isolate")))
  })
  observeEvent(input$meta_col, {
    vals <- if (input$meta_col == "(none)") character() else sort(unique(na.omit(data()$meta[[input$meta_col]])))
    updateSelectizeInput(session, "meta_val", choices = vals, selected = character())
  })
  output$status <- renderText(sprintf("%d loci, %d isolates", nrow(data()$master), n_distinct(data()$master$isolate)))

  # isolates allowed by the isolate + metadata filters (used by every tab)
  keep_iso <- reactive({
    iso <- unique(c(data()$master$isolate, data()$genomes$isolate))
    if (length(input$iso)) iso <- intersect(iso, input$iso)
    if (input$meta_col != "(none)" && length(input$meta_val))
      iso <- intersect(iso, data()$meta$isolate[data()$meta[[input$meta_col]] %in% input$meta_val])
    iso
  })
  group_of <- function(df) if (input$meta_col != "(none)" && input$meta_col %in% names(df)) df[[input$meta_col]] else NULL

  filtered <- reactive({
    m <- data()$master |>
      filter(isolate %in% keep_iso(), novelty %in% input$nov, priority_score >= input$score, n_tools >= input$ntools)
    if (input$resist) m <- filter(m, !resistance_gene_in_locus %in% c("none", NA))
    if (nzchar(input$search))
      m <- m[apply(m, 1, function(r) any(grepl(input$search, r, ignore.case = TRUE, fixed = TRUE))), ]
    m
  })
  genomes <- reactive({
    g <- data()$genomes
    validate(need(!is.null(g), "No genome_summary.tsv in this folder"))
    filter(g, isolate %in% keep_iso())
  })

  # ---- BGC loci ----
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
      layout(xaxis = list(title = "Similarity to closest MIBiG BGC (%)"), yaxis = list(title = "Priority score"),
             legend = list(title = list(text = "Novelty")),
             shapes = list(list(type = "line", x0 = 0, x1 = 100, y0 = 5, y1 = 5, line = list(dash = "dash", color = "grey"))))
  })
  output$loci <- renderDT({
    m <- filtered() |> select(any_of(c(loci_cols, setdiff(names(data()$meta), "isolate")))) |> arrange(desc(priority_score))
    datatable(m, filter = "top", rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", buttons = c("csv", "excel"))) |>
      formatStyle("novelty", color = styleEqual(names(novelty_col), novelty_col), fontWeight = "bold")
  })

  # ---- per isolate ----
  output$per_isolate <- renderPlotly({
    m <- filtered() |> count(isolate, novelty)
    validate(need(nrow(m) > 0, "No loci match the filters"))
    plot_ly(m, y = ~isolate, x = ~n, color = ~novelty, colors = novelty_col, type = "bar", orientation = "h") |>
      layout(barmode = "stack", title = "Loci per isolate (filtered)", xaxis = list(title = "BGC loci"),
             yaxis = list(title = "", categoryorder = "total ascending"))
  })
  output$per_tool <- renderPlotly({
    m <- filtered() |> separate_rows(tools, sep = ";") |> filter(tools %in% names(tool_lab)) |>
      count(isolate, tool = tool_lab[tools])
    validate(need(nrow(m) > 0, ""))
    plot_ly(m, x = ~isolate, y = ~n, color = ~tool, type = "bar") |>
      layout(barmode = "group", title = "Which tools found the loci", xaxis = list(title = ""),
             yaxis = list(title = "Loci the tool contributes to"))
  })

  # ---- QC ----
  qc_plot <- function(g, x, y, xlab, ylab, title, shapes = list(), logx = FALSE) {
    grp <- group_of(g)
    g$hover <- paste0("<b>", g$isolate, "</b><br>", xlab, ": ", g[[x]], "<br>", ylab, ": ", g[[y]])
    p <- if (is.null(grp)) plot_ly(g, x = g[[x]], y = g[[y]], text = ~hover, marker = list(color = "#2271B2", size = 10))
         else plot_ly(g, x = g[[x]], y = g[[y]], color = grp, text = ~hover, marker = list(size = 10))
    p |> add_markers(hoverinfo = "text") |>
      layout(title = title, xaxis = list(title = xlab, type = if (logx) "log" else "linear"),
             yaxis = list(title = ylab), shapes = shapes)
  }
  output$qc_checkm <- renderPlotly({
    g <- genomes()
    qc_plot(g, "completeness", "contamination", "Completeness (%)", "Contamination (%)", "CheckM2 (guides: 95 % / 5 %)",
            shapes = list(list(type = "line", x0 = 95, x1 = 95, y0 = 0, y1 = max(c(g$contamination, 6), na.rm = TRUE),
                               line = list(dash = "dash", color = "grey")),
                          list(type = "line", x0 = min(c(g$completeness, 80), na.rm = TRUE), x1 = 100, y0 = 5, y1 = 5,
                               line = list(dash = "dash", color = "grey"))))
  })
  output$qc_assembly <- renderPlotly(qc_plot(genomes(), "genome_length", "contigs", "Genome length (bp)", "Contigs",
                                             "Assembly (QUAST)", logx = TRUE))
  output$qc_reads <- renderPlotly({
    g <- genomes()
    validate(need(any(!is.na(g$raw_read_N50)), "No read statistics (genome-only inputs, e.g. controls)"))
    qc_plot(g, "raw_read_N50", "raw_meanQ", "Read N50 (bp)", "Mean read quality (Q)", "Reads (SeqKit)")
  })
  output$qc_depth <- renderPlotly({
    g <- genomes() |> filter(!is.na(mean_depth))
    validate(need(nrow(g) > 0, "No assembly depth (genome-only inputs, e.g. controls)"))
    plot_ly(g, x = ~isolate, y = ~mean_depth, type = "bar", marker = list(color = "#2271B2")) |>
      layout(title = "Mean assembly depth (Flye)", xaxis = list(title = ""), yaxis = list(title = "Depth (x)"))
  })
  output$genomes <- renderDT({
    g <- genomes() |> mutate(QC = if_else(completeness >= 95 & contamination <= 5, "PASS", "CHECK", missing = "?"),
                             .after = isolate)
    datatable(g, filter = "top", rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("csv", "excel"))) |>
      formatStyle("QC", color = styleEqual(c("PASS", "CHECK"), c("#009E73", "#D55E00")), fontWeight = "bold")
  })

  # ---- metadata ----
  output$meta_note <- renderUI({
    if (!file.exists(file.path(project, "metadata.tsv")))
      p("Add your own sample information as ", code(file.path(project, "metadata.tsv")),
        ": tab-separated, first column ", code("isolate"), " (same names as in raw_data/), any other columns ",
        "(e.g. site, soil_type, date, inhibition_assay). It is then shown here, joined to every table, ",
        "and usable in 'Colour / filter by metadata'.")
  })
  output$metadata <- renderDT(datatable(data()$meta |> filter(isolate %in% keep_iso()), filter = "top",
                                        rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)))

  # ---- controls ----
  output$ctrl_note <- renderUI({
    if (is.null(data()$recovery))
      p("No control tables in this folder. Run ", code("bash setup/bgc_pipeline.sh check-controls"),
        " and choose results_controls.")
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
    datatable(data()$ctrl_loci, rownames = FALSE, options = list(pageLength = 15, dom = "t", scrollX = TRUE)) |>
      formatStyle("verdict", color = styleEqual(c("PASS", "CHECK", "NOT_IN_RUN"), c("#009E73", "#D55E00", "grey")),
                  fontWeight = "bold")
  })
  output$ctrl_expected <- renderDT({
    validate(need(!is.null(data()$expected), ""))
    datatable(data()$expected, filter = "top", rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })
}

shinyApp(ui, server)
