#!/usr/bin/env Rscript
# Figures of the positive/negative control run, for the report.
#   1) bash setup/bgc_pipeline.sh controls        (wait until the jobs are done)
#   2) bash setup/bgc_pipeline.sh check-controls  (writes the controls_*.tsv tables)
#   3) Rscript controls/plot_controls.R [SUMMARY_DIR] [OUT_DIR]
# SUMMARY_DIR default: results_controls/_combined/summary   OUT_DIR default: SUMMARY_DIR/figures
# In RStudio: set summ_dir / controls_dir below and press Source.
# Output (PDF + PNG each):
#   fig_controls_overview  A recovery of known BGCs, B loci per control by novelty, C priority scores
#   fig_recovery           A alone (main text)
#   fig_expected_bgcs      every expected MIBiG BGC and which tool found it (supplementary)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
summ_dir     <- if (length(args) >= 1) args[1] else "results_controls/_combined/summary"
out_dir      <- if (length(args) >= 2) args[2] else file.path(summ_dir, "figures")
script       <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
controls_dir <- if (length(script) == 1) dirname(script) else "controls"   # folder with control_panel.tsv

need <- file.path(summ_dir, c("controls_expected_bgcs.tsv", "controls_loci.tsv", "master_summary.tsv"))
if (!all(file.exists(need)))
  stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "),
       "\nRun first: bash setup/bgc_pipeline.sh check-controls")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- style (theme_user from your usual R plots) --------------------------------------------
theme_user <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(colour = "grey40"),
          axis.text = element_text(colour = "grey20"), axis.ticks = element_line(colour = "grey50"),
          legend.title = element_text(face = "bold"), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), plot.margin = margin(8, 10, 6, 8))
}
theme_set(theme_user())
novelty_col <- c(known = "#2271B2", related = "#E69F00", novel_candidate = "#D55E00", novel_other = "#CC79A7",
                 low_support = "#BBBBBB")
tool_lab <- c(antismash = "antiSMASH", gecco = "GECCO", deepbgc = "DeepBGC", bgc_prophet = "BGC-Prophet")

# Short strain names for the axes (falls back to the isolate name for anything else)
strain_lab <- c(
  "pos_Bvelezensis_FZB42" = "B. velezensis FZB42", "pos_Scoelicolor_A3-2" = "S. coelicolor A3(2)",
  "pos_Pprotegens_Pf-5" = "P. protegens Pf-5", "pos_Svenezuelae_ATCC10712" = "S. venezuelae ATCC 10712",
  "pos_Sgriseus_NBRC13350" = "S. griseus NBRC 13350", "pos_Serythraea_NRRL2338" = "S. erythraea NRRL 2338",
  "pos_Ecoli_Nissle1917" = "E. coli Nissle 1917", "neg_Ecoli_K12_MG1655" = "E. coli K-12 (ONT reads)",
  "neg_Ecoli_K12_MG1655_ref" = "E. coli K-12 (reference)", "neg_Buchnera_APS" = "Buchnera APS",
  "neg_Rprowazekii_MadridE" = "R. prowazekii Madrid E", "neg_Cglutamicum_ATCC13032" = "C. glutamicum ATCC 13032")
lab_of <- function(x) ifelse(x %in% names(strain_lab), strain_lab[x], x)

# ---- data ----------------------------------------------------------------------------------
read_tab <- function(f) read_tsv(file.path(summ_dir, f), col_types = cols(.default = "c"), progress = FALSE)
expected <- read_tab("controls_expected_bgcs.tsv")
loci     <- read_tab("controls_loci.tsv")
master   <- read_tab("master_summary.tsv") |> mutate(priority_score = as.numeric(priority_score))

role <- loci |> select(isolate, control) |>
  mutate(group = recode(control, positive = "Positive", negative = "Negative",
                        technical = "Negative", .default = "Other"),
         group = factor(group, levels = c("Positive", "Negative", "Other")))
controls <- role |> filter(group != "Other")        # leaves out any non-control isolates

# ---- A: recovery of the expected MIBiG BGCs ------------------------------------------------
hits <- expected |>
  separate_rows(found_by, sep = ",") |>
  filter(found_by %in% names(tool_lab)) |>
  distinct(isolate, mibig_accession, tool = found_by)
recovery <- expected |>
  distinct(isolate, mibig_accession) |>
  crossing(tool = c(names(tool_lab), "any")) |>
  left_join(hits |> mutate(found = TRUE), by = c("isolate", "mibig_accession", "tool")) |>
  group_by(isolate, mibig_accession) |>
  mutate(found = if_else(tool == "any", any(found %in% TRUE), found %in% TRUE)) |>
  group_by(isolate, tool) |>
  summarise(k = sum(found), n = n(), .groups = "drop") |>
  mutate(pct = 100 * k / n,
         tool = factor(recode(tool, !!!tool_lab, any = "Any tool"), levels = c(tool_lab, "Any tool")),
         strain = fct_rev(factor(lab_of(isolate))))

p_recovery <- ggplot(recovery, aes(tool, strain, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(k, "/", n), colour = pct > 55), size = 3.3, show.legend = FALSE) +
  scale_fill_gradient(low = "#F2F5F9", high = "#2271B2", limits = c(0, 100), name = "Recovered (%)") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  geom_vline(xintercept = 4.5, colour = "grey30", linewidth = 0.5) +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Known BGCs recovered",
       subtitle = "Expected MIBiG BGCs found / expected", x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0))

# ---- B: loci per control, by novelty -------------------------------------------------------
novelty_counts <- master |>
  filter(isolate %in% controls$isolate) |>
  count(isolate, novelty) |>
  right_join(controls |> select(isolate, group), by = "isolate") |>
  mutate(n = replace_na(n, 0L), novelty = replace_na(novelty, "known"),
         novelty = factor(novelty, levels = rev(names(novelty_col))),
         strain = lab_of(isolate))
order_strain <- novelty_counts |> group_by(strain) |> summarise(t = sum(n)) |> arrange(t) |> pull(strain)
verdicts <- loci |> filter(isolate %in% controls$isolate, verdict != "-") |>
  left_join(novelty_counts |> group_by(isolate) |> summarise(total = sum(n)), by = "isolate") |>
  left_join(controls |> select(isolate, group), by = "isolate") |>
  mutate(strain = factor(lab_of(isolate), levels = order_strain))

p_novelty <- ggplot(novelty_counts |> mutate(strain = factor(strain, levels = order_strain)),
                    aes(n, strain, fill = novelty)) +
  geom_col(width = 0.7) +
  geom_text(data = verdicts, aes(total, strain, label = verdict), inherit.aes = FALSE,
            hjust = -0.2, size = 2.9, colour = "grey30") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = novelty_col, breaks = names(novelty_col),
                    labels = c("Known (80 %+)", "Related", "Novel candidate", "Novel, non-antibiotic type", "Low support (1 ML tool)"),
                    name = "Novelty") +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "BGC loci per control", subtitle = "All tools merged; label = verdict",
       x = "Number of BGC loci", y = NULL) +
  theme(legend.position = "bottom")

# ---- C: priority scores ----------------------------------------------------------------------
scores <- master |>
  filter(isolate %in% controls$isolate) |>
  left_join(controls |> select(isolate, group), by = "isolate") |>
  mutate(strain = factor(lab_of(isolate), levels = order_strain),
         novelty = factor(novelty, levels = names(novelty_col)))

p_scores <- ggplot(scores, aes(priority_score, strain, colour = novelty)) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.15, height = 0.2, size = 1.8, alpha = 0.8) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = novelty_col, drop = FALSE,
                      labels = c("Known", "Related", "Novel candidate", "Novel other", "Low support"), name = "Novelty") +
  guides(colour = guide_legend(nrow = 2)) +
  labs(title = "Priority scores", subtitle = "Dashed: score 5 (strong candidate)",
       x = "Priority score", y = NULL) +
  theme(legend.position = "bottom")

# ---- supplementary: every expected BGC × tool --------------------------------------------------
per_bgc <- expected |>
  mutate(label = paste0(substr(sub(";.*", "", compound), 1, 28), " (", mibig_accession, ")")) |>
  crossing(tool = names(tool_lab)) |>
  mutate(found = mapply(function(t, f) t %in% strsplit(f, ",")[[1]], tool, found_by),
         tool = factor(tool_lab[tool], levels = tool_lab), strain = lab_of(isolate))

p_per_bgc <- ggplot(per_bgc, aes(tool, fct_rev(label), fill = found)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  facet_grid(strain ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(`TRUE` = "#2271B2", `FALSE` = "grey90"),
                    labels = c(`TRUE` = "Found", `FALSE` = "Not found"), name = NULL) +
  scale_x_discrete(position = "top") +
  labs(title = "Expected MIBiG BGCs per positive control", x = NULL, y = NULL) +
  theme(plot.title.position = "plot", axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 7),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        legend.position = "bottom")

# ---- save ------------------------------------------------------------------------------------
save_both <- function(p, name, w, h) {
  for (ext in c("pdf", "png"))
    ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = w, height = h, units = "mm",
           dpi = 300, bg = "white")
}
overview <- (p_recovery | (p_novelty / p_scores)) +
  plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "A")
save_both(overview & theme(legend.text = element_text(size = 8)), "fig_controls_overview", 270, 240)
save_both(p_recovery, "fig_recovery", 180, 110)
save_both(p_per_bgc, "fig_expected_bgcs", 180, 40 + 4.2 * nrow(distinct(expected, mibig_accession)))
message("Figures written to ", normalizePath(out_dir))
