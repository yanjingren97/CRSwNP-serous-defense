options(stringsAsFactors = FALSE)
options(error = function() { traceback(2); q(status = 1) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdist)
  library(ggbeeswarm)
  library(patchwork)
  library(ragg)
  library(scales)
  library(grid)
})

root <- "."
out <- Sys.getenv(
  "FIG7_HIGHIMPACT_OUT",
  unset = file.path(root, "deliverables", "Figure7_highimpact_v2_90mm_and_180mm")
)
tiff_dir <- file.path(out, "01_TIFF")
qa_dir <- file.path(out, "02_QA_final_size_96dpi")
code_dir <- file.path(out, "03_R_code")
data_dir <- file.path(out, "04_plot_data")
dir.create(tiff_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#23313D"
navy <- "#173F5F"
blue <- "#245D78"
teal <- "#2A948C"
orange <- "#D9793D"
grey <- "#74828C"
midgrey <- "#AEBAC1"
light <- "#DDE5E9"
pale <- "#F5F7F8"

theme_sci <- function(base_size = 8.2) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      axis.line = element_line(colour = ink, linewidth = 0.32),
      axis.ticks = element_line(colour = ink, linewidth = 0.28),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(colour = ink, face = "bold", size = 8.5),
      plot.title = element_blank(),
      plot.subtitle = element_text(colour = ink, face = "bold", size = 8.0),
      plot.caption = element_text(colour = grey, size = 7.2, hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(colour = ink, size = 7.2),
      legend.key = element_blank(),
      panel.grid.major = element_line(colour = "#EDF1F3", linewidth = 0.24),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3.0, 3.0, 3.0, 3.0, "mm")
    )
}

save_raster <- function(plot, stem, width_mm, height_mm) {
  ragg::agg_tiff(
    file.path(tiff_dir, paste0(stem, ".tiff")),
    width = width_mm, height = height_mm, units = "mm", res = 600,
    compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(qa_dir, paste0(stem, "_QA96.png")),
    width = width_mm, height = height_mm, units = "mm", res = 96,
    background = "white"
  )
  print(plot)
  dev.off()
}

fmt_signed <- function(x, digits = 3) {
  sub("^-", "−", sprintf(paste0("%+.", digits, "f"), x))
}

hedges_g <- function(x_crs, x_ctrl) {
  n1 <- length(x_crs)
  n0 <- length(x_ctrl)
  pooled_sd <- sqrt(((n1 - 1) * var(x_crs) + (n0 - 1) * var(x_ctrl)) /
                      (n1 + n0 - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (n1 + n0 - 2) - 1)
  correction * (mean(x_crs) - mean(x_ctrl)) / pooled_sd
}

# -----------------------------------------------------------------------------
# Figure 7A: residual ECM coupling, with two square cohort panels.
# -----------------------------------------------------------------------------
ecm <- read.csv(file.path(
  root, "deliverables", "section_3_6_boundary_evidence",
  "ecm_coupling_residualized_sample_scores.csv"
))
ecm_summary <- read.csv(file.path(
  root, "results", "bulk_crs_validation", "serous_ecm_coupling.csv"
))
ecm_summary <- ecm_summary[ecm_summary$group == "within_group_residual", ]
ecm$dataset <- factor(ecm$dataset, levels = c("GSE36830", "GSE136825"))
ecm_summary$dataset <- factor(ecm_summary$dataset, levels = levels(ecm$dataset))
ecm_summary$label <- paste0(
  "n = ", ecm_summary$n,
  "  |  ρ = ", sprintf("%.3f", ecm_summary$rho),
  "  |  p = ", sprintf("%.3f", ecm_summary$p)
)
ecm_limit <- max(abs(c(ecm$serous_residual, ecm$ecm_residual)), na.rm = TRUE) * 1.07

pA <- ggplot(ecm, aes(serous_residual, ecm_residual)) +
  geom_hline(yintercept = 0, colour = midgrey, linewidth = 0.34) +
  geom_vline(xintercept = 0, colour = midgrey, linewidth = 0.34) +
  geom_point(colour = navy, size = 1.45, alpha = 0.58, stroke = 0) +
  geom_label(
    data = ecm_summary,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.05, vjust = 1.15,
    family = font_family, fontface = "plain", size = 2.45,
    label.padding = unit(0.65, "mm"), linewidth = 0,
    fill = alpha("white", 0.88), colour = ink
  ) +
  facet_wrap(~dataset, nrow = 1) +
  coord_fixed(ratio = 1, xlim = c(-ecm_limit, ecm_limit),
              ylim = c(-ecm_limit, ecm_limit), clip = "off") +
  labs(
    x = "Residualized serous-defense score",
    y = "Residualized 11-gene ECM score"
  ) +
  theme_sci(8.3) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(5.0, "mm"),
    axis.text = element_text(size = 7.2),
    axis.title = element_text(size = 8.0),
    strip.text = element_text(size = 8.6),
    plot.margin = margin(2.5, 3.0, 2.5, 3.0, "mm")
  )

# -----------------------------------------------------------------------------
# Figure 7B: GeoMx raw patient data plus an estimation panel.
# The point estimate is the frozen Hedges' g; 95% CIs use its conventional
# sampling variance. Exact rank-test p values remain separate annotations.
# -----------------------------------------------------------------------------
spatial <- read.csv(file.path(root, "results", "spatial_validation", "patient_module_scores.csv"))
spatial$group <- ifelse(spatial$group == "CTRL", "Control", "CRSwNP")
spatial$group <- factor(spatial$group, levels = c("Control", "CRSwNP"))
spatial$compartment <- factor(spatial$compartment, levels = c("PanCK", "CD45"))
sp_summary <- read.csv(file.path(root, "results", "spatial_validation", "spatial_module_summary.csv"))
sp_summary$compartment <- factor(sp_summary$compartment, levels = levels(spatial$compartment))

sp_effect <- sp_summary
sp_effect$effect_variance <- with(
  sp_effect,
  (n_crs + n_ctrl) / (n_crs * n_ctrl) +
    hedges_g_crs_minus_ctrl^2 / (2 * (n_crs + n_ctrl - 2))
)
sp_effect$effect_se <- sqrt(sp_effect$effect_variance)
sp_effect$ci_low <- sp_effect$hedges_g_crs_minus_ctrl - qnorm(0.975) * sp_effect$effect_se
sp_effect$ci_high <- sp_effect$hedges_g_crs_minus_ctrl + qnorm(0.975) * sp_effect$effect_se
sp_effect$compartment <- factor(sp_effect$compartment, levels = levels(spatial$compartment))
sp_effect$p_label <- paste0("exact p = ", sprintf("%.3f", sp_effect$p_exact))

pB_raw <- ggplot(spatial, aes(group, score, colour = group)) +
  geom_hline(yintercept = 0, colour = midgrey, linewidth = 0.34) +
  ggbeeswarm::geom_quasirandom(
    width = 0.15, size = 2.25, alpha = 0.82, stroke = 0,
    groupOnX = TRUE
  ) +
  stat_summary(
    fun = median, geom = "crossbar", aes(group = group),
    width = 0.48, linewidth = 0.78, colour = ink
  ) +
  facet_wrap(~compartment, nrow = 1) +
  scale_colour_manual(values = c("Control" = navy, "CRSwNP" = orange), guide = "none") +
  scale_x_discrete(labels = c("Control\nn = 4", "CRSwNP\nn = 6")) +
  labs(
    subtitle = "GeoMx proxy · 4/40 locked genes measured",
    x = NULL,
    y = "Four-gene spatial proxy score"
  ) +
  theme_sci(8.2) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 7.1, face = "bold"),
    axis.text.y = element_text(size = 7.0),
    axis.title.y = element_text(size = 7.8),
    strip.text = element_text(size = 8.4),
    plot.subtitle = element_text(size = 7.2, colour = grey, face = "plain", hjust = 1),
    panel.spacing.x = unit(4.0, "mm"),
    plot.margin = margin(2.5, 3.0, 0.5, 3.0, "mm")
  )

effect_lim <- range(c(sp_effect$ci_low, sp_effect$ci_high), finite = TRUE)
effect_pad <- diff(effect_lim) * 0.18
pB_effect <- ggplot(sp_effect, aes(hedges_g_crs_minus_ctrl, compartment)) +
  geom_vline(xintercept = 0, colour = midgrey, linewidth = 0.42) +
  ggdist::geom_pointinterval(
    aes(xmin = ci_low, xmax = ci_high),
    point_size = 2.7, interval_size = 0.85,
    colour = teal, fill = teal
  ) +
  geom_text(
    aes(x = Inf, label = p_label), hjust = 1.03,
    family = font_family, size = 2.35, colour = ink
  ) +
  scale_x_continuous(
    limits = c(effect_lim[1] - effect_pad, effect_lim[2] + effect_pad * 2.2),
    breaks = pretty_breaks(4)
  ) +
  labs(
    x = "Hedges’ g (CRSwNP − control)\nWald 95% CI",
    y = NULL,
    caption = "Patient-level GeoMx proxy\nExact rank-test p shown separately"
  ) +
  theme_sci(8.0) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.2, face = "bold"),
    axis.text.x = element_text(size = 6.9),
    axis.title.x = element_text(size = 7.5),
    plot.caption = element_text(size = 6.7),
    plot.margin = margin(1.0, 3.0, 2.5, 3.0, "mm")
  )

pB <- pB_raw / pB_effect + patchwork::plot_layout(heights = c(1.8, 1.0))

# -----------------------------------------------------------------------------
# Figure 7C: compact two-line proteomic coverage strip plus descriptive points.
# -----------------------------------------------------------------------------
protein_summary <- read.csv(file.path(
  root, "results", "protein_validation", "PXD013330_module_summary.csv"
))
protein_scores <- read.csv(file.path(
  root, "results", "protein_validation", "PXD013330_sample_module_scores.csv"
))
protein_scores$group <- factor(
  protein_scores$group,
  levels = c("CON", "CRS", "CRSwNP"),
  labels = c("Control", "CRS", "CRSwNP")
)

coverage_counts <- data.frame(
  metric = factor(
    c("Detected ≥1 column", "Complete in all 9"),
    levels = c("Complete in all 9", "Detected ≥1 column")
  ),
  count = c(
    protein_summary$module_genes_detected_any,
    protein_summary$module_genes_complete_9_samples
  )
)
coverage_strip <- merge(
  expand.grid(metric = levels(coverage_counts$metric), index = seq_len(40)),
  coverage_counts,
  by = "metric", all.x = TRUE
)
coverage_strip$covered <- coverage_strip$index <= coverage_strip$count
coverage_strip$metric <- factor(coverage_strip$metric, levels = levels(coverage_counts$metric))

pC_coverage <- ggplot(coverage_strip, aes(index, metric, fill = covered)) +
  geom_tile(width = 0.76, height = 0.58, colour = "white", linewidth = 0.16) +
  geom_text(
    data = coverage_counts, aes(x = 42.0, y = metric, label = paste0(count, "/40")),
    inherit.aes = FALSE, hjust = 1, family = font_family,
    fontface = "bold", size = 2.35, colour = ink
  ) +
  scale_fill_manual(values = c("TRUE" = teal, "FALSE" = "#E5EBEE"), guide = "none") +
  scale_x_continuous(limits = c(0.5, 42.3), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = font_family, base_size = 8) +
  theme(
    axis.text.y = element_text(colour = ink, face = "bold", size = 7.0, hjust = 1),
    plot.margin = margin(3.0, 3.0, 0.5, 3.0, "mm")
  )

pC_scores <- ggplot(protein_scores, aes(group, score, colour = group)) +
  geom_hline(yintercept = 0, colour = midgrey, linewidth = 0.34) +
  ggbeeswarm::geom_quasirandom(
    width = 0.14, size = 2.75, alpha = 0.88, stroke = 0,
    groupOnX = TRUE
  ) +
  stat_summary(
    fun = median, geom = "crossbar", aes(group = group),
    width = 0.48, linewidth = 0.85, colour = ink
  ) +
  scale_colour_manual(
    values = c("Control" = navy, "CRS" = teal, "CRSwNP" = orange),
    guide = "none"
  ) +
  scale_x_discrete(labels = c(
    "Control\n3 columns", "CRS\n3 columns", "CRSwNP\n3 columns"
  )) +
  labs(
    subtitle = "PXD013330 · pooled technical columns\nDescriptive only",
    x = NULL,
    y = "Seven-protein proxy score",
    caption = paste0(
      "Medians shown; Δ(CRSwNP − control) = ",
      fmt_signed(protein_summary$median_difference, 3),
      "\nPooled technical columns\nDescriptive only; no patient-level inference"
    )
  ) +
  theme_sci(8.2) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 7.0, face = "bold"),
    axis.text.y = element_text(size = 7.0),
    axis.title.y = element_text(size = 7.8),
    plot.subtitle = element_text(size = 7.8, face = "bold"),
    plot.caption = element_text(size = 6.6),
    plot.margin = margin(1.0, 3.0, 2.5, 3.0, "mm")
  )

pC <- pC_coverage / pC_scores + patchwork::plot_layout(heights = c(0.42, 1.0))

# Individual 90-mm panels for manual assembly.
save_raster(pA, "Fig7A_ECM_residual_correlations_90mm", 90, 54)
save_raster(pB, "Fig7B_GeoMx_patient_estimation_90mm", 90, 104)
save_raster(pC, "Fig7C_PXD013330_descriptive_90mm", 90, 104)

# Full 180-mm figure: A across the top, B and C below.
p_full <- pA / (pB | pC) +
  patchwork::plot_layout(heights = c(0.78, 1.0)) &
  theme(plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm"))
save_raster(p_full, "Figure7_complete_180mm", 180, 184)

# Lettered full figure for direct manuscript placement. Wrapping each composite
# panel makes patchwork assign exactly three tags rather than tagging the
# internal raw-data and estimation subplots separately.
p_full_abc <-
  patchwork::wrap_elements(full = pA) /
  (patchwork::wrap_elements(full = pB) | patchwork::wrap_elements(full = pC)) +
  patchwork::plot_layout(heights = c(0.78, 1.0)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm"),
    plot.tag = element_text(
      family = font_family, face = "bold", size = 14,
      colour = ink, hjust = 0, vjust = 1
    ),
    plot.tag.position = c(0, 1)
  )
save_raster(p_full_abc, "Figure7_complete_180mm_ABC", 180, 184)

write.csv(ecm, file.path(data_dir, "Fig7A_ECM_residualized_sample_scores.csv"), row.names = FALSE)
write.csv(ecm_summary, file.path(data_dir, "Fig7A_ECM_statistics.csv"), row.names = FALSE)
write.csv(spatial, file.path(data_dir, "Fig7B_GeoMx_patient_scores.csv"), row.names = FALSE)
write.csv(sp_effect, file.path(data_dir, "Fig7B_GeoMx_effects_95CI.csv"), row.names = FALSE)
write.csv(coverage_strip, file.path(data_dir, "Fig7C_PXD013330_coverage_strip.csv"), row.names = FALSE)
write.csv(protein_scores, file.path(data_dir, "Fig7C_PXD013330_technical_columns.csv"), row.names = FALSE)

notes <- c(
  "FIGURE 7 HIGH-IMPACT REDESIGN",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Individual panels: 90 mm wide. Combined figure: 180 mm wide.",
  "All TIFF files: RGB, 600 dpi, LZW compression, Arial typography.",
  "Figure 7A: no confidence ellipse; rho and p are Spearman statistics after group residualization.",
  "Figure 7B: all available patient points are displayed; Hedges' g is the frozen effect estimate.",
  "Figure 7B intervals: Wald 95% CIs from the conventional Hedges' g sampling variance.",
  "GeoMx coverage: four of 40 locked genes; the plotted score is a four-gene proxy.",
  "Figure 7C: PXD013330 columns are pooled technical measurements, not independent patients.",
  "No patient-level p value is shown for PXD013330; the panel is descriptive only.",
  "Individual panels and the standard combined export remain unlettered.",
  "A separate 180-mm complete export embeds A/B/C in bold 14-pt Arial.",
  "",
  capture.output(sessionInfo())
)
writeLines(notes, file.path(out, "R_sessionInfo_and_design_notes.txt"))

file.copy(
  from = file.path(root, "code", "export_figure7_highimpact_v2.R"),
  to = file.path(code_dir, "export_figure7_highimpact_v2.R"),
  overwrite = TRUE
)

message(normalizePath(out, winslash = "/", mustWork = TRUE))
