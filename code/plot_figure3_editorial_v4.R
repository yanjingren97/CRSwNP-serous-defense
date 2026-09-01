options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(scales)
  library(grid)
})

project_dir <- "."
out_dir <- file.path(project_dir, "deliverables", "figure3_intrinsic_v4_editorial")
raster_dir <- Sys.getenv("FIG3_V4_RASTER_DIR", unset = tempdir())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#263442"
muted <- "#7C8993"
grid_col <- "#DCE4E8"
navy <- "#245873"
teal <- "#258E88"
orange <- "#D96C28"
warm_light <- "#F4D2B5"

theme_editorial <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = ink),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(p, stem, width, height) {
  agg_tiff(file.path(raster_dir, paste0(stem, ".tiff")),
           width = width, height = height, units = "in", res = 600,
           compression = "lzw", background = "white")
  print(p)
  dev.off()

  agg_png(file.path(raster_dir, paste0(stem, ".png")),
          width = width, height = height, units = "in", res = 220,
          background = "white")
  print(p)
  dev.off()

  cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")),
            width = width, height = height, family = "sans")
  print(p)
  dev.off()
}

# -----------------------------------------------------------------------------
# 1. Cell-type specificity as polar evidence glyphs.
# Each petal is scaled to the maximum observed value within that metric.
# This preserves cross-cell-type rank while keeping incomparable units separate.
# -----------------------------------------------------------------------------
spec <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                           "tables", "celltype_crossstudy_specificity.csv"))
spec$facet_label <- spec$celltype

spec_long <- rbind(
  data.frame(celltype = spec$celltype, facet_label = spec$facet_label,
             metric = "Cross-study\nrho", value = spec$spearman_rho),
  data.frame(celltype = spec$celltype, facet_label = spec$facet_label,
             metric = "Direction\nagreement", value = spec$strong_direction_rate),
  data.frame(celltype = spec$celltype, facet_label = spec$facet_label,
             metric = "Strong-gene\nburden", value = spec$strong_both)
)
spec_long$metric <- factor(spec_long$metric,
                           levels = c("Cross-study\nrho", "Direction\nagreement",
                                      "Strong-gene\nburden"))
spec_long$scaled <- ave(spec_long$value, spec_long$metric, FUN = function(z) z / max(z))
spec_long$celltype <- factor(spec_long$celltype,
                             levels = c("Secretory", "Basal", "Fibroblast", "Goblet"))
facet_levels <- spec$facet_label[match(levels(spec_long$celltype), spec$celltype)]
spec_long$facet_label <- factor(spec_long$facet_label, levels = facet_levels)

cell_cols <- c("Secretory" = orange, "Basal" = navy,
               "Fibroblast" = teal, "Goblet" = "#7C8E9A")

p_spec <- ggplot(spec_long, aes(metric, scaled, fill = celltype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.55, alpha = 0.96) +
  geom_hline(yintercept = c(0.25, 0.50, 0.75, 1.00),
             colour = grid_col, linewidth = 0.35) +
  coord_polar(start = -pi / 3, clip = "off") +
  facet_wrap(~facet_label, ncol = 2) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0.5, 1), labels = NULL,
                     expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = "Petal length is scaled within each metric; the outer ring denotes the cohort maximum.") +
  theme_editorial(8.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 7.3, face = "bold", colour = ink),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 9.2, lineheight = 1.15, margin = margin(0, 0, 5, 0)),
    plot.caption = element_text(size = 7.2, colour = muted, hjust = 0.5,
                                margin = margin(8, 0, 0, 0)),
    plot.margin = margin(10, 16, 8, 16)
  )
save_plot(p_spec, "celltype_specificity_polar_profile_v4", 7.2, 6.0)

# -----------------------------------------------------------------------------
# 2. Leave-one-patient-out stability as an empirical polar envelope.
# The ribbon is the min-max influence envelope of all patient omissions relative
# to the full paired cohort; the orange polygon is the full-cohort reference.
# -----------------------------------------------------------------------------
loo <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                          "tables", "leave_one_out_summary.csv"))
full <- loo[loo$omitted_donor == "none_full", ][1, ]
omit <- loo[loo$omitted_donor != "none_full", ]

ratio_frame <- data.frame(
  metric = rep(c("Module-effect\nmagnitude", "External\nrho",
                 "Direction\nagreement", "Strong-gene\nburden"), each = nrow(omit)),
  ratio = c(abs(omit$locked_median_paired_effect) / abs(full$locked_median_paired_effect),
            omit$rho_external / full$rho_external,
            omit$strong_same_direction / full$strong_same_direction,
            omit$strong_genes / full$strong_genes)
)
metric_levels <- c("Module-effect\nmagnitude", "External\nrho",
                   "Direction\nagreement", "Strong-gene\nburden")
envelope <- aggregate(ratio ~ metric, ratio_frame,
                      FUN = function(z) c(min = min(z), max = max(z), median = median(z)))
envelope <- data.frame(
  metric = envelope$metric,
  ymin = envelope$ratio[, "min"],
  ymax = envelope$ratio[, "max"],
  median = envelope$ratio[, "median"]
)
envelope$metric <- factor(envelope$metric, levels = metric_levels)
envelope <- envelope[order(envelope$metric), ]
envelope$x <- 0:3
envelope_closed <- rbind(envelope, transform(envelope[1, ], x = 4))
envelope_closed$median_low <- envelope_closed$median - 0.010
envelope_closed$median_high <- envelope_closed$median + 0.010
reference <- data.frame(x = 0:4, ymin = 0.986, ymax = 1.014)

p_loo <- ggplot(envelope_closed, aes(x)) +
  geom_hline(yintercept = c(0.75, 1.00, 1.25),
             colour = c(grid_col, "#A7B2BA", grid_col),
             linewidth = c(0.35, 0.65, 0.35), linetype = c(1, 2, 1)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, group = 1),
              fill = teal, alpha = 0.30, colour = NA) +
  geom_ribbon(data = reference, aes(x = x, ymin = ymin, ymax = ymax, group = 1),
              inherit.aes = FALSE, fill = orange, alpha = 0.92, colour = NA) +
  geom_ribbon(aes(ymin = median_low, ymax = median_high, group = 1),
              fill = navy, alpha = 0.96, colour = NA) +
  coord_polar(start = -pi / 4, clip = "off") +
  scale_x_continuous(limits = c(0, 4), breaks = 0:3, labels = metric_levels,
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.62, 1.38), breaks = c(0.75, 1, 1.25),
                     labels = c("0.75x", "1.00x", "1.25x"), expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       subtitle = "All patient omissions retained decreases in 40/40 locked genes",
       caption = "Teal ribbon: omission range   |   Navy band: omission median   |   Orange band: full cohort") +
  theme_editorial(9) +
  theme(
    panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.4),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 8.2, face = "bold", colour = ink),
    axis.text.y = element_text(size = 7.2, colour = muted),
    axis.ticks = element_blank(),
    plot.subtitle = element_text(size = 8.7, face = "bold", hjust = 0.5,
                                 colour = ink, margin = margin(0, 0, 10, 0)),
    plot.caption = element_text(size = 7.2, colour = muted, hjust = 0.5,
                                margin = margin(10, 0, 0, 0)),
    plot.margin = margin(12, 20, 10, 20)
  )
save_plot(p_loo, "leave_one_patient_out_stability_envelope_v4", 5.8, 5.8)

# -----------------------------------------------------------------------------
# 3. Distributional meta-evidence plot for representative genes.
# Each filled ridge is the normal sampling distribution implied by the estimate
# and standard error, replacing point-and-whisker forest glyphs.
# -----------------------------------------------------------------------------
models <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                             "tables", "locked_gene_model_estimates.csv"))
representative <- c("STATH", "DMBT1", "AZGP1", "ODAM", "GJC3", "PPP1R1B",
                    "CA2", "KCNN4", "CHRM3", "TCN1", "ANO1", "SOD3")
forest <- models[models$gene %in% representative, ]
forest <- forest[order(forest$meta_effect), ]
gene_levels <- rev(forest$gene)

layers <- list(
  `Paired discovery` = c("effect_paired", "se_paired"),
  `External cohort` = c("effect_external", "se_external"),
  `Fixed meta-effect` = c("meta_effect", "meta_se")
)
ridge_parts <- list()
z <- 1L
for (layer_name in names(layers)) {
  cols <- layers[[layer_name]]
  for (i in seq_len(nrow(forest))) {
    mu <- forest[[cols[1]]][i]
    se <- forest[[cols[2]]][i]
    xs <- seq(mu - 3.5 * se, mu + 3.5 * se, length.out = 180)
    dens <- dnorm(xs, mean = mu, sd = se)
    amp <- 0.37 * dens / max(dens)
    base <- match(forest$gene[i], gene_levels)
    ridge_parts[[z]] <- data.frame(
      gene = forest$gene[i], layer = layer_name,
      x = c(xs, rev(xs)),
      y = c(base + amp, rev(base - amp)),
      order_id = seq_len(2 * length(xs)),
      group = paste(layer_name, forest$gene[i], sep = "__"),
      I2 = forest$I2_percent[i]
    )
    z <- z + 1L
  }
}
ridge <- do.call(rbind, ridge_parts)
ridge$layer <- factor(ridge$layer, levels = names(layers))

layer_cols <- c("Paired discovery" = navy,
                "External cohort" = teal,
                "Fixed meta-effect" = orange)
ridge$fill_col <- layer_cols[as.character(ridge$layer)]
# In the pooled layer, heterogeneity subtly modulates saturation without adding
# another geometric symbol.
meta_rows <- ridge$layer == "Fixed meta-effect"
ridge$fill_col[meta_rows] <- colorRampPalette(c(warm_light, orange))(101)[
  pmin(100, pmax(0, round(ridge$I2[meta_rows]))) + 1
]

p_density <- ggplot(ridge, aes(x, y, group = group, fill = fill_col)) +
  geom_polygon(colour = "white", linewidth = 0.28, alpha = 0.96) +
  geom_vline(xintercept = 0, colour = "#98A4AD", linewidth = 0.55,
             linetype = "22") +
  facet_grid(. ~ layer) +
  scale_fill_identity() +
  scale_y_continuous(breaks = seq_along(gene_levels), labels = gene_levels,
                     expand = expansion(add = c(0.55, 0.55))) +
  scale_x_continuous(limits = c(-12, 1.1),
                     breaks = c(-10, -8, -6, -4, -2, 0),
                     expand = c(0, 0)) +
  labs(x = "Effect distribution: polyp - non-polyp", y = NULL,
       caption = "Ridge width reflects standard error; pooled ridge colour intensity reflects between-study heterogeneity.") +
  theme_editorial(9) +
  theme(
    panel.grid.major.y = element_line(colour = "#E8EDF0", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(colour = ink, linewidth = 0.45),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8.4, face = "bold"),
    strip.text = element_text(size = 9.0, face = "bold", margin = margin(0, 0, 7, 0)),
    panel.spacing.x = unit(4.5, "mm"),
    plot.caption = element_text(size = 7.2, colour = muted, hjust = 0,
                                margin = margin(7, 0, 0, 0)),
    plot.margin = margin(8, 10, 8, 8)
  )
save_plot(p_density, "representative_gene_distributional_meta_v4", 8.2, 5.4)

write.csv(spec_long, file.path(out_dir, "celltype_specificity_polar_profile_data.csv"), row.names = FALSE)
write.csv(ratio_frame, file.path(out_dir, "leave_one_patient_out_ratio_data.csv"), row.names = FALSE)
write.csv(envelope, file.path(out_dir, "leave_one_patient_out_envelope_data.csv"), row.names = FALSE)
write.csv(ridge, file.path(out_dir, "representative_gene_distributional_meta_data.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
