options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(scales)
})

project_dir <- "."
out_dir <- file.path(project_dir, "deliverables", "figure3_intrinsic_v1")
preview_dir <- file.path(out_dir, "previews")
raster_dir <- Sys.getenv("FIG3_RASTER_DIR", unset = tempdir())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#263442"
muted <- "#758596"
grid_col <- "#DCE3E8"
accent <- "#E76F00"
blue <- "#137CBD"
teal <- "#159D8C"

theme_sci <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink),
      axis.line = element_line(colour = ink, linewidth = 0.45),
      axis.ticks = element_line(colour = ink, linewidth = 0.4),
      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.35),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#F2F5F7", colour = NA),
      strip.text = element_text(face = "bold", colour = ink),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.8),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_plot <- function(p, stem, width, height) {
  tiff_file <- file.path(raster_dir, paste0(stem, ".tiff"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(raster_dir, paste0(stem, ".png"))
  agg_tiff(tiff_file, width = width, height = height, units = "in",
           res = 600, compression = "lzw", background = "white")
  print(p)
  dev.off()
  cairo_pdf(pdf_file, width = width, height = height, family = font_family)
  print(p)
  dev.off()
  agg_png(png_file, width = width, height = height, units = "in",
          res = 220, background = "white")
  print(p)
  dev.off()
}

# Frozen non-definition-gene effects in Serous glandular cells.
effects <- read.csv(gzfile(file.path(project_dir, "results", "secretory_subtypes",
                                    "subtype_nondefinition_gene_effects.csv.gz")))
effects <- effects[effects$cell_min == 20 & effects$subtype == "Serous_glandular", ]
effects <- effects[is.finite(effects$effect_paired) & is.finite(effects$effect_external), ]

locked <- read.csv(file.path(project_dir, "results", "locked_40_gene_module.csv"))$gene
effects$locked <- effects$gene %in% locked
effects$strong <- abs(effects$effect_paired) >= 0.5 & abs(effects$effect_external) >= 0.5
effects$same_direction <- sign(effects$effect_paired) == sign(effects$effect_external)

rho <- suppressWarnings(cor(effects$effect_paired, effects$effect_external,
                            method = "spearman", use = "complete.obs"))
strong_n <- sum(effects$strong)
strong_agreement <- mean(effects$same_direction[effects$strong])

# Panel 1: cross-study effect concordance using two-dimensional binning.
locked_effects <- effects[effects$locked, ]
label_genes <- c("STATH", "DMBT1", "AZGP1", "CA2")
lab <- locked_effects[locked_effects$gene %in% label_genes, ]

p_concordance <- ggplot(effects, aes(effect_paired, effect_external)) +
  geom_hline(yintercept = 0, colour = "#AAB2B9", linewidth = 0.45) +
  geom_vline(xintercept = 0, colour = "#AAB2B9", linewidth = 0.45) +
  geom_bin_2d(bins = 72, aes(fill = after_stat(count))) +
  scale_fill_gradientn(
    colours = c("#EEF3F6", "#B8CEDD", "#5D8EAF", "#244A68"),
    trans = "sqrt", name = "Genes/bin"
  ) +
  geom_point(data = locked_effects, aes(effect_paired, effect_external),
             inherit.aes = FALSE, shape = 21, size = 2.15, stroke = 0.45,
             fill = accent, colour = "#1D252C") +
  geom_text(data = lab, aes(effect_paired, effect_external, label = gene),
            inherit.aes = FALSE, family = font_family, size = 2.45,
            colour = "#1D252C", nudge_y = -0.28, check_overlap = TRUE) +
  annotate("label", x = min(effects$effect_paired, na.rm = TRUE) + 0.4,
           y = max(effects$effect_external, na.rm = TRUE) - 0.2,
           hjust = 0, vjust = 1,
           label = sprintf("15,445 non-definition genes\nSpearman rho = %.3f\nStrong genes: %s; concordant = %.1f%%",
                           rho, comma(strong_n), 100 * strong_agreement),
           family = font_family, size = 3.0, colour = ink, fill = alpha("white", 0.9)) +
  coord_cartesian(clip = "off") +
  labs(x = "Paired polyp - ethmoid effect (GSE235711)",
       y = "Polyp - inferior turbinate effect (GSE276503)") +
  theme_sci(9) +
  theme(panel.grid.major.y = element_blank(), legend.position = c(0.83, 0.18),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))
save_plot(p_concordance, "crossstudy_concordance", 6.4, 5.7)

# Panel 2: threshold-free rank-rank hypergeometric overlap landscape.
rrho_grid <- function(x, y, decreasing, label) {
  n <- length(x)
  ord_x <- order(x, decreasing = decreasing, na.last = NA)
  ord_y <- order(y, decreasing = decreasing, na.last = NA)
  steps <- unique(round(seq(0.03, 0.50, length.out = 34) * n))
  rank_y <- integer(n)
  rank_y[ord_y] <- seq_len(n)
  out <- vector("list", length(steps) * length(steps))
  z <- 1L
  for (i in steps) {
    prefix_x <- ord_x[seq_len(i)]
    ry <- rank_y[prefix_x]
    for (j in steps) {
      overlap <- sum(ry <= j)
      p <- phyper(overlap - 1, i, n - i, j, lower.tail = FALSE)
      out[[z]] <- data.frame(
        discovery_pct = 100 * i / n,
        external_pct = 100 * j / n,
        logp = min(-log10(max(p, .Machine$double.xmin)), 120),
        direction = label
      )
      z <- z + 1L
    }
  }
  do.call(rbind, out)
}

rrho <- rbind(
  rrho_grid(effects$effect_paired, effects$effect_external, FALSE, "Concordant decreases"),
  rrho_grid(effects$effect_paired, effects$effect_external, TRUE, "Concordant increases")
)

p_rrho <- ggplot(rrho, aes(discovery_pct, external_pct, fill = logp)) +
  geom_tile() +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_gradientn(colours = c("#F5F7F8", "#BFD6DE", "#4A9A9A", "#F0B85B", "#B94132"),
                       name = expression(-log[10](italic(P))), limits = c(0, 120),
                       oob = squish) +
  scale_x_continuous(expand = c(0, 0), breaks = c(10, 25, 40, 50), labels = label_percent(scale = 1)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(10, 25, 40, 50), labels = label_percent(scale = 1)) +
  labs(x = "Top-ranked genes included: paired discovery",
       y = "Top-ranked genes included: external cohort") +
  theme_sci(9) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), legend.position = "right")
save_plot(p_rrho, "rank_rank_overlap_landscape", 7.0, 3.5)

# Panel 3: identically structured cell-type specificity audit.
spec <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                           "tables", "celltype_crossstudy_specificity.csv"))
spec$celltype <- factor(spec$celltype, levels = c("Goblet", "Fibroblast", "Basal", "Secretory"))
spec_long <- rbind(
  data.frame(celltype = spec$celltype, metric = "Spearman rho",
             value = spec$spearman_rho, label = sprintf("%.3f", spec$spearman_rho),
             strong_both = spec$strong_both),
  data.frame(celltype = spec$celltype, metric = "Strong-gene direction agreement",
             value = spec$strong_direction_rate,
             label = percent(spec$strong_direction_rate, accuracy = 0.1),
             strong_both = spec$strong_both)
)
spec_long$metric <- factor(spec_long$metric,
                           levels = c("Spearman rho", "Strong-gene direction agreement"))
spec_long$highlight <- ifelse(spec_long$celltype == "Secretory", "Secretory", "Other compartments")

p_specificity <- ggplot(spec_long, aes(value, celltype)) +
  geom_segment(aes(x = 0, xend = value, yend = celltype, colour = highlight),
               linewidth = 1.0, alpha = 0.45) +
  geom_point(aes(size = strong_both, fill = highlight), shape = 21,
             colour = "#1D252C", stroke = 0.5) +
  geom_text(aes(label = label), hjust = -0.35, family = font_family,
            size = 3.0, colour = ink) +
  facet_grid(. ~ metric, scales = "free_x") +
  scale_colour_manual(values = c("Secretory" = accent, "Other compartments" = muted), guide = "none") +
  scale_fill_manual(values = c("Secretory" = accent, "Other compartments" = "white"), guide = "none") +
  scale_size_continuous(range = c(3.3, 5.2), name = "Strong genes") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.22))) +
  labs(x = NULL, y = NULL) +
  theme_sci(9) +
  theme(legend.position = "bottom", panel.grid.major.y = element_line(colour = grid_col),
        axis.line.x = element_line(colour = ink), axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
save_plot(p_specificity, "celltype_specificity_audit", 6.8, 3.2)

# Panel 4: complete locked-module study-effect heatmap.
hm <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                         "tables", "locked_gene_study_effect_heatmap_data.csv"))
models <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                             "tables", "locked_gene_model_estimates.csv"))
ord <- models$gene[order(models$meta_effect)]
hm$gene <- factor(hm$gene, levels = rev(ord))
hm$cohort <- factor(hm$cohort, levels = c("Paired discovery", "External cohort"))

p_heatmap <- ggplot(hm, aes(cohort, gene, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.42) +
  geom_text(aes(label = sprintf("%.1f", effect)), family = font_family,
            size = 2.35, colour = ink) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-8.5, 8.5), oob = squish,
                       name = "Effect") +
  labs(x = NULL, y = NULL) +
  theme_sci(8.3) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold", size = 8.5),
        legend.position = "right", plot.margin = margin(5, 8, 5, 5))
save_plot(p_heatmap, "locked_40gene_effect_heatmap", 4.8, 7.7)

# Panel 5: biologically representative gene-level fixed-effect estimates.
representative <- c("STATH", "DMBT1", "AZGP1", "ODAM", "GJC3", "CHRM3",
                    "CA2", "TCN1", "ANO1", "KCNN4", "SOD3", "PPP1R1B")
forest <- models[models$gene %in% representative, ]
forest <- forest[order(forest$meta_effect), ]
forest$gene <- factor(forest$gene, levels = rev(forest$gene))

p_forest <- ggplot(forest, aes(meta_effect, gene)) +
  geom_vline(xintercept = 0, colour = "#8D969E", linewidth = 0.55) +
  geom_errorbarh(aes(xmin = meta_ci_low, xmax = meta_ci_high, colour = I2_percent),
                 height = 0, linewidth = 0.85) +
  geom_point(aes(fill = I2_percent), shape = 21, size = 3.1,
             colour = "#1D252C", stroke = 0.5) +
  geom_text(aes(label = sprintf("%.2f", meta_effect)), x = 0.18,
            hjust = 0, family = font_family, size = 2.65, colour = ink) +
  scale_colour_gradient(low = blue, high = accent, limits = c(0, 70), oob = squish, guide = "none") +
  scale_fill_gradient(low = blue, high = accent, limits = c(0, 70), oob = squish,
                      name = expression(I^2~"(%)")) +
  scale_x_continuous(limits = c(-10.8, 1.55),
                     breaks = c(-10, -7.5, -5, -2.5, 0),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Cross-study fixed-effect estimate", y = NULL) +
  theme_sci(9) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.key.width = grid::unit(1.0, "cm"),
        legend.margin = margin(0, 0, 0, 0))
save_plot(p_forest, "representative_gene_meta_forest", 5.3, 4.8)

# Panel 6: leave-one-patient-out stability of the locked program.
loo <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                          "tables", "leave_one_out_summary.csv"))
loo$label <- factor(loo$label, levels = rev(loo$label))
loo$is_full <- loo$omitted_donor == "none_full"
loo$annotation <- sprintf("rho %.3f   |   %.1f%% concordant",
                          loo$rho_external, 100 * loo$strong_same_direction)

p_loo <- ggplot(loo, aes(locked_median_paired_effect, label)) +
  geom_vline(xintercept = 0, colour = "#8D969E", linewidth = 0.55) +
  geom_segment(aes(x = locked_median_paired_effect, xend = 0, yend = label),
               colour = "#AEB7BE", linewidth = 0.75) +
  geom_point(aes(fill = is_full), shape = 21, size = 4.0,
             colour = "#1D252C", stroke = 0.55) +
  geom_text(aes(x = 0.13, label = annotation), hjust = 0,
            family = font_family, size = 2.75, colour = ink) +
  scale_fill_manual(values = c(`TRUE` = accent, `FALSE` = "white"), guide = "none") +
  scale_x_continuous(limits = c(-3.2, 1.70),
                     breaks = c(-3, -2.5, -2, -1.5, -1, -0.5, 0)) +
  labs(x = "Median paired effect across the locked 40-gene module", y = NULL) +
  theme_sci(9) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(colour = grid_col),
        plot.margin = margin(8, 8, 8, 8))
save_plot(p_loo, "leave_one_patient_out_stability", 7.2, 4.0)

write.csv(effects, file.path(out_dir, "crossstudy_nondefinition_gene_effects.csv"), row.names = FALSE)
write.csv(rrho, file.path(out_dir, "rank_rank_overlap_plot_data.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
