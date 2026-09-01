options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

project_dir <- "."
final_dir <- Sys.getenv("FIG3_FINAL_RASTER_DIR", unset = tempdir())
qa_dir <- file.path(final_dir, "QA_final_size_96dpi")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#263442"
muted <- "#778692"
navy <- "#245873"
teal <- "#268E88"
orange <- "#D96C28"

save_gg_final <- function(p, stem, width_mm, height_mm) {
  agg_tiff(file.path(final_dir, paste0(stem, ".tiff")),
           width = width_mm, height = height_mm, units = "mm", res = 600,
           compression = "lzw", background = "white")
  print(p)
  dev.off()

  agg_png(file.path(qa_dir, paste0(stem, "_QA96.png")),
          width = width_mm, height = height_mm, units = "mm", res = 96,
          background = "white")
  print(p)
  dev.off()
}

save_ht_final <- function(ht, stem, width_mm, height_mm, gaps = NULL) {
  draw_once <- function() {
    draw_args <- list(
      object = ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legends = TRUE,
      padding = unit(c(4, 4, 4, 4), "mm")
    )
    if (!is.null(gaps)) draw_args$ht_gap <- gaps
    do.call(draw, draw_args)
  }

  agg_tiff(file.path(final_dir, paste0(stem, ".tiff")),
           width = width_mm, height = height_mm, units = "mm", res = 600,
           compression = "lzw", background = "white")
  draw_once()
  dev.off()

  agg_png(file.path(qa_dir, paste0(stem, "_QA96.png")),
          width = width_mm, height = height_mm, units = "mm", res = 96,
          background = "white")
  draw_once()
  dev.off()
}

# Build the two ggplot panels from their audited analysis script.
old_v4_dir <- Sys.getenv("FIG3_V4_RASTER_DIR", unset = "")
Sys.setenv(FIG3_V4_RASTER_DIR = file.path(tempdir(), "fig3_v4_source_render"))
source(file.path(project_dir, "code", "plot_figure3_editorial_v4.R"), local = FALSE)
if (nzchar(old_v4_dir)) Sys.setenv(FIG3_V4_RASTER_DIR = old_v4_dir)

p_spec_90 <- p_spec +
  labs(caption = NULL) +
  facet_wrap(~facet_label, ncol = 2) +
  scale_x_discrete(labels = c(
    "Cross-study\nrho" = "rho",
    "Direction\nagreement" = "Agree.",
    "Strong-gene\nburden" = "Genes"
  )) +
  theme(
    text = element_text(family = font_family),
    axis.text.x = element_text(family = font_family, size = 5.3, face = "bold"),
    strip.text = element_text(family = font_family, size = 7.1, face = "bold",
                              margin = margin(0, 0, 2, 0)),
    plot.margin = margin(3, 7, 3, 7)
  )

p_spec_180 <- p_spec +
  labs(caption = NULL) +
  facet_wrap(~facet_label, ncol = 4) +
  scale_x_discrete(labels = c(
    "Cross-study\nrho" = "rho",
    "Direction\nagreement" = "Agree.",
    "Strong-gene\nburden" = "Genes"
  )) +
  theme(
    text = element_text(family = font_family),
    axis.text.x = element_text(family = font_family, size = 5.8, face = "bold"),
    strip.text = element_text(family = font_family, size = 8.0, face = "bold",
                              margin = margin(0, 0, 2, 0)),
    plot.margin = margin(3, 7, 3, 7)
  )

p_density_90 <- p_density +
  labs(caption = NULL) +
  theme(
    text = element_text(family = font_family),
    axis.text.y = element_text(family = font_family, size = 5.7, face = "bold"),
    axis.text.x = element_text(family = font_family, size = 5.3),
    axis.title.x = element_text(family = font_family, size = 6.4),
    strip.text = element_text(family = font_family, size = 6.5, face = "bold",
                              margin = margin(0, 0, 3, 0)),
    panel.spacing.x = unit(2.0, "mm"),
    plot.margin = margin(3, 4, 3, 3)
  )

p_density_180 <- p_density +
  labs(caption = NULL) +
  theme(
    text = element_text(family = font_family),
    axis.text.y = element_text(family = font_family, size = 7.6, face = "bold"),
    axis.text.x = element_text(family = font_family, size = 7.0),
    axis.title.x = element_text(family = font_family, size = 8.2),
    strip.text = element_text(family = font_family, size = 8.5, face = "bold"),
    panel.spacing.x = unit(4.0, "mm"),
    plot.margin = margin(4, 5, 4, 4)
  )

save_gg_final(p_spec_90, "Fig3C_celltype_specificity_90mm", 90, 86)
save_gg_final(p_spec_180, "Fig3C_celltype_specificity_180mm", 180, 62)
save_gg_final(p_density_90, "Fig3E_representative_meta_density_90mm", 90, 82)
save_gg_final(p_density_180, "Fig3E_representative_meta_density_180mm", 180, 105)

# Rebuild the two original cross-study panels so Figure 3 contains the complete
# A-F result sequence at the same final physical sizes and typography.
old_v1_dir <- Sys.getenv("FIG3_RASTER_DIR", unset = "")
Sys.setenv(FIG3_RASTER_DIR = file.path(tempdir(), "fig3_v1_source_render"))
source(file.path(project_dir, "code", "plot_figure3_intrinsic_v1.R"), local = FALSE)
if (nzchar(old_v1_dir)) Sys.setenv(FIG3_RASTER_DIR = old_v1_dir)

p_concordance_clean <- p_concordance
# Remove the original long boxed statistics annotation. A shorter unboxed
# summary remains readable at 90 mm and avoids a report/dashboard appearance.
p_concordance_clean$layers <- p_concordance_clean$layers[-length(p_concordance_clean$layers)]
summary_label <- sprintf("15,445 genes\nrho = %.3f | concordant = %.1f%%",
                         rho, 100 * strong_agreement)
summary_x <- min(effects$effect_paired, na.rm = TRUE) + 0.35
summary_y <- max(effects$effect_external, na.rm = TRUE) - 0.20

p_A90 <- p_concordance_clean +
  annotate("text", x = summary_x, y = summary_y, label = summary_label,
           hjust = 0, vjust = 1, family = font_family, size = 2.25,
           fontface = "bold", colour = ink, lineheight = 1.05) +
  labs(x = "Paired effect (GSE235711)", y = "External effect (GSE276503)") +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(family = font_family, size = 5.8),
    axis.title = element_text(family = font_family, size = 6.4, face = "bold"),
    legend.title = element_text(family = font_family, size = 5.8),
    legend.text = element_text(family = font_family, size = 5.4),
    plot.margin = margin(4, 5, 4, 4)
  )

p_A180 <- p_concordance_clean +
  annotate("text", x = summary_x, y = summary_y, label = summary_label,
           hjust = 0, vjust = 1, family = font_family, size = 3.0,
           fontface = "bold", colour = ink, lineheight = 1.08) +
  labs(x = "Paired effect (GSE235711)", y = "External effect (GSE276503)") +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(family = font_family, size = 7.4),
    axis.title = element_text(family = font_family, size = 8.2, face = "bold"),
    legend.title = element_text(family = font_family, size = 7.2),
    legend.text = element_text(family = font_family, size = 6.8),
    plot.margin = margin(5, 6, 5, 5)
  )

p_B90 <- p_rrho +
  labs(x = "Top-ranked genes: paired cohort",
       y = "Top-ranked genes: external cohort") +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(family = font_family, size = 5.5),
    axis.title = element_text(family = font_family, size = 6.1, face = "bold"),
    strip.text = element_text(family = font_family, size = 6.4, face = "bold"),
    legend.title = element_text(family = font_family, size = 5.7),
    legend.text = element_text(family = font_family, size = 5.3),
    plot.margin = margin(4, 5, 4, 4),
    panel.spacing = unit(2.0, "mm")
  )

p_B180 <- p_rrho +
  labs(x = "Top-ranked genes included: paired discovery",
       y = "Top-ranked genes included: external cohort") +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(family = font_family, size = 7.2),
    axis.title = element_text(family = font_family, size = 8.0, face = "bold"),
    strip.text = element_text(family = font_family, size = 8.2, face = "bold"),
    legend.title = element_text(family = font_family, size = 7.1),
    legend.text = element_text(family = font_family, size = 6.7),
    plot.margin = margin(5, 6, 5, 5),
    panel.spacing = unit(3.0, "mm")
  )

save_gg_final(p_A90, "Fig3A_crossstudy_concordance_90mm", 90, 78)
save_gg_final(p_A180, "Fig3A_crossstudy_concordance_180mm", 180, 142)
save_gg_final(p_B90, "Fig3B_RRHO_landscape_90mm", 90, 62)
save_gg_final(p_B180, "Fig3B_RRHO_landscape_180mm", 180, 90)

# 40-gene multi-evidence heatmap, rebuilt for each final physical width.
models <- read.csv(file.path(project_dir, "deliverables", "stage_03_full_figure_set",
                             "tables", "locked_gene_model_estimates.csv"))
ord <- order(models$meta_effect)
models <- models[ord, ]
key_genes <- c("STATH", "DMBT1", "AZGP1", "ODAM", "GJC3", "PPP1R1B",
               "CHRM3", "CA2", "TCN1", "ANO1")

effect_col <- colorRamp2(c(-8.5, -6, -4, -2, 0),
                         c("#123B55", "#1F6783", "#5D9BB1", "#C4DCE4", "#FAFBFC"))
i2_col <- colorRamp2(c(0, 35, 70), c("#F2F4F5", "#E7B98F", "#C8642E"))
q_score <- pmin(-log10(pmax(models$meta_fdr, .Machine$double.xmin)), 50)
q_col <- colorRamp2(c(0, 10, 30, 50), c("#F4F6F6", "#BBD9D5", "#59A9A1", "#187C79"))

make_40gene_final <- function(mode = c("90", "180")) {
  mode <- match.arg(mode)
  small <- mode == "90"
  rn_size <- if (small) 6.2 else 7.3
  cn_size <- if (small) 6.2 else 8.0
  lg_size <- if (small) 5.9 else 7.2
  strip_w <- if (small) unit(4.0, "mm") else unit(6.0, "mm")

  main_mat <- cbind(`Paired` = models$effect_paired,
                    `External` = models$effect_external)
  rownames(main_mat) <- models$gene
  pooled_mat <- matrix(models$meta_effect, ncol = 1,
                       dimnames = list(models$gene, "Meta"))
  i2_mat <- matrix(models$I2_percent, ncol = 1,
                   dimnames = list(models$gene, "I2"))
  q_mat <- matrix(q_score, ncol = 1,
                  dimnames = list(models$gene, "FDR"))

  rn_col <- ifelse(models$gene %in% key_genes, "#C95D25", "#2A3742")
  rn_face <- ifelse(models$gene %in% key_genes, "bold", "plain")

  h1 <- Heatmap(
    main_mat, name = "Cohort effect", col = effect_col,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.75), border = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontfamily = font_family, fontsize = rn_size,
                        col = rn_col, fontface = rn_face),
    row_names_max_width = unit(if (small) 19 else 27, "mm"),
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    heatmap_legend_param = list(
      title = "Polyp - non-polyp effect", direction = "horizontal",
      title_position = "topleft", legend_width = unit(if (small) 27 else 42, "mm"),
      at = c(-8, -4, 0),
      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                      fontface = "bold")
    ), use_raster = FALSE
  )
  h2 <- Heatmap(
    pooled_mat, name = "Pooled effect", col = effect_col,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.75), border = FALSE,
    width = strip_w,
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    show_heatmap_legend = FALSE, use_raster = FALSE
  )
  h3 <- Heatmap(
    i2_mat, name = "I2 (%)", col = i2_col,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.75), border = FALSE,
    width = strip_w,
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    heatmap_legend_param = list(
      title = "I2 (%)", direction = "horizontal", title_position = "topleft",
      legend_width = unit(if (small) 18 else 27, "mm"), at = c(0, 40, 80),
      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                      fontface = "bold")
    ), use_raster = FALSE
  )
  h4 <- Heatmap(
    q_mat, name = "-log10(q)", col = q_col,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.75), border = FALSE,
    width = strip_w,
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    heatmap_legend_param = list(
      title = "-log10(q)", direction = "horizontal", title_position = "topleft",
      legend_width = unit(if (small) 19 else 28, "mm"), at = c(0, 25, 50),
      labels = c("0", "25", "50+"),
      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                      fontface = "bold")
    ), use_raster = FALSE
  )
  h1 + h2 + h3 + h4
}

ht40_90 <- make_40gene_final("90")
ht40_180 <- make_40gene_final("180")
save_ht_final(ht40_90, "Fig3D_locked40_multievidence_90mm", 90, 160,
              gaps = unit(c(2.0, 2.0, 1.6), "mm"))
save_ht_final(ht40_180, "Fig3D_locked40_multievidence_180mm", 180, 178,
              gaps = unit(c(4.0, 3.0, 3.0), "mm"))

# Leave-one-patient-out gene landscape, compact at 90 mm and expanded at 180 mm.
loo <- read.csv(file.path(project_dir, "results", "robustness",
                          "serous_patient_leave_one_out_genes.csv"))
scenario_order <- c("none_full", "CRSwNP_1", "CRSwNP_2", "CRSwNP_3", "CRSwNP_6")
scenario_labels <- c("Full", "-P1", "-P2", "-P3", "-P6")
gene_full <- loo[loo$omitted_donor == "none_full", c("gene", "paired", "external")]
genes_loo <- gene_full$gene[order(gene_full$paired)]
paired_mat <- sapply(scenario_order, function(s) {
  z <- loo[loo$omitted_donor == s, c("gene", "paired")]
  setNames(z$paired, z$gene)[genes_loo]
})
colnames(paired_mat) <- scenario_labels
rownames(paired_mat) <- genes_loo
delta_mat <- sweep(paired_mat[, -1, drop = FALSE], 1, paired_mat[, 1], "-")
colnames(delta_mat) <- scenario_labels[-1]
ext_vec <- setNames(gene_full$external, gene_full$gene)[genes_loo]
ext_mat <- matrix(ext_vec, ncol = 1, dimnames = list(genes_loo, "Ext."))
max_change <- apply(abs(delta_mat), 1, max)
max_mat <- matrix(max_change, ncol = 1, dimnames = list(genes_loo, "Max"))
delta_lim <- max(2, ceiling(max(abs(delta_mat)) * 2) / 2)
delta_col <- colorRamp2(c(-delta_lim, 0, delta_lim),
                        c("#167D78", "#F7F8F8", "#D96C28"))
max_col <- colorRamp2(c(0, 1, 2), c("#F2F5F6", "#92B2BF", "#334F63"))

make_loo_final <- function(mode = c("90", "180")) {
  mode <- match.arg(mode)
  small <- mode == "90"
  rn_size <- if (small) 6.0 else 7.0
  cn_size <- if (small) 6.2 else 7.5
  title_size <- if (small) 6.4 else 8.5
  lg_size <- if (small) 5.9 else 7.0
  strip_w <- if (small) unit(3.6, "mm") else unit(6.0, "mm")
  rn_col <- ifelse(genes_loo %in% key_genes, "#C95D25", "#2A3742")
  rn_face <- ifelse(genes_loo %in% key_genes, "bold", "plain")

  aclass <- factor(c("Full", rep("LOPO", 4)), levels = c("Full", "LOPO"))
  tha <- HeatmapAnnotation(
    Analysis = aclass,
    col = list(Analysis = c(Full = orange, LOPO = teal)),
    simple_anno_size = unit(if (small) 2.4 else 3.3, "mm"),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Analysis = list(direction = "horizontal", nrow = 1,
                      title_position = "topleft",
                      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
                      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                                      fontface = "bold"))
    )
  )
  h1 <- Heatmap(
    paired_mat, name = "Paired effect", col = effect_col,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.72), border = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontfamily = font_family, fontsize = rn_size,
                        col = rn_col, fontface = rn_face),
    row_names_max_width = unit(if (small) 18 else 27, "mm"),
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    column_title = "Paired effects",
    column_title_gp = gpar(fontfamily = font_family, fontsize = title_size,
                           fontface = "bold", col = ink),
    top_annotation = tha,
    heatmap_legend_param = list(
      title = "Polyp - non-polyp effect", direction = "horizontal",
      title_position = "topleft", legend_width = unit(if (small) 26 else 40, "mm"),
      at = c(-8, -4, 0),
      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                      fontface = "bold")
    ), use_raster = FALSE
  )
  h2 <- Heatmap(
    delta_mat, name = "Omission change", col = delta_col,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.72), border = FALSE,
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    column_title = "Change from full",
    column_title_gp = gpar(fontfamily = font_family, fontsize = title_size,
                           fontface = "bold", col = ink),
    heatmap_legend_param = list(
      title = "Effect change", direction = "horizontal",
      title_position = "topleft", legend_width = unit(if (small) 22 else 34, "mm"),
      at = c(-delta_lim, 0, delta_lim),
      labels = c(paste0("-", delta_lim), "0", paste0("+", delta_lim)),
      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                      fontface = "bold")
    ), use_raster = FALSE
  )
  h3 <- Heatmap(
    ext_mat, name = "External effect", col = effect_col,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.72), border = FALSE,
    width = strip_w, show_column_names = TRUE,
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    show_heatmap_legend = FALSE, use_raster = FALSE
  )
  h4 <- Heatmap(
    max_mat, name = "Max |change|", col = max_col,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.72), border = FALSE,
    width = strip_w, show_column_names = TRUE,
    column_names_side = "top", column_names_rot = 0,
    column_names_gp = gpar(fontfamily = font_family, fontsize = cn_size,
                           fontface = "bold", col = ink),
    heatmap_legend_param = list(
      title = "Max |change|", direction = "horizontal",
      title_position = "topleft", legend_width = unit(if (small) 18 else 27, "mm"),
      at = c(0, 1, 2),
      labels_gp = gpar(fontfamily = font_family, fontsize = lg_size),
      title_gp = gpar(fontfamily = font_family, fontsize = lg_size + 0.4,
                      fontface = "bold")
    ), use_raster = FALSE
  )
  h1 + h2 + h3 + h4
}

htloo_90 <- make_loo_final("90")
htloo_180 <- make_loo_final("180")
save_ht_final(htloo_90, "Fig3F_LOPO_gene_landscape_90mm", 90, 160,
              gaps = unit(c(1.8, 1.5, 1.5), "mm"))
save_ht_final(htloo_180, "Fig3F_LOPO_gene_landscape_180mm", 180, 170,
              gaps = unit(c(3.0, 2.0, 2.0), "mm"))

writeLines(c(
  "Figure 3 final TIFF export session",
  "All TIFF files: RGB, 600 dpi, LZW compression",
  "Font: Arial",
  paste0("Generated: ", Sys.time()),
  capture.output(sessionInfo())
), file.path(final_dir, "R_sessionInfo_and_export_notes.txt"))
