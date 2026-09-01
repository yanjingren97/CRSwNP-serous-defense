options(stringsAsFactors = FALSE)
options(error = function() { traceback(2); q(status = 1) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ragg)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})

root <- Sys.getenv("ENT_ROOT", unset = ".")
out <- Sys.getenv("FIG4_FINAL_RASTER_DIR",
                  unset = file.path(root, "deliverables", "Figure4_final_90mm"))
qa_dir <- file.path(out, "QA_final_size_96dpi")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
navy <- "#183B56"
blue <- "#2B6F92"
teal <- "#2A9D8F"
cyan <- "#65B8C4"
orange <- "#E4772A"
gold <- "#E6B655"
red <- "#BC493F"
ink <- "#243442"
grey <- "#768692"
light <- "#DCE6EB"

theme_final <- function(base_size = 7) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = ink),
      axis.ticks = element_line(linewidth = 0.30, colour = ink),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.6),
      strip.background = element_rect(fill = "#F1F5F7", colour = NA),
      strip.text = element_text(face = "bold", colour = ink),
      plot.margin = margin(3.0, 3.5, 3.0, 3.5, "mm")
    )
}

save_gg_final <- function(p, stem, height_mm) {
  tf <- file.path(out, paste0(stem, ".tiff"))
  ragg::agg_tiff(tf, width = 90, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white")
  print(p)
  dev.off()
  qf <- file.path(qa_dir, paste0(stem, "_QA96.png"))
  ragg::agg_png(qf, width = 90, height = height_mm, units = "mm", res = 96,
                background = "white")
  print(p)
  dev.off()
}

save_ht_final <- function(draw_fun, stem, height_mm) {
  tf <- file.path(out, paste0(stem, ".tiff"))
  ragg::agg_tiff(tf, width = 90, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white")
  draw_fun()
  dev.off()
  qf <- file.path(qa_dir, paste0(stem, "_QA96.png"))
  ragg::agg_png(qf, width = 90, height = height_mm, units = "mm", res = 96,
                background = "white")
  draw_fun()
  dev.off()
}

wrap_term <- function(x, width = 18) {
  x <- gsub(" \\(GO:[0-9]+\\)$", "", x)
  vapply(x, function(z) paste(utils::head(strwrap(z, width = width), 2), collapse = "\n"), character(1))
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) return(0)
  length(intersect(a, b)) / length(u)
}

select_nonredundant_terms <- function(d, direction, n = 17) {
  q <- d[d$direction == direction & is.finite(d$fdr) & d$fdr < 0.05 &
           d$overlap >= 4 & d$term_size <= 700, , drop = FALSE]
  q$score <- -log10(pmax(q$fdr, 1e-50))
  q <- q[order(-q$score, -q$overlap, q$term_size), , drop = FALSE]
  q <- utils::head(q, 140)
  gs <- strsplit(as.character(q$genes), ";", fixed = TRUE)
  keep <- integer(0)
  for (i in seq_len(nrow(q))) {
    if (!length(keep) || all(vapply(keep, function(k) jaccard(gs[[i]], gs[[k]]) < 0.58,
                                     logical(1)))) {
      keep <- c(keep, i)
    }
    if (length(keep) >= n) break
  }
  q <- q[keep, , drop = FALSE]
  q$gene_set <- I(gs[keep])
  q
}

make_term_map <- function(q, direction = c("down", "up")) {
  direction <- match.arg(direction)
  n <- nrow(q)
  jm <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) jm[i, j] <- jaccard(q$gene_set[[i]], q$gene_set[[j]])
  diag(jm) <- 1
  set.seed(ifelse(direction == "down", 114, 221))
  xy <- tryCatch(cmdscale(as.dist(1 - jm), k = 2, add = TRUE), error = function(e) NULL)
  if (is.null(xy) || ncol(as.matrix(xy)) < 2) {
    a <- seq(0, 2 * pi, length.out = n + 1)[-1]
    xy <- cbind(cos(a), sin(a))
  }
  xy <- as.matrix(xy)
  xy[, 1] <- rescale(xy[, 1], c(-0.82, 0.82))
  xy[, 2] <- rescale(xy[, 2], c(-0.82, 0.82))
  nodes <- transform(q, id = seq_len(n), x = xy[, 1], y = xy[, 2])

  ep <- which(upper.tri(jm) & jm > 0.04, arr.ind = TRUE)
  edges <- if (nrow(ep)) data.frame(
    x = nodes$x[ep[, 1]], y = nodes$y[ep[, 1]],
    xend = nodes$x[ep[, 2]], yend = nodes$y[ep[, 2]],
    overlap = jm[ep]
  ) else data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric(), overlap = numeric())
  if (nrow(edges) > 34) edges <- edges[order(-edges$overlap), ][seq_len(34), ]

  lab <- nodes[order(-nodes$score, -nodes$overlap), ][seq_len(min(8, n)), ]
  lab$side <- ifelse(lab$x >= median(nodes$x), "right", "left")
  lab$lx <- ifelse(lab$side == "right", 1.00, -1.00)
  lab$ly <- NA_real_
  for (s in c("left", "right")) {
    ii <- which(lab$side == s)
    ii <- ii[order(-lab$y[ii])]
    lab$ly[ii] <- seq(0.82, -0.82, length.out = length(ii))
  }
  lab$hjust <- ifelse(lab$side == "right", 0, 1)
  lab$label <- wrap_term(lab$term)

  pal <- if (direction == "down") c("#DDEEF1", cyan, teal, navy) else c("#F8E7C4", gold, orange, red)
  strip_label <- if (direction == "down") "Concordant decreases" else "Concordant increases"

  ggplot() +
    geom_segment(data = edges,
                 aes(x, y, xend = xend, yend = yend, alpha = overlap),
                 colour = "#9CB0BC", linewidth = 0.45, lineend = "round") +
    geom_segment(data = lab, aes(x, y, xend = lx, yend = ly),
                 colour = "#A8B8C1", linewidth = 0.30) +
    geom_point(data = nodes,
               aes(x, y, size = overlap, fill = score, shape = library),
               colour = ink, stroke = 0.35) +
    geom_text(data = lab, aes(lx, ly, label = label, hjust = hjust),
              family = font_family, size = 2.05, lineheight = 0.90, colour = ink) +
    annotate("label", x = -1.55, y = 1.06, label = strip_label,
             hjust = 0, family = font_family, fontface = "bold", size = 2.25,
             linewidth = 0, fill = if (direction == "down") "#E6F2F2" else "#FAECD9",
             colour = if (direction == "down") navy else "#8A3D1F") +
    scale_fill_gradientn(colours = pal, name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(2.2, 7.2), name = "Genes") +
    scale_shape_manual(values = c("GO_BP_2025" = 21, "Reactome_2024" = 24),
                       labels = c("GO biological process", "Reactome"), name = NULL) +
    scale_alpha_continuous(range = c(0.18, 0.65), guide = "none") +
    coord_equal(xlim = c(-1.68, 1.68), ylim = c(-1.02, 1.10), clip = "off") +
    theme_void(base_family = font_family, base_size = 6.2) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-2, 0, 0, 0),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(7.0, "mm"),
      legend.text = element_text(size = 6.0, colour = ink),
      legend.title = element_text(size = 6.1, colour = ink),
      plot.margin = margin(3, 5, 1.5, 5, "mm")
    ) +
    guides(
      fill = guide_colourbar(order = 1, direction = "horizontal", title.position = "top",
                             barwidth = unit(19, "mm"), barheight = unit(2.0, "mm")),
      size = guide_legend(order = 2, direction = "horizontal", title.position = "top", nrow = 1),
      shape = guide_legend(order = 3, direction = "horizontal", nrow = 1,
                           override.aes = list(size = 3, fill = "white"))
    )
}

# A–B: functional architecture from concordant genes.
enr <- read.csv(file.path(root, "results", "formal_models",
                          "serous_concordant_functional_enrichment.csv"), check.names = FALSE)
down_terms <- select_nonredundant_terms(enr, "down", 17)
up_terms <- select_nonredundant_terms(enr, "up", 17)
save_gg_final(make_term_map(down_terms, "down"), "Fig4A_downregulated_functional_map_90mm", 92)
save_gg_final(make_term_map(up_terms, "up"), "Fig4B_upregulated_functional_map_90mm", 92)

# C: descriptive Serous diffusion manifold coloured by the audited state coordinate.
st <- read.csv(file.path(root, "results", "advanced_singlecell", "serous_cell_state_axis.csv"))
dc <- read.csv(file.path(root, "results", "advanced_singlecell", "serous_diffusion_coordinates.csv"))
set.seed(20260831)
dc_plot <- if (nrow(dc) > 18000) dc[sample(seq_len(nrow(dc)), 18000), ] else dc
ctr <- c(mean(dc_plot$DC1), mean(dc_plot$DC2))
mat_hi <- dc_plot[dc_plot$mature_score >= quantile(dc_plot$mature_score, 0.94), ]
cil_hi <- dc_plot[dc_plot$cilia_score >= quantile(dc_plot$cilia_score, 0.94), ]
ends <- data.frame(
  x = c(mean(mat_hi$DC1), mean(cil_hi$DC1)),
  y = c(mean(mat_hi$DC2), mean(cil_hi$DC2)),
  label = c("Maturation-rich", "Cilia-rich"),
  colour = c(orange, teal)
)
pC <- ggplot(dc_plot, aes(DC1, DC2)) +
  geom_point(aes(colour = failure_axis), size = 0.52, alpha = 0.58, stroke = 0) +
  stat_density_2d(colour = alpha("white", 0.60), linewidth = 0.28, bins = 5) +
  scale_colour_gradientn(colours = c(navy, "#7FB8C1", "#F2E7C9", orange),
                         limits = c(0, 1), name = "State\ncoordinate") +
  geom_segment(data = ends,
               aes(x = ctr[1], y = ctr[2], xend = x, yend = y),
               inherit.aes = FALSE, linewidth = 0.72,
               colour = c(orange, teal),
               arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  geom_text(data = ends, aes(x, y, label = label), inherit.aes = FALSE,
            family = font_family, fontface = "bold", size = 2.05, vjust = -0.75,
            colour = c("#9A451E", "#176B62")) +
  labs(x = "Diffusion component 1", y = "Diffusion component 2") +
  coord_equal(xlim = quantile(dc_plot$DC1, c(0.003, 0.997)),
              ylim = quantile(dc_plot$DC2, c(0.003, 0.997))) + theme_final(7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = unit(24, "mm"), legend.key.height = unit(2.2, "mm"),
        panel.border = element_rect(fill = NA, colour = light, linewidth = 0.35)) +
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5))
save_gg_final(pC, "Fig4C_serous_state_continuum_90mm", 78)

# D: cohort-wise ridge audit; cell densities are descriptive and vertical ticks mark analysis-unit medians.
sd <- read.csv(file.path(root, "results", "advanced_singlecell", "serous_donor_state_axis.csv"))
defs <- data.frame(
  cohort = c("Paired discovery | GSE235711", "Paired discovery | GSE235711",
             "External cohort | GSE276503", "External cohort | GSE276503"),
  dataset = c("GSE235711", "GSE235711", "GSE276503", "GSE276503"),
  tissue = c("Ethmoid", "Nasal_polyp", "Inferior_turbinate", "Nasal_polyp"),
  label = c("Ethmoid", "Nasal polyp", "Inferior turbinate", "Nasal polyp"),
  base = c(1, 2, 1, 2),
  stringsAsFactors = FALSE
)
ridge_parts <- list()
tick_parts <- list()
for (i in seq_len(nrow(defs))) {
  z <- defs[i, ]
  if (z$dataset == "GSE235711") {
    vals <- st$failure_axis[st$dataset == z$dataset & st$disease == "CRSwNP" & st$tissue == z$tissue]
    ds <- sd[sd$dataset == z$dataset & sd$disease == "CRSwNP" & sd$tissue == z$tissue, ]
  } else {
    vals <- st$failure_axis[st$dataset == z$dataset & st$disease == "CRSwNP" & st$tissue == z$tissue]
    ds <- sd[sd$dataset == z$dataset & sd$disease == "CRSwNP" & sd$tissue == z$tissue, ]
  }
  den <- density(vals, from = 0, to = 1, n = 256, bw = 0.045)
  amp <- den$y / max(den$y) * 0.62
  ridge_parts[[i]] <- data.frame(cohort = z$cohort, tissue = z$tissue, label = z$label,
                                 base = z$base, x = den$x, ymin = z$base, ymax = z$base + amp)
  tick_parts[[i]] <- data.frame(cohort = z$cohort, tissue = z$tissue, label = z$label,
                                base = z$base, x = ds$median_failure_axis)
}
ridges <- do.call(rbind, ridge_parts)
ticks <- do.call(rbind, tick_parts)
cohort_order <- c("Paired discovery | GSE235711", "External cohort | GSE276503")
ridges$cohort <- factor(ridges$cohort, levels = cohort_order)
ticks$cohort <- factor(ticks$cohort, levels = cohort_order)
ann <- data.frame(
  cohort = c("Paired discovery | GSE235711", "External cohort | GSE276503"),
  x = c(0.98, 0.98), y = c(2.70, 2.70),
  label = c("Patient-level Δ = +0.107  |  P = 0.625",
            "Biopsy-level Δ = +0.112  |  P = 0.850")
)
ann$cohort <- factor(ann$cohort, levels = cohort_order)
ridge_labs <- unique(ridges[, c("cohort", "label", "base")])
tissue_cols <- c("Ethmoid" = blue, "Inferior_turbinate" = teal, "Nasal_polyp" = orange)
pD <- ggplot(ridges, aes(x = x)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = tissue), alpha = 0.82, colour = NA) +
  geom_line(aes(y = ymax, colour = tissue), linewidth = 0.55) +
  geom_segment(data = ticks, aes(x = x, xend = x, y = base - 0.08, yend = base + 0.10),
               inherit.aes = FALSE, colour = ink, linewidth = 0.38) +
  geom_text(data = ridge_labs, aes(x = 0.015, y = base + 0.12, label = label),
            inherit.aes = FALSE, hjust = 0, family = font_family, fontface = "bold",
            size = 2.05, colour = ink) +
  geom_text(data = ann, aes(x, y, label = label), inherit.aes = FALSE,
            hjust = 1, family = font_family, size = 2.00, colour = ink) +
  facet_wrap(~cohort, ncol = 1) +
  scale_fill_manual(values = tissue_cols, guide = "none") +
  scale_colour_manual(values = tissue_cols, guide = "none") +
  scale_y_continuous(breaks = NULL,
                     limits = c(0.78, 2.78), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Descriptive Serous state coordinate", y = NULL,
       caption = "Filled curves: cell distributions  |  vertical ticks: analysis-unit medians") +
  theme_final(6.7) +
  theme(panel.spacing = unit(3.2, "mm"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 6.5), plot.caption = element_text(size = 6.0, colour = grey, hjust = 0))
save_gg_final(pD, "Fig4D_state_coordinate_donor_audit_90mm", 92)

# E: all-activity cross-cohort screening landscape.
rg <- read.csv(file.path(root, "results", "advanced_regulatory", "regulatory_cross_cohort_summary.csv"))
rg <- rg[is.finite(rg$effect_paired) & is.finite(rg$effect_external), ]
rank_pct <- function(x) (rank(x, ties.method = "average") - 0.5) / length(x)
rg$paired_percentile <- ave(rg$effect_paired, rg$network, FUN = rank_pct)
rg$external_percentile <- ave(rg$effect_external, rg$network, FUN = rank_pct)
rg$priority <- tolower(as.character(rg$replicated_priority)) == "true"
priority <- rg[rg$priority, ]
priority <- priority[order(-priority$meta_effect, priority$source), ]
priority$key <- seq_len(nrow(priority))
priority$lx <- NA_real_; priority$ly <- NA_real_
top_ids <- seq_len(min(5, nrow(priority)))
right_ids <- setdiff(seq_len(nrow(priority)), top_ids)
priority$lx[top_ids] <- seq(0.54, 0.94, length.out = length(top_ids))
priority$ly[top_ids] <- 0.985
if (length(right_ids)) {
  priority$lx[right_ids] <- 0.985
  priority$ly[right_ids] <- seq(0.84, 0.56, length.out = length(right_ids))
}
rg_bg <- rg[!rg$priority, ]
pE <- ggplot(rg_bg, aes(paired_percentile, external_percentile)) +
  stat_density_2d(colour = "#C2D0D8", linewidth = 0.35, bins = 5) +
  geom_point(aes(shape = network), size = 1.25, alpha = 0.38, colour = "#617784") +
  geom_hline(yintercept = 0.5, colour = light, linewidth = 0.35) +
  geom_vline(xintercept = 0.5, colour = light, linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9CAFBA", linewidth = 0.40) +
  geom_point(data = priority, aes(paired_percentile, external_percentile), inherit.aes = FALSE,
             shape = 21, size = 2.8, fill = orange, colour = ink, stroke = 0.50) +
  geom_segment(data = priority,
               aes(x = paired_percentile, y = external_percentile, xend = lx, yend = ly),
               inherit.aes = FALSE, colour = "#9B6A4E", linewidth = 0.30) +
  geom_point(data = priority, aes(lx, ly), inherit.aes = FALSE,
             shape = 21, size = 3.3, fill = orange, colour = ink, stroke = 0.50) +
  geom_text(data = priority, aes(lx, ly, label = key), inherit.aes = FALSE,
            family = font_family, fontface = "bold", colour = "white", size = 2.05) +
  annotate("label", x = 0.03, y = 0.97,
           label = paste0(nrow(rg), " activities screened\n", nrow(priority), " replicated priorities"),
           hjust = 0, vjust = 1, family = font_family, size = 2.0, lineheight = 0.95,
           linewidth = 0, fill = alpha("white", 0.88), colour = ink) +
  scale_shape_manual(values = c("CollecTRI" = 16, "PROGENy" = 17), name = NULL) +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1), expand = c(0.015, 0.015)) +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1), expand = c(0.015, 0.015)) +
  labs(x = "Activity-effect percentile | paired cohort",
       y = "Activity-effect percentile | external cohort") +
  coord_equal() + theme_final(6.7) +
  theme(legend.position = c(0.80, 0.14), legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
        panel.border = element_rect(fill = NA, colour = light, linewidth = 0.35))
save_gg_final(pE, "Fig4E_regulatory_screening_landscape_90mm", 84)

# F: row-standardized priority activity fingerprint across anatomical groups.
act <- rbind(
  read.csv(file.path(root, "results", "advanced_regulatory", "collectri_sample_activities.csv")),
  read.csv(file.path(root, "results", "advanced_regulatory", "progeny_sample_activities.csv"))
)
act <- act[act$source %in% priority$source, ]
act$group_id <- paste(act$dataset, act$disease, act$tissue, sep = "|")
group_defs <- data.frame(
  group_id = c(
    "GSE235711|Healthy|Ethmoid", "GSE235711|CRSsNP|Ethmoid",
    "GSE235711|CRSwNP|Ethmoid", "GSE235711|CRSwNP|Nasal_polyp",
    "GSE276503|Healthy|Inferior_turbinate", "GSE276503|CRSwNP|Inferior_turbinate",
    "GSE276503|CRSwNP|Middle_turbinate", "GSE276503|CRSwNP|Nasal_polyp"
  ),
  label = c("Healthy ET", "CRSsNP ET", "CRSwNP ET", "CRSwNP polyp",
            "Healthy IT", "CRSwNP IT", "CRSwNP MT", "CRSwNP polyp"),
  dataset = rep(c("GSE235711", "GSE276503"), each = 4),
  stringsAsFactors = FALSE
)
ag <- aggregate(activity ~ source + group_id, act, mean)
mat <- matrix(NA_real_, nrow = nrow(priority), ncol = nrow(group_defs),
              dimnames = list(priority$source, group_defs$label))
for (i in seq_len(nrow(ag))) {
  rr <- match(ag$source[i], rownames(mat)); cc <- match(ag$group_id[i], group_defs$group_id)
  if (!is.na(rr) && !is.na(cc)) mat[rr, cc] <- ag$activity[i]
}
mat_z <- t(apply(mat, 1, function(x) {
  z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  pmax(-2.5, pmin(2.5, z))
}))
rownames(mat_z) <- paste(priority$key, priority$source, sep = "  ")
col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5), c(blue, "#F7F7F4", red))
col_split <- factor(group_defs$dataset, levels = c("GSE235711", "GSE276503"))
htF <- Heatmap(
  mat_z, name = "Row z-score", col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split, cluster_column_slices = FALSE,
  column_gap = unit(2.1, "mm"),
  row_names_gp = gpar(fontfamily = font_family, fontsize = 6.2, col = ink),
  column_names_gp = gpar(fontfamily = font_family, fontsize = 6.1, col = ink),
  column_labels = group_defs$label,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_centered = TRUE,
  column_names_max_height = unit(12, "mm"),
  column_names_rot = 42,
  row_names_side = "left",
  rect_gp = gpar(col = "white", lwd = 0.75),
  heatmap_legend_param = list(
    direction = "horizontal", title_position = "topcenter",
    at = c(-2, 0, 2), labels = c("-2", "0", "+2"),
    legend_width = unit(29, "mm"),
    title_gp = gpar(fontfamily = font_family, fontsize = 6.2),
    labels_gp = gpar(fontfamily = font_family, fontsize = 5.7)
  ),
  column_title_gp = gpar(fontfamily = font_family, fontsize = 6.5, fontface = "bold", col = ink)
)
draw_F <- function() {
  grid.newpage()
  draw(htF, heatmap_legend_side = "bottom", merge_legends = TRUE,
       padding = unit(c(3, 3, 2, 3), "mm"))
}
save_ht_final(draw_F, "Fig4F_priority_activity_fingerprint_90mm", 92)

notes <- c(
  "FIGURE 4 FINAL 90-mm EXPORT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "All panels: 90 mm wide, RGB, 600 dpi, LZW TIFF, Arial.",
  "QA PNG files were rendered at 96 dpi at the same physical dimensions.",
  "Panel D cell densities are descriptive; annotations use cohort-specific analysis units (paired patients or biopsies).",
  "The state coordinate is not pseudotime and does not establish lineage conversion.",
  "Regulatory activities are target-response inferences and do not establish causal TF activation.",
  "",
  capture.output(sessionInfo())
)
writeLines(notes, file.path(out, "R_sessionInfo_and_export_notes.txt"))
