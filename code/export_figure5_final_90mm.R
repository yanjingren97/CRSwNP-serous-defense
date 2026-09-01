options(stringsAsFactors = FALSE)
options(error = function() { traceback(2); q(status = 1) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ragg)
  library(grid)
  library(gridExtra)
})

root <- "."
out <- Sys.getenv("FIG5_FINAL_RASTER_DIR",
                  unset = file.path(root, "deliverables", "Figure5_final_90mm"))
qa_dir <- file.path(out, "QA_final_size_96dpi")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#243442"
navy <- "#173B57"
blue <- "#367FA1"
teal <- "#2A9D8F"
orange <- "#E4772A"
red <- "#BE4B42"
gold <- "#E6B655"
grey <- "#758692"
light <- "#DCE6EB"

theme_final <- function(base_size = 7) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = ink),
      axis.ticks = element_line(linewidth = 0.30, colour = ink),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_rect(fill = "#F0F5F7", colour = NA),
      strip.text = element_text(face = "bold", colour = ink),
      plot.margin = margin(3.0, 3.5, 3.0, 3.5, "mm")
    )
}

save_final <- function(p, stem, height_mm) {
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

save_two_panel_final <- function(p_top, p_bottom, stem, height_mm,
                                 heights = c(0.60, 0.40)) {
  draw_panels <- function() {
    grid.arrange(p_top, p_bottom, ncol = 1, heights = heights)
  }
  tf <- file.path(out, paste0(stem, ".tiff"))
  ragg::agg_tiff(tf, width = 90, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white")
  draw_panels()
  dev.off()
  qf <- file.path(qa_dir, paste0(stem, "_QA96.png"))
  ragg::agg_png(qf, width = 90, height = height_mm, units = "mm", res = 96,
                background = "white")
  draw_panels()
  dev.off()
}

format_figure_p <- function(x, threshold = 1e-4) {
  if (is.finite(x) && x < threshold) {
    "p < 1 × 10⁻⁴"
  } else {
    paste0("p = ", format.pval(x, digits = 2, eps = threshold))
  }
}

validation <- read.csv(file.path(root, "results", "bulk_crs_validation",
                                 "serous_module_validation.csv"))
sensitivity <- read.csv(file.path(root, "results", "bulk_crs_validation",
                                  "serous_module_sensitivity.csv"))
interaction <- read.csv(file.path(root, "results", "formal_models",
                                  "anatomy_disease_contrast_of_contrasts.csv"))
interaction$t_critical <- qt(0.975, df = interaction$df)
interaction$lower <- interaction$estimate - interaction$t_critical * interaction$se
interaction$upper <- interaction$estimate + interaction$t_critical * interaction$se
w <- 1 / interaction$se^2
pooled <- data.frame(
  dataset = "Pooled fixed effect",
  estimate = sum(w * interaction$estimate) / sum(w),
  se = sqrt(1 / sum(w))
)
pooled$lower <- pooled$estimate - 1.96 * pooled$se
pooled$upper <- pooled$estimate + 1.96 * pooled$se
pooled$p <- 2 * pnorm(-abs(pooled$estimate / pooled$se))

# A–B: density silhouettes and sample ticks for the two independent bulk cohorts.
bulk_ridge_plot <- function(dataset, group_order, group_labels, stem, height_mm,
                            add_paired_panel = FALSE) {
  d <- read.csv(file.path(root, "results", "bulk_crs_validation", paste0(dataset, "_scores.csv")))
  d <- d[d$group %in% group_order, ]
  defs <- data.frame(
    group = group_order,
    label = group_labels,
    base = rev(seq_along(group_order)),
    colour = c(navy, teal, orange),
    stringsAsFactors = FALSE
  )
  curves <- list(); ticks <- list(); meds <- list()
  for (i in seq_len(nrow(defs))) {
    q <- d$score[d$group == defs$group[i]]
    den <- density(q, n = 256, bw = "nrd0")
    amp <- den$y / max(den$y) * 0.56
    curves[[i]] <- data.frame(x = den$x, ymin = defs$base[i], ymax = defs$base[i] + amp,
                              group = defs$group[i], label = defs$label[i], base = defs$base[i])
    ticks[[i]] <- data.frame(x = q, base = defs$base[i], group = defs$group[i])
    meds[[i]] <- data.frame(x = median(q), base = defs$base[i], group = defs$group[i])
  }
  curves <- do.call(rbind, curves); ticks <- do.call(rbind, ticks); meds <- do.call(rbind, meds)
  curves$group <- factor(curves$group, levels = group_order)
  ticks$group <- factor(ticks$group, levels = group_order)
  meds$group <- factor(meds$group, levels = group_order)
  y_breaks <- defs$base
  y_labels <- paste0(defs$label, "  (n=", as.integer(table(factor(d$group, levels = group_order))), ")")

  main_contrast <- validation[validation$dataset == dataset &
                                grepl("CRSwNP_NP", validation$case) &
                                grepl("CRSwNP_(UT|IT)", validation$control), ][1, ]
  ann_prefix <- if (add_paired_panel) {
    "All samples | unpaired sensitivity"
  } else {
    "Polyp vs CRSwNP non-polyp"
  }
  ann <- paste0(ann_prefix, "\nHedges g = ", sprintf("%.2f", main_contrast$hedges_g),
                "  |  ", format_figure_p(main_contrast$mannwhitney_p))
  cols <- setNames(defs$colour, defs$group)
  xr <- range(curves$x)
  p <- ggplot(curves, aes(x)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = group), alpha = 0.82, colour = NA) +
    geom_line(aes(y = ymax, colour = group), linewidth = 0.58) +
    geom_segment(data = ticks, aes(x = x, xend = x, y = base - 0.09, yend = base + 0.10),
                 inherit.aes = FALSE, colour = ink, linewidth = 0.34, alpha = 0.80) +
    geom_segment(data = meds, aes(x = x, xend = x, y = base - 0.13, yend = base + 0.20),
                 inherit.aes = FALSE, colour = "white", linewidth = 1.15) +
    geom_segment(data = meds, aes(x = x, xend = x, y = base - 0.13, yend = base + 0.20),
                 inherit.aes = FALSE, colour = ink, linewidth = 0.48) +
    annotate("label", x = xr[2], y = max(defs$base) + 0.82, label = ann,
             hjust = 1, vjust = 1, family = font_family, size = 2.0, lineheight = 0.95,
             linewidth = 0, fill = alpha("white", 0.90), colour = ink) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_colour_manual(values = cols, guide = "none") +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       limits = c(min(defs$base) - 0.22, max(defs$base) + 0.86), expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    labs(x = "Locked Serous-defense module score", y = NULL) +
    theme_final(6.8) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 6.2), panel.grid.major.x = element_line(colour = "#E8EEF1", linewidth = 0.25))
  if (!add_paired_panel) {
    save_final(p, stem, height_mm)
    return(invisible(NULL))
  }

  paired <- read.csv(file.path(root, "results", "bulk_crs_validation",
                               "GSE136825_paired_module_scores.csv"))
  paired_stats <- read.csv(file.path(root, "results", "bulk_crs_validation",
                                     "GSE136825_paired_module_statistics.csv"))
  paired_long <- rbind(
    data.frame(pair = paste(paired$pair_series, paired$pair_id, sep = "-"),
               tissue = "CRSwNP IT", score = paired$inferior_turbinate_score,
               direction = paired$direction),
    data.frame(pair = paste(paired$pair_series, paired$pair_id, sep = "-"),
               tissue = "CRSwNP polyp", score = paired$polyp_score,
               direction = paired$direction)
  )
  paired_long$tissue <- factor(paired_long$tissue,
                               levels = c("CRSwNP IT", "CRSwNP polyp"))
  paired_subtitle <- paste0(
    paired_stats$n_lower_in_polyp, "/", paired_stats$n_pairs,
    " lower in polyp  |  median Δ = ",
    sprintf("%+.3f", paired_stats$median_paired_difference_polyp_minus_it),
    "  |  exact Wilcoxon p = ",
    formatC(paired_stats$wilcoxon_exact_two_sided_p / 1e-5,
            format = "f", digits = 3), " × 10⁻⁵"
  )
  p_pair <- ggplot(paired_long, aes(tissue, score, group = pair, colour = direction)) +
    geom_line(linewidth = 0.34, alpha = 0.62) +
    geom_point(shape = 21, fill = "white", stroke = 0.42, size = 1.35, alpha = 0.88) +
    stat_summary(aes(group = 1), fun = median, geom = "line",
                 linewidth = 1.05, colour = ink) +
    stat_summary(aes(group = 1), fun = median, geom = "point",
                 shape = 21, size = 2.45,
                 fill = gold, colour = ink, stroke = 0.55) +
    scale_colour_manual(values = c("Lower in polyp" = teal,
                                   "Higher in polyp" = orange,
                                   "No change" = grey), guide = "none") +
    labs(title = "Matched specimens | paired primary analysis",
         subtitle = paired_subtitle,
         x = NULL, y = "Module score") +
    theme_final(6.8) +
    theme(
      plot.title = element_text(size = 6.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 5.9, colour = ink, hjust = 0),
      axis.text.x = element_text(size = 6.2, face = "bold"),
      axis.text.y = element_text(size = 5.8),
      axis.title.y = element_text(size = 6.1),
      panel.grid.major.y = element_line(colour = "#E8EEF1", linewidth = 0.25),
      plot.margin = margin(2.0, 3.5, 3.0, 3.5, "mm")
    )
  save_two_panel_final(p, p_pair, stem, 106)
  invisible(NULL)
}

bulk_ridge_plot(
  "GSE36830",
  c("Healthy_UT", "CRSwNP_UT", "CRSwNP_NP"),
  c("Healthy non-polyp", "CRSwNP non-polyp", "Nasal polyp"),
  "Fig5A_GSE36830_anatomy_module_landscape_90mm", 76
)
bulk_ridge_plot(
  "GSE136825",
  c("Healthy_IT", "CRSwNP_IT", "CRSwNP_NP"),
  c("Healthy non-polyp", "CRSwNP non-polyp", "Nasal polyp"),
  "Fig5B_GSE136825_anatomy_module_landscape_90mm", 76,
  add_paired_panel = TRUE
)

# C: model-based sampling distributions for the contrast of contrasts.
ix <- rbind(
  interaction[, c("dataset", "estimate", "se", "lower", "upper", "p")],
  pooled[, c("dataset", "estimate", "se", "lower", "upper", "p")]
)
ix$dataset <- factor(ix$dataset, levels = c("GSE36830", "GSE136825", "Pooled fixed effect"))
ix$base <- c(3, 2, 1)[match(ix$dataset, levels(ix$dataset))]
ix$colour <- c(blue, teal, orange)[match(ix$dataset, levels(ix$dataset))]
xgrid <- seq(min(ix$lower) - 0.45, max(ix$upper) + 0.35, length.out = 500)
curves <- do.call(rbind, lapply(seq_len(nrow(ix)), function(i) {
  yy <- dnorm(xgrid, ix$estimate[i], ix$se[i]); yy <- yy / max(yy) * 0.64
  data.frame(x = xgrid, ymin = ix$base[i], ymax = ix$base[i] + yy,
             dataset = ix$dataset[i], base = ix$base[i])
}))
ix$summary <- paste0(sprintf("%+.2f", ix$estimate), "  [", sprintf("%+.2f", ix$lower), ", ",
                     sprintf("%+.2f", ix$upper), "]\np = ", format.pval(ix$p, digits = 2, eps = 0.001))
pC <- ggplot(curves, aes(x)) +
  geom_vline(xintercept = 0, colour = grey, linetype = "dashed", linewidth = 0.42) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = dataset), alpha = 0.82, colour = NA) +
  geom_line(aes(y = ymax, colour = dataset), linewidth = 0.58) +
  geom_segment(data = ix, aes(x = estimate, xend = estimate, y = base - 0.08, yend = base + 0.18),
               inherit.aes = FALSE, colour = ink, linewidth = 0.55) +
  geom_text(data = ix, aes(x = max(xgrid), y = base + 0.48, label = summary), inherit.aes = FALSE,
            hjust = 1, family = font_family, size = 1.95, lineheight = 0.94, colour = ink) +
  scale_fill_manual(values = c("GSE36830" = blue, "GSE136825" = teal,
                               "Pooled fixed effect" = orange), guide = "none") +
  scale_colour_manual(values = c("GSE36830" = blue, "GSE136825" = teal,
                                 "Pooled fixed effect" = orange), guide = "none") +
  scale_y_continuous(breaks = ix$base, labels = as.character(ix$dataset),
                     limits = c(0.78, 3.82), expand = c(0, 0)) +
  labs(x = "Polyp-specific decrement | contrast of contrasts", y = NULL,
       caption = "Gaussian model curves  |  brackets: 95% confidence intervals") +
  theme_final(6.8) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 6.2), plot.caption = element_text(size = 5.9, colour = grey, hjust = 0))
save_final(pC, "Fig5C_anatomy_disease_interaction_distributions_90mm", 80)

# D: six key bulk contrasts as an effect/evidence glyph matrix.
atlas <- validation[!(validation$dataset == "GSE36830" & validation$case == "CRSsNP_UT"), ]
atlas$contrast <- ifelse(grepl("CRSwNP_NP", atlas$case) & grepl("Healthy", atlas$control),
                         "Polyp vs healthy non-polyp",
                  ifelse(grepl("CRSwNP_(UT|IT)", atlas$case) & grepl("Healthy", atlas$control),
                         "CRSwNP non-polyp vs healthy",
                         "Polyp vs CRSwNP non-polyp"))
atlas$contrast <- factor(atlas$contrast,
                         levels = c("Polyp vs healthy non-polyp",
                                    "CRSwNP non-polyp vs healthy",
                                    "Polyp vs CRSwNP non-polyp"))
atlas$dataset <- factor(atlas$dataset, levels = c("GSE36830", "GSE136825"))
atlas$evidence <- -log10(atlas$mannwhitney_p)
pD <- ggplot(atlas, aes(dataset, contrast)) +
  geom_hline(yintercept = seq(0.5, 3.5, 1), colour = "#EEF2F4", linewidth = 0.4) +
  geom_vline(xintercept = seq(0.5, 2.5, 1), colour = "#EEF2F4", linewidth = 0.4) +
  geom_point(aes(size = evidence, fill = hedges_g), shape = 21, colour = ink, stroke = 0.55) +
  scale_fill_gradient2(low = navy, mid = "#F7F7F4", high = red, midpoint = 0,
                       limits = c(-2.8, 0.2), oob = squish, name = "Hedges g") +
  scale_size_continuous(range = c(4.8, 10.5), breaks = c(2, 5, 8), name = expression(-log[10](italic(p)))) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 0.82, clip = "off") +
  theme_final(6.8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold", size = 6.5), axis.text.y = element_text(size = 6.1),
        legend.position = "bottom", legend.box = "vertical", legend.key.height = unit(2.4, "mm"),
        legend.spacing.y = unit(0.8, "mm"),
        panel.border = element_rect(fill = NA, colour = light, linewidth = 0.35)) +
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top",
                                barwidth = unit(22, "mm"), barheight = unit(2.1, "mm"), order = 1),
         size = guide_legend(direction = "horizontal", title.position = "top", nrow = 1, order = 2))
save_final(pD, "Fig5D_bulk_contrast_evidence_atlas_90mm", 80)

# E: pre-specified module-definition robustness.
defs <- sensitivity[sensitivity$module %in% c("top10", "top20", "top40", "drop_top5"), ]
defs$module_label <- factor(defs$module,
                            levels = c("top10", "top20", "top40", "drop_top5"),
                            labels = c("Top 10", "Top 20", "Locked 40", "Drop top 5"))
defs$dataset <- factor(defs$dataset, levels = c("GSE36830", "GSE136825"))
pE <- ggplot(defs, aes(dataset, module_label)) +
  geom_hline(yintercept = seq(0.5, 4.5, 1), colour = "#EEF2F4", linewidth = 0.4) +
  geom_vline(xintercept = seq(0.5, 2.5, 1), colour = "#EEF2F4", linewidth = 0.4) +
  geom_point(aes(size = genes_detected, fill = hedges_g), shape = 21, colour = ink, stroke = 0.55) +
  scale_fill_gradientn(colours = c(navy, blue, "#BFDCE2"), limits = c(-2.1, -1.35),
                       oob = squish, name = "Hedges g") +
  scale_size_continuous(range = c(5.2, 12.0), breaks = c(10, 20, 40), name = "Genes measured") +
  labs(x = NULL, y = NULL) + coord_fixed(ratio = 0.72, clip = "off") +
  theme_final(6.8) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(face = "bold", size = 6.5), axis.text.y = element_text(size = 6.3),
        legend.position = "bottom", legend.box = "horizontal",
        panel.border = element_rect(fill = NA, colour = light, linewidth = 0.35)) +
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top",
                                barwidth = unit(22, "mm"), barheight = unit(2.1, "mm"), order = 1),
         size = guide_legend(direction = "horizontal", title.position = "top", nrow = 1, order = 2))
save_final(pE, "Fig5E_module_definition_robustness_90mm", 76)

# F: leave-one-gene-out effect distributions.
loo <- sensitivity[grepl("^LOO:", sensitivity$module), ]
loo$dataset <- factor(loo$dataset, levels = c("GSE36830", "GSE136825"))
full <- sensitivity[sensitivity$module == "top40", ]
full$dataset <- factor(full$dataset, levels = levels(loo$dataset))
bases <- c("GSE36830" = 2, "GSE136825" = 1)
cols <- c("GSE36830" = blue, "GSE136825" = teal)
parts <- lapply(levels(loo$dataset), function(ds) {
  q <- loo$hedges_g[loo$dataset == ds]
  den <- density(q, n = 256, bw = "nrd0")
  amp <- den$y / max(den$y) * 0.58
  data.frame(x = den$x, ymin = bases[[ds]], ymax = bases[[ds]] + amp, dataset = ds, base = bases[[ds]])
})
curves <- do.call(rbind, parts)
ticks <- transform(loo, base = unname(bases[as.character(dataset)]))
full$base <- unname(bases[as.character(full$dataset)])
ranges <- aggregate(hedges_g ~ dataset, loo, function(x) c(lo = min(x), hi = max(x)))
ranges <- data.frame(dataset = ranges$dataset, lo = ranges$hedges_g[, 1], hi = ranges$hedges_g[, 2])
ranges$base <- unname(bases[as.character(ranges$dataset)])
ranges$label <- paste0("LOO range  ", sprintf("%.2f", ranges$lo), " to ", sprintf("%.2f", ranges$hi))
pF <- ggplot(curves, aes(x)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = dataset), alpha = 0.82, colour = NA) +
  geom_line(aes(y = ymax, colour = dataset), linewidth = 0.58) +
  geom_segment(data = ticks, aes(x = hedges_g, xend = hedges_g, y = base - 0.08, yend = base + 0.08),
               inherit.aes = FALSE, colour = ink, linewidth = 0.30, alpha = 0.65) +
  geom_segment(data = full, aes(x = hedges_g, xend = hedges_g, y = base - 0.13, yend = base + 0.22),
               inherit.aes = FALSE, colour = orange, linewidth = 0.90) +
  geom_text(data = ranges, aes(x = max(curves$x), y = base + 0.43, label = label), inherit.aes = FALSE,
            hjust = 1, family = font_family, size = 2.0, colour = ink) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_y_continuous(breaks = unname(bases), labels = names(bases),
                     limits = c(0.78, 2.72), expand = c(0, 0)) +
  labs(x = "Polyp vs CRSwNP non-polyp | Hedges g", y = NULL,
       caption = "40 leave-one-gene-out modules  |  orange: full module") +
  theme_final(6.8) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 6.3),
        plot.caption = element_text(size = 5.9, colour = grey, hjust = 0))
save_final(pF, "Fig5F_leave_one_gene_out_bulk_stability_90mm", 70)

notes <- c(
  "FIGURE 5 FINAL 90-mm EXPORT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "All panels: 90 mm wide, RGB, 600 dpi, LZW TIFF, Arial.",
  "QA PNG files were rendered at 96 dpi at the same physical dimensions.",
  paste0("Pooled interaction estimate: ", sprintf("%.6f", pooled$estimate),
         "; 95% CI ", sprintf("%.6f", pooled$lower), " to ", sprintf("%.6f", pooled$upper),
         "; p=", signif(pooled$p, 5), "."),
  "Bulk validation is tissue-level orthogonal evidence and does not assign the signal exclusively to Serous cells.",
  "",
  capture.output(sessionInfo())
)
writeLines(notes, file.path(out, "R_sessionInfo_and_export_notes.txt"))
