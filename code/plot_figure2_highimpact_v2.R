options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ragg)
  library(grid)
})

root <- Sys.getenv("ENT_ROOT", unset = ".")
out <- file.path(root, "deliverables/figure2_highimpact_v2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
skip_exports <- identical(Sys.getenv("ENT_SKIP_PANEL_EXPORT", unset = "0"), "1")

font <- "Arial"
ink <- "#17212B"
muted <- "#667786"
grid_col <- "#DCE3E8"
paper <- "#FFFFFF"
blue <- "#0072B2"
sky <- "#56B4E9"
green <- "#009E73"
orange <- "#E69F00"
vermillion <- "#D55E00"
purple <- "#7656A8"
gold <- "#C7A76C"
grey <- "#8C98A4"

state_cols <- c(
  "Serous_glandular" = blue,
  "Mucous_glandular" = orange,
  "Surface_club" = green,
  "Goblet" = gold,
  "Inflammatory_secretory" = vermillion,
  "Transitional_secretory" = purple,
  "Ambiguous_secretory" = grey
)
state_order <- names(state_cols)
state_labels <- c(
  "Serous_glandular" = "Serous glandular",
  "Mucous_glandular" = "Mucous glandular",
  "Surface_club" = "Surface club",
  "Goblet" = "Goblet",
  "Inflammatory_secretory" = "Inflammatory secretory",
  "Transitional_secretory" = "Transitional secretory",
  "Ambiguous_secretory" = "Ambiguous secretory"
)

theme_hi <- function(base_size = 9) {
  theme_classic(base_family = font, base_size = base_size) +
    theme(
      plot.background = element_rect(fill = paper, colour = NA),
      panel.background = element_rect(fill = paper, colour = NA),
      axis.line = element_line(colour = ink, linewidth = 0.45),
      axis.ticks = element_line(colour = ink, linewidth = 0.40),
      axis.ticks.length = unit(1.6, "mm"),
      axis.text = element_text(colour = ink),
      axis.title = element_text(colour = ink),
      legend.text = element_text(colour = ink),
      legend.title = element_text(colour = ink),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#F2F5F7", colour = NA),
      strip.text = element_text(face = "bold", colour = ink),
      plot.margin = margin(6, 7, 6, 7, "mm")
    )
}

save_panel <- function(plot, stem, width_mm, height_mm) {
  if (skip_exports) return(invisible(NULL))
  ggsave(file.path(out, paste0(stem, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm",
         device = cairo_pdf, family = font, bg = paper)
  ggsave(file.path(out, paste0(stem, ".svg")), plot,
         width = width_mm, height = height_mm, units = "mm",
         device = svglite::svglite, bg = paper)
  tmp_tiff <- tempfile(fileext = ".tiff")
  ggsave(tmp_tiff, plot, width = width_mm, height = height_mm, units = "mm",
         device = ragg::agg_tiff, dpi = 600, compression = "lzw", background = paper)
  if (!file.copy(tmp_tiff, file.path(out, paste0(stem, ".tiff")), overwrite = TRUE)) {
    stop("TIFF copy failed: ", stem)
  }
}

# Integrated secretory atlas: opaque high-contrast points, rare states drawn last.
emb <- read.csv(file.path(root, "results/advanced_singlecell/integrated_secretory_embedding.csv"),
                check.names = FALSE)
counts <- table(emb$subtype)
draw_levels <- names(sort(counts, decreasing = TRUE))
emb$draw_order <- match(emb$subtype, draw_levels)
emb <- emb[order(emb$draw_order), ]
emb$subtype <- factor(emb$subtype, levels = state_order)
legend_labels <- paste0(state_labels[state_order], "  (n=", comma(as.numeric(counts[state_order])), ")")
xr <- range(emb$tSNE1, na.rm = TRUE)
yr <- range(emb$tSNE2, na.rm = TRUE)
p <- ggplot(emb, aes(tSNE1, tSNE2, colour = subtype)) +
  geom_point(size = 0.48, alpha = 0.88, stroke = 0) +
  annotate("label", x = xr[1] + 0.015 * diff(xr), y = yr[2] - 0.015 * diff(yr),
           label = paste0(comma(nrow(emb)), " cells  |  2 cohorts"),
           hjust = 0, vjust = 1, family = font, size = 2.65,
           linewidth = 0, fill = alpha("white", 0.88), colour = ink) +
  scale_colour_manual(values = state_cols, breaks = state_order,
                      labels = legend_labels, name = NULL, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.035)) +
  scale_y_continuous(expand = expansion(mult = 0.035)) +
  coord_equal(clip = "off") +
  theme_void(base_family = font, base_size = 9) +
  theme(
    plot.background = element_rect(fill = paper, colour = NA),
    legend.position = "right",
    legend.justification = "center",
    legend.key.height = unit(4.8, "mm"),
    legend.text = element_text(size = 8.5, colour = ink),
    plot.margin = margin(6, 7, 6, 7, "mm")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3.0, alpha = 1)))
pA <- p
save_panel(p, "secretory_atlas_v2", 190, 128)

# Bidirectional label transfer: emphasize the diagonal and principal-state recall.
cf <- read.csv(file.path(root, "results/advanced_singlecell/cross_dataset_label_transfer_confusion.csv"))
short <- c(
  "Serous_glandular" = "Serous", "Mucous_glandular" = "Mucous",
  "Surface_club" = "Surface club", "Goblet" = "Goblet",
  "Inflammatory_secretory" = "Inflammatory",
  "Transitional_secretory" = "Transitional", "Ambiguous_secretory" = "Ambiguous"
)
short_order <- unname(short[state_order])
cf$reference <- factor(unname(short[cf$truth]), levels = short_order)
cf$transferred <- factor(unname(short[cf$predicted]), levels = rev(short_order))
cf$direction <- ifelse(
  cf$train_dataset == "GSE235711",
  "GSE235711 reference  ->  GSE276503 query\nSerous recall: 81.7%",
  "GSE276503 reference  ->  GSE235711 query\nSerous recall: 89.2%"
)
diag_cf <- cf[cf$truth == cf$predicted, ]
p <- ggplot(cf, aes(reference, transferred, fill = fraction)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_tile(data = diag_cf, fill = NA, colour = ink, linewidth = 0.78) +
  geom_text(aes(label = ifelse(fraction >= 0.015, percent(fraction, accuracy = 1), ""),
                colour = fraction >= 0.48), family = font, size = 2.45, fontface = "bold") +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       limits = c(0, 1), labels = label_percent(), name = "Within-reference\nfraction") +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = ink), guide = "none") +
  coord_fixed() +
  labs(x = "Reference identity", y = "Transferred identity") +
  theme_hi(8.2) +
  theme(
    axis.text.x = element_text(angle = 42, hjust = 1, size = 7.4),
    axis.text.y = element_text(size = 7.4),
    axis.ticks = element_blank(), axis.line = element_blank(),
    panel.spacing = unit(7, "mm"), strip.text = element_text(size = 8.2, lineheight = 1.05),
    legend.key.height = unit(18, "mm")
  )
pB <- p
save_panel(p, "bidirectional_label_transfer_v2", 230, 118)

# Donor-level composition heatmap with anatomy-aware row grouping.
comp <- read.csv(file.path(root, "results/secretory_subtypes/donor_subtype_composition.csv"))
hm <- comp[comp$dataset == "GSE276503" & comp$subtype %in% state_order, ]
hm$group <- ifelse(hm$disease == "Healthy", "Healthy | inferior turbinate",
                   ifelse(hm$tissue == "Inferior_turbinate", "CRSwNP | inferior turbinate",
                          ifelse(hm$tissue == "Middle_turbinate", "CRSwNP | middle turbinate",
                                 "CRSwNP | nasal polyp")))
group_order <- c("Healthy | inferior turbinate", "CRSwNP | inferior turbinate",
                 "CRSwNP | middle turbinate", "CRSwNP | nasal polyp")
row_key <- unique(hm[, c("donor", "group")])
row_key$group_rank <- match(row_key$group, group_order)
row_key <- row_key[order(row_key$group_rank, row_key$donor), ]
hm$donor <- factor(hm$donor, levels = rev(row_key$donor))
hm$subtype <- factor(hm$subtype, levels = state_order)
hm$group <- factor(hm$group, levels = group_order)
p <- ggplot(hm, aes(subtype, donor, fill = proportion)) +
  geom_tile(colour = "white", linewidth = 0.62) +
  geom_text(aes(label = ifelse(proportion >= 0.10, percent(proportion, accuracy = 1), ""),
                colour = proportion >= 0.50), family = font, size = 2.05, fontface = "bold") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       limits = c(0, 1), labels = label_percent(), name = "Cell fraction") +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = ink), guide = "none") +
  scale_x_discrete(labels = state_labels[state_order]) +
  labs(x = NULL, y = NULL) +
  theme_hi(7.2) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 7.0),
    axis.text.y = element_text(size = 6.2),
    axis.ticks = element_blank(), axis.line = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7.0),
    panel.spacing.y = unit(1.2, "mm")
  )
pC <- p
save_panel(p, "GSE276503_secretory_composition_heatmap_v2", 185, 178)

# External-cohort Serous fractions: distribution, raw donors, medians and sample sizes.
ser <- hm[hm$subtype == "Serous_glandular", ]
ser$display_group <- factor(as.character(ser$group), levels = group_order,
                            labels = c("Healthy\nIT", "CRSwNP\nIT", "CRSwNP\nMT", "CRSwNP\npolyp"))
group_cols <- c("Healthy\nIT" = "#AFCBE3", "CRSwNP\nIT" = "#73C1AE",
                "CRSwNP\nMT" = "#E4C66C", "CRSwNP\npolyp" = "#E9855A")
group_n <- table(ser$display_group)
med <- aggregate(proportion ~ display_group, ser, median)
p <- ggplot(ser, aes(display_group, proportion, fill = display_group)) +
  geom_violin(width = 0.78, scale = "width", trim = TRUE, alpha = 0.25,
              colour = NA) +
  geom_boxplot(width = 0.42, outlier.shape = NA, linewidth = 0.58,
               alpha = 0.78, colour = ink) +
  geom_point(position = position_jitter(width = 0.10, seed = 20260831),
             shape = 21, size = 2.6, stroke = 0.48, colour = ink, alpha = 0.94) +
  geom_point(data = med, shape = 23, size = 3.2, fill = "white", colour = ink,
             stroke = 0.65) +
  annotate("segment", x = 2, xend = 4, y = 1.035, yend = 1.035, colour = ink, linewidth = 0.50) +
  annotate("segment", x = 2, xend = 2, y = 1.005, yend = 1.035, colour = ink, linewidth = 0.50) +
  annotate("segment", x = 4, xend = 4, y = 1.005, yend = 1.035, colour = ink, linewidth = 0.50) +
  annotate("text", x = 3, y = 1.075, label = "P = 0.0062  |  FDR = 0.0109",
           family = font, size = 3.0, colour = ink) +
  annotate("text", x = 1:4, y = -0.055,
           label = paste0("n = ", as.numeric(group_n)), family = font, size = 2.7, colour = muted) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_y_continuous(limits = c(-0.08, 1.12), breaks = seq(0, 1, 0.25),
                     labels = label_percent(), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Serous fraction among secretory cells") +
  theme_hi(9) +
  theme(axis.text.x = element_text(size = 8.2, face = "bold"))
pD <- p
save_panel(p, "GSE276503_serous_fraction_v2", 148, 112)

# State-specific composition effects with a visual no-change band and effect labels.
eff <- read.csv(file.path(root, "results/frozen_v1/formal_composition_effects.csv"))
seff <- eff[eff$dataset == "GSE276503" & eff$contrast == "NP_vs_CRSwNP_IT", ]
seff$state_label <- unname(state_labels[seff$subtype])
seff$state_label <- factor(seff$state_label, levels = rev(unname(state_labels[state_order])))
seff$direction <- ifelse(seff$median_difference < 0, "Lower in polyp", "Higher in polyp")
seff$effect_text <- sprintf("%+.1f pp", 100 * seff$median_difference)
p <- ggplot(seff, aes(median_difference, state_label)) +
  annotate("rect", xmin = -0.05, xmax = 0.05, ymin = -Inf, ymax = Inf,
           fill = "#F0F3F5", colour = NA) +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high, colour = direction),
                orientation = "y", width = 0.17, linewidth = 0.78) +
  geom_point(aes(fill = subtype), shape = 21, size = 3.6,
             colour = ink, stroke = 0.55) +
  geom_text(aes(x = 0.645, label = effect_text, colour = direction),
            hjust = 1, family = font, size = 2.75, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = state_cols, guide = "none") +
  scale_colour_manual(values = c("Lower in polyp" = blue, "Higher in polyp" = vermillion),
                      name = NULL) +
  scale_x_continuous(limits = c(-0.69, 0.66), breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels = label_percent(accuracy = 1), expand = expansion(mult = 0)) +
  labs(x = "Median fraction difference: polyp minus CRSwNP IT", y = NULL) +
  theme_hi(8.5) +
  theme(
    legend.position = "top", legend.justification = "left",
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.35),
    axis.text.y = element_text(size = 8.2)
  )
pE <- p
save_panel(p, "GSE276503_subtype_effects_v2", 172, 112)

# Paired donor slope graph with donor-specific direction and percentage-point change.
pair <- read.csv(file.path(root, "deliverables/stage_03_full_figure_set/tables/paired_serous_composition_data.csv"))
pair$tissue <- factor(pair$tissue, levels = c("Ethmoid", "Polyp"),
                      labels = c("Ethmoid", "Nasal polyp"))
wide <- reshape(pair[, c("donor", "tissue", "proportion_of_candidates")],
                idvar = "donor", timevar = "tissue", direction = "wide")
names(wide) <- sub("proportion_of_candidates\\.", "", names(wide))
wide$delta <- wide$`Nasal polyp` - wide$Ethmoid
wide$direction <- ifelse(wide$delta < 0, "Decreased", "Increased")
pair <- merge(pair, wide[, c("donor", "delta", "direction")], by = "donor")
pair$label <- ifelse(pair$tissue == "Nasal polyp",
                     paste0(pair$donor, "  ", sprintf("%+.1f pp", 100 * pair$delta)), "")
pair$label_y <- pair$proportion_of_candidates
pair$label_y[pair$tissue == "Nasal polyp" & pair$donor == "CRSwNP_1"] <- 0.218
pair$label_y[pair$tissue == "Nasal polyp" & pair$donor == "CRSwNP_6"] <- 0.162
p <- ggplot(pair, aes(tissue, proportion_of_candidates, group = donor, colour = direction)) +
  geom_line(linewidth = 1.05, alpha = 0.88) +
  geom_point(aes(fill = tissue), shape = 21, size = 3.7,
             stroke = 0.58, colour = ink) +
  geom_text(aes(y = label_y, label = label), hjust = -0.08, family = font, size = 2.75,
            colour = ink, fontface = "bold", show.legend = FALSE) +
  annotate("label", x = 1.02, y = 0.97,
           label = "Median paired change: -17.4 pp\nWilcoxon P = 0.625",
           hjust = 0, vjust = 1, family = font, size = 2.85,
           linewidth = 0, fill = alpha("white", 0.90), colour = ink) +
  scale_colour_manual(values = c("Decreased" = blue, "Increased" = vermillion), name = NULL) +
  scale_fill_manual(values = c("Ethmoid" = "#AFCBE3", "Nasal polyp" = "#E9855A"), guide = "none") +
  scale_x_discrete(expand = expansion(add = c(0.18, 0.72))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = label_percent()) +
  labs(x = NULL, y = "Serous fraction among secretory cells") +
  theme_hi(9) +
  theme(
    legend.position = "top", legend.justification = "left",
    axis.text.x = element_text(size = 8.5, face = "bold"),
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.35)
  ) +
  coord_cartesian(clip = "off")
pF <- p
save_panel(p, "GSE235711_paired_serous_fraction_v2", 155, 115)

# Serous effect forest: highlight the prespecified primary contrast and show estimates.
forest <- eff[eff$subtype == "Serous_glandular", ]
contrast_labels <- c(
  "CRSsNP_vs_Healthy_Ethmoid" = "CRSsNP vs healthy | ethmoid",
  "CRSwNP_vs_Healthy_Ethmoid" = "CRSwNP vs healthy | ethmoid",
  "CRSwNP_vs_Healthy_IT" = "CRSwNP vs healthy | inferior turbinate",
  "AR_vs_Healthy_IT" = "AR vs healthy | inferior turbinate",
  "NP_vs_CRSwNP_IT" = "Polyp vs CRSwNP IT | external",
  "paired_NP_vs_Ethmoid" = "Polyp vs ethmoid | paired"
)
ord <- c("CRSsNP_vs_Healthy_Ethmoid", "CRSwNP_vs_Healthy_Ethmoid",
         "CRSwNP_vs_Healthy_IT", "AR_vs_Healthy_IT",
         "paired_NP_vs_Ethmoid", "NP_vs_CRSwNP_IT")
forest <- forest[match(ord, forest$contrast), ]
forest$y <- seq_len(nrow(forest))
forest$display <- paste0(forest$dataset, "  |  ", unname(contrast_labels[forest$contrast]))
forest$primary <- forest$contrast == "NP_vs_CRSwNP_IT"
forest$effect_text <- sprintf("%+.1f pp", 100 * forest$median_difference)
p <- ggplot(forest, aes(median_difference, y)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = forest$y[forest$primary] - 0.42, ymax = forest$y[forest$primary] + 0.42,
           fill = "#FFF2E8", colour = NA) +
  geom_vline(xintercept = 0, colour = muted, linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high, colour = primary),
                orientation = "y", width = 0.17, linewidth = 0.78) +
  geom_point(aes(fill = primary), shape = 21, size = 3.5,
             colour = ink, stroke = 0.55) +
  geom_text(aes(x = 0.405, label = effect_text),
            hjust = 1, family = font, size = 2.75, colour = ink, fontface = "bold") +
  scale_colour_manual(values = c("FALSE" = "#59758A", "TRUE" = vermillion), guide = "none") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = vermillion), guide = "none") +
  scale_y_continuous(breaks = forest$y, labels = forest$display, expand = expansion(add = 0.55)) +
  scale_x_continuous(limits = c(-0.70, 0.42), breaks = c(-0.5, -0.25, 0, 0.25),
                     labels = label_percent(accuracy = 1), expand = expansion(mult = 0)) +
  labs(x = "Median difference in Serous fraction", y = NULL) +
  theme_hi(8.3) +
  theme(
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.35),
    axis.text.y = element_text(size = 7.7),
    axis.ticks.y = element_blank()
  )
pG <- p
save_panel(p, "serous_composition_effect_forest_v2", 192, 112)

# Complete Figure 2 layout. Panel letters are added only to this assembled version.
panel_grob <- function(plot, label) {
  grobTree(
    rectGrob(gp = gpar(fill = paper, col = NA)),
    ggplotGrob(plot),
    textGrob(label, x = unit(1.6, "mm"), y = unit(1, "npc") - unit(1.4, "mm"),
             just = c("left", "top"),
             gp = gpar(fontfamily = font, fontsize = 14, fontface = "bold", col = ink)),
    vp = viewport(clip = "on")
  )
}

pA_complete <- pA +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.key.height = unit(3.4, "mm"), legend.text = element_text(size = 7.1)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 2.6, alpha = 1)))
pB_complete <- pB + coord_cartesian(clip = "off") +
  theme(legend.position = "none")
pF_complete <- pF
pF_complete$layers[[3]] <- NULL

layout_matrix <- rbind(
  rep(1, 10),
  rep(2, 10),
  c(rep(3, 6), rep(4, 4)),
  c(rep(5, 6), rep(6, 4)),
  rep(7, 10)
)
complete <- gridExtra::arrangeGrob(
  grobs = list(panel_grob(pA_complete, "A"), panel_grob(pB_complete, "B"),
               panel_grob(pC, "C"), panel_grob(pD, "D"),
               panel_grob(pE, "E"), panel_grob(pF_complete, "F"),
               panel_grob(pG, "G")),
  layout_matrix = layout_matrix,
  widths = unit(rep(22, 10), "mm"),
  heights = unit(c(90, 60, 80, 55, 45), "mm"),
  padding = unit(1.2, "mm")
)

if (!skip_exports) {
  ggsave(file.path(out, "Figure2_complete_layout.pdf"), complete,
         width = 220, height = 330, units = "mm", device = cairo_pdf,
         family = font, bg = paper)
  ggsave(file.path(out, "Figure2_complete_layout.svg"), complete,
         width = 220, height = 330, units = "mm", device = svglite::svglite,
         bg = paper)
  tmp_complete <- tempfile(fileext = ".tiff")
  ggsave(tmp_complete, complete, width = 220, height = 330, units = "mm",
         device = ragg::agg_tiff, dpi = 600, compression = "lzw", background = paper)
  if (!file.copy(tmp_complete, file.path(out, "Figure2_complete_layout.tiff"), overwrite = TRUE)) {
    stop("Complete Figure 2 TIFF copy failed")
  }
  tmp_preview <- tempfile(fileext = ".png")
  ggsave(tmp_preview, complete, width = 220, height = 330, units = "mm",
         device = ragg::agg_png, dpi = 220, background = paper)
  if (!file.copy(tmp_preview, file.path(out, "Figure2_complete_layout_preview.png"), overwrite = TRUE)) {
    stop("Complete Figure 2 preview copy failed")
  }
}

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
