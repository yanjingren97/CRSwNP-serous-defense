options(stringsAsFactors = FALSE)
options(error = function() { traceback(2); q(status = 1) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ragg)
  library(grid)
  library(igraph)
  library(ggrepel)
})

root <- "."
out <- Sys.getenv("FIG6_FINAL_RASTER_DIR",
                  unset = file.path(root, "deliverables", "Figure6_final_90mm"))
qa_dir <- file.path(out, "QA_final_size_96dpi")
data_dir <- file.path(out, "panel_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#243442"
navy <- "#173B57"
blue <- "#357FA1"
teal <- "#2A9D8F"
cyan <- "#69B7C4"
orange <- "#E4772A"
red <- "#BE4B42"
gold <- "#E6B655"
grey <- "#758692"
light <- "#DCE6EB"

theme_final <- function(base_size = 8.2) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = ink),
      axis.ticks = element_line(linewidth = 0.30, colour = ink),
      axis.text = element_text(size = base_size - 0.2, colour = ink),
      axis.title = element_text(size = base_size, face = "bold", colour = ink),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.3, face = "bold"),
      legend.text = element_text(size = base_size - 0.6),
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

save_ht_final <- function(draw_fun, stem, height_mm) {
  tf <- file.path(out, paste0(stem, ".tiff"))
  ragg::agg_tiff(tf, width = 90, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white")
  draw_fun(); dev.off()
  qf <- file.path(qa_dir, paste0(stem, "_QA96.png"))
  ragg::agg_png(qf, width = 90, height = height_mm, units = "mm", res = 96,
                background = "white")
  draw_fun(); dev.off()
}

summary_tab <- read.csv(file.path(root, "results", "advanced_coexpression",
                                  "module_preservation_summary.csv"))
nodes <- read.csv(file.path(root, "results", "advanced_coexpression",
                            "preserved_module_nodes.csv"))
edges <- read.csv(file.path(root, "results", "advanced_coexpression",
                            "preserved_module_edges.csv"))
all_pairs <- read.csv(file.path(root, "results", "advanced_coexpression",
                                "all_module_gene_pairs.csv"))
all_pairs$is_consensus <- tolower(as.character(all_pairs$is_consensus)) == "true"
perm <- read.csv(file.path(root, "results", "advanced_coexpression",
                           "module_density_permutation.csv"))
nodes$locked <- tolower(as.character(nodes$is_locked)) == "true"
nodes$class <- ifelse(nodes$locked, "Locked defense gene", "Expanded neighbor")
write.csv(summary_tab, file.path(data_dir, "Fig6_summary_statistics.csv"), row.names = FALSE)
write.csv(nodes, file.path(data_dir, "Fig6A_gene_module_data.csv"), row.names = FALSE)
write.csv(all_pairs, file.path(data_dir, "Fig6B_all_19503_gene_pairs.csv"), row.names = FALSE)
write.csv(edges, file.path(data_dir, "Fig6B_consensus_7389_edges.csv"), row.names = FALSE)
write.csv(perm, file.path(data_dir, "Fig6C_density_permutation_data.csv"), row.names = FALSE)

# A: gene-module correlation preservation. Degree is encoded by point area but
# intentionally has no second legend: at 90 mm, preserving the data field is
# more valuable than a redundant size key.
pA <- ggplot(nodes, aes(discovery_module_cor, validation_module_cor)) +
  geom_hline(yintercept = 0, colour = light, linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = light, linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9FB0BA", linewidth = 0.40) +
  stat_density_2d(colour = "#B9CBD4", linewidth = 0.32, bins = 5) +
  geom_point(aes(size = consensus_degree, fill = class), shape = 21,
             colour = ink, stroke = 0.35, alpha = 0.86) +
  annotate("label", x = 0.235, y = 0.905,
           label = paste0(nrow(nodes), " module genes\n",
                          percent(summary_tab$validation_same_positive_fraction, accuracy = 0.1),
                          " retained positive kME\nSpearman rho = ",
                          sprintf("%.3f", summary_tab$gene_module_correlation_rho)),
           hjust = 0, vjust = 1, family = font_family, fontface = "bold",
           size = 2.55, lineheight = 1.00,
           linewidth = 0, fill = alpha("white", 0.90), colour = ink) +
  scale_fill_manual(values = c("Locked defense gene" = orange, "Expanded neighbor" = cyan), name = NULL) +
  scale_size_continuous(range = c(1.8, 5.4), guide = "none") +
  labs(x = "Gene-module correlation | GSE136825",
       y = "Gene-module correlation | GSE36830") +
  coord_cartesian(xlim = c(0.20, 0.96), ylim = c(-0.22, 0.94), clip = "off") +
  theme_final(8.2) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-2, 0, 0, 0),
        aspect.ratio = 0.90,
        panel.border = element_rect(fill = NA, colour = light, linewidth = 0.35)) +
  guides(fill = guide_legend(direction = "horizontal", nrow = 1,
                             override.aes = list(size = 3.2)))
save_final(pA, "Fig6A_gene_module_reproducibility_90mm", 78)

# B: cross-cohort co-expression reproducibility across the complete edge
# universe. The orange overlay is the thresholded consensus subset used for
# network construction; rho and p are calculated across all 19,503 pairs.
pB <- ggplot(all_pairs, aes(cor_discovery, cor_validation)) +
  geom_bin_2d(bins = 52) +
  geom_point(data = all_pairs[all_pairs$is_consensus, ],
             colour = orange, size = 0.34, alpha = 0.12) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "#8FA2AD", linewidth = 0.42) +
  geom_vline(xintercept = 0.55, linetype = "dotted", colour = orange, linewidth = 0.52) +
  geom_hline(yintercept = 0.30, linetype = "dotted", colour = orange, linewidth = 0.52) +
  scale_fill_gradientn(colours = c("#F6F9FA", "#C9E0E5", cyan, blue, navy),
                       trans = "sqrt", name = "All pairs/bin") +
  annotate("label", x = -0.95, y = 0.95,
           label = paste0("All pairs: ", comma(nrow(all_pairs)), "\n",
                          "Consensus subset: ", comma(sum(all_pairs$is_consensus)),
                          "\nSpearman ρ = ",
                          sprintf("%.3f", summary_tab$edge_correlation_preservation_rho),
                          " | p = 2.31 × 10⁻²⁹⁴"),
           hjust = 0, vjust = 1, family = font_family, fontface = "bold",
           size = 2.62, lineheight = 1.04,
           linewidth = 0, fill = alpha("white", 0.90), colour = ink) +
  annotate("label", x = -0.95, y = -0.93,
           label = "Orange: consensus edges\nr(discovery) ≥ 0.55 and r(validation) ≥ 0.30",
           hjust = 0, vjust = 0, family = font_family, fontface = "bold",
           size = 2.28, lineheight = 1.00, label.padding = unit(0.9, "mm"),
           linewidth = 0, fill = alpha("white", 0.90), colour = "#A9531F") +
  labs(x = "Gene-pair correlation | GSE136825",
       y = "Gene-pair correlation | GSE36830") +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), clip = "off") +
  theme_final(8.6) +
  theme(legend.position = "bottom",
        legend.margin = margin(-1, 0, 0, 0),
        legend.key.width = unit(22, "mm"),
        legend.key.height = unit(2.2, "mm"),
        panel.border = element_rect(fill = NA, colour = light, linewidth = 0.35)) +
  guides(fill = guide_colorbar(direction = "horizontal", title.position = "left",
                                  title.hjust = 0.5, barwidth = unit(22, "mm"),
                                  barheight = unit(2.2, "mm"), ticks = FALSE))
save_final(pB, "Fig6B_edge_reproducibility_landscape_90mm", 84)

# C: co-expression-density reproducibility audit against same-sized random sets.
obs_density <- summary_tab$validation_absolute_density[1]
den <- density(perm$permuted_validation_density, n = 512, bw = "nrd0")
null_curve <- data.frame(x = den$x, ymin = 0, ymax = den$y / max(den$y))
null_ticks <- data.frame(x = perm$permuted_validation_density)
pC <- ggplot(null_curve, aes(x)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#BFDCE2", alpha = 0.95) +
  geom_line(aes(y = ymax), colour = blue, linewidth = 0.65) +
  geom_segment(data = null_ticks, aes(x = x, xend = x, y = -0.025, yend = 0.035),
               inherit.aes = FALSE, colour = ink, linewidth = 0.22, alpha = 0.25) +
  annotate("segment", x = obs_density, xend = obs_density, y = -0.04, yend = 0.78,
           colour = orange, linewidth = 1.05) +
  annotate("label", x = obs_density, y = 0.89,
           label = paste0("Observed density = ", sprintf("%.3f", obs_density),
                          "\nEmpirical p = ", sprintf("%.4f", summary_tab$validation_density_empirical_p)),
           hjust = 1, vjust = 1, family = font_family, fontface = "bold",
           size = 2.75, lineheight = 1.02,
           linewidth = 0, fill = alpha("white", 0.92), colour = ink) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(-0.05, 1.02), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Absolute module density in validation cohort",
       caption = "Null: 500 same-sized random gene sets  |  orange: frozen module") +
  theme_final(8.6) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.caption = element_text(size = 7.2, colour = "#5E707D", hjust = 0))
save_final(pC, "Fig6C_coexpression_density_same_sized_null_90mm", 66)

# D: topological map of the consensus module.
ord_nodes <- nodes[order(-nodes$consensus_degree, -nodes$consensus_strength), ]
selected_genes <- unique(c(head(ord_nodes$gene, 26), head(ord_nodes$gene[ord_nodes$locked], 8)))
net_nodes <- nodes[nodes$gene %in% selected_genes, ]
net_edges <- edges[edges$source %in% selected_genes & edges$target %in% selected_genes, ]
ng <- nrow(net_nodes)
idx <- setNames(seq_len(ng), net_nodes$gene)

# Layout uses the complete induced graph, whereas only the strongest local
# edges are drawn. A deterministic force layout prevents artificial isolates.
g_layout <- graph_from_data_frame(
  net_edges[, c("source", "target", "consensus_weight")], directed = FALSE,
  vertices = data.frame(name = net_nodes$gene, stringsAsFactors = FALSE)
)
set.seed(20260831)
xy <- layout_with_fr(g_layout, weights = E(g_layout)$consensus_weight^3,
                     niter = 2500, grid = "nogrid")
xy <- norm_coords(xy, xmin = -0.86, xmax = 0.86, ymin = -0.82, ymax = 0.82)
net_nodes$x <- xy[match(net_nodes$gene, V(g_layout)$name), 1]
net_nodes$y <- xy[match(net_nodes$gene, V(g_layout)$name), 2]

net_edges <- net_edges[order(-net_edges$consensus_weight), ]
keep_rows <- seq_len(min(62, nrow(net_edges)))
for (g in selected_genes) {
  incident <- which(net_edges$source == g | net_edges$target == g)
  keep_rows <- unique(c(keep_rows, head(incident, 2)))
}
net_edges <- net_edges[sort(keep_rows), ]
net_edges$x <- net_nodes$x[match(net_edges$source, net_nodes$gene)]
net_edges$y <- net_nodes$y[match(net_edges$source, net_nodes$gene)]
net_edges$xend <- net_nodes$x[match(net_edges$target, net_nodes$gene)]
net_edges$yend <- net_nodes$y[match(net_edges$target, net_nodes$gene)]
label_genes <- unique(c(head(ord_nodes$gene, 8), head(ord_nodes$gene[ord_nodes$locked], 4)))
net_lab <- net_nodes[net_nodes$gene %in% label_genes, ]
pD <- ggplot() +
  geom_segment(data = net_edges,
               aes(x, y, xend = xend, yend = yend, alpha = consensus_weight,
                   linewidth = consensus_weight),
               colour = "#8EA8B6", lineend = "round") +
  geom_point(data = net_nodes,
             aes(x, y, size = consensus_degree, fill = class),
             shape = 21, colour = ink, stroke = 0.42) +
  geom_text_repel(data = net_lab, aes(x, y, label = gene),
                  family = font_family, fontface = "bold", size = 2.55,
                  colour = ink, box.padding = 0.34, point.padding = 0.22,
                  min.segment.length = 0, segment.colour = "#AFC0C8",
                  segment.size = 0.22, max.overlaps = Inf, force = 2.4,
                  force_pull = 0.25, max.time = 3, max.iter = 20000,
                  seed = 20260831) +
  scale_fill_manual(values = c("Locked defense gene" = orange, "Expanded neighbor" = cyan), name = NULL) +
  scale_size_continuous(range = c(3.0, 7.5), guide = "none") +
  scale_alpha_continuous(range = c(0.10, 0.58), guide = "none") +
  scale_linewidth_continuous(range = c(0.25, 1.00), guide = "none") +
  coord_equal(xlim = c(-1.14, 1.14), ylim = c(-1.05, 1.05), clip = "off") +
  theme_void(base_family = font_family, base_size = 8.0) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = 7.4, colour = ink),
        plot.margin = margin(5, 6, 3, 6, "mm")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 4)))
save_final(pD, "Fig6D_consensus_coexpression_network_90mm", 96)
write.csv(net_nodes, file.path(data_dir, "Fig6D_network_nodes.csv"), row.names = FALSE)
write.csv(net_edges, file.path(data_dir, "Fig6D_network_edges.csv"), row.names = FALSE)

# E: multi-metric fingerprint of the top hubs. A native ggplot matrix is used
# so font rendering, legends and physical sizing match the other five panels.
hub <- head(ord_nodes, 18)
hub_mat <- cbind(
  `Discovery kME` = hub$discovery_module_cor,
  `Validation kME` = hub$validation_module_cor,
  `Degree pct.` = rank(nodes$consensus_degree, ties.method = "average")[match(hub$gene, nodes$gene)] / nrow(nodes),
  `Strength pct.` = rank(nodes$consensus_strength, ties.method = "average")[match(hub$gene, nodes$gene)] / nrow(nodes)
)
rownames(hub_mat) <- hub$gene
hub_long <- data.frame(
  gene = rep(rownames(hub_mat), times = ncol(hub_mat)),
  metric = rep(colnames(hub_mat), each = nrow(hub_mat)),
  value = as.vector(hub_mat),
  stringsAsFactors = FALSE
)
hub_long$gene <- factor(hub_long$gene, levels = rev(hub$gene))
hub_long$metric_x <- match(hub_long$metric, colnames(hub_mat))
hub_anno <- data.frame(
  gene = factor(hub$gene, levels = rev(hub$gene)),
  class = factor(ifelse(hub$locked, "Locked", "Expanded"),
                 levels = c("Locked", "Expanded")),
  x = 0.28
)
pE <- ggplot(hub_long, aes(metric_x, gene)) +
  geom_tile(aes(fill = value), width = 0.96, height = 0.94,
            colour = "white", linewidth = 0.32) +
  geom_point(data = hub_anno, aes(x, gene, colour = class),
             inherit.aes = FALSE, size = 3.0) +
  scale_fill_gradientn(colours = c("#F5F7F8", cyan, navy),
                       limits = c(0.45, 0.95), oob = squish,
                       breaks = c(0.5, 0.7, 0.9), name = "Scaled metric") +
  scale_colour_manual(values = c("Locked" = orange, "Expanded" = cyan), name = "Class") +
  scale_x_continuous(breaks = c(0.28, 1:4),
                     labels = c("Class", colnames(hub_mat)),
                     limits = c(0.03, 4.52), expand = c(0, 0)) +
  scale_y_discrete(position = "right", expand = expansion(add = 0.1)) +
  theme_minimal(base_family = font_family, base_size = 8.0) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.y = element_text(colour = ink, size = 7.4),
        axis.text.x = element_text(colour = ink, size = 7.0, face = "bold", angle = 42,
                                   hjust = 1, vjust = 1),
        axis.ticks = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.justification = "center", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.6, "mm"),
        legend.title = element_text(size = 7.2, face = "bold"),
        legend.text = element_text(size = 6.9),
        legend.key.height = unit(2.2, "mm"),
        legend.key.width = unit(16, "mm"),
        plot.margin = margin(3, 4, 2, 3, "mm")) +
  guides(fill = guide_colorbar(direction = "horizontal", title.position = "top",
                               barwidth = unit(25, "mm"), barheight = unit(2.2, "mm"),
                               ticks = FALSE, order = 2),
         colour = guide_legend(direction = "horizontal", title.position = "top",
                               nrow = 1, order = 1,
                               override.aes = list(size = 3.1)))
save_final(pE, "Fig6E_hub_multimetric_fingerprint_90mm", 108)
write.csv(hub_long, file.path(data_dir, "Fig6E_hub_metrics.csv"), row.names = FALSE)
write.csv(hub_anno, file.path(data_dir, "Fig6E_hub_class.csv"), row.names = FALSE)

# F: direct frozen-prior bridges from replicated activities to consensus-module genes.
rg <- read.csv(file.path(root, "results", "advanced_regulatory",
                         "regulatory_cross_cohort_summary.csv"))
priority <- rg$source[tolower(as.character(rg$replicated_priority)) == "true"]
ct <- read.csv(gzfile(file.path(root, "results", "advanced_regulatory",
                               "collectri_network_frozen.csv.gz")))
pg <- read.csv(gzfile(file.path(root, "results", "advanced_regulatory",
                               "progeny_network_frozen.csv.gz")))
e1 <- ct[ct$source %in% priority & ct$target %in% nodes$gene, c("source", "target", "weight")]
e2 <- pg[pg$source %in% priority & pg$target %in% nodes$gene, c("source", "target", "weight")]
bridge <- rbind(e1, e2)
bridge$sign <- ifelse(bridge$weight >= 0, "Positive prior", "Negative prior")
regs <- unique(bridge$source)
targs <- unique(bridge$target)
reg_nodes <- data.frame(node = regs, layer = "Regulator/pathway", x = 0,
                        y = seq(length(regs), 1, length.out = length(regs)),
                        stringsAsFactors = FALSE)
tar_nodes <- data.frame(node = targs, layer = "Consensus module target", x = 1,
                        y = seq(length(targs), 1, length.out = length(targs)),
                        stringsAsFactors = FALSE)
bridge_nodes <- rbind(reg_nodes, tar_nodes)
bridge_nodes <- merge(bridge_nodes,
                      nodes[, c("gene", "locked", "consensus_degree")],
                      by.x = "node", by.y = "gene", all.x = TRUE, sort = FALSE)
bridge_nodes$node_class <- ifelse(bridge_nodes$layer == "Regulator/pathway", "Priority activity",
                                  ifelse(bridge_nodes$locked %in% TRUE, "Locked target", "Expanded target"))
bridge_nodes$consensus_degree[is.na(bridge_nodes$consensus_degree)] <- max(nodes$consensus_degree) * 0.72
bridge <- merge(bridge, bridge_nodes[, c("node", "x", "y")], by.x = "source", by.y = "node")
names(bridge)[names(bridge) %in% c("x", "y")] <- c("x1", "y1")
bridge <- merge(bridge, bridge_nodes[, c("node", "x", "y")], by.x = "target", by.y = "node")
names(bridge)[names(bridge) %in% c("x", "y")] <- c("x2", "y2")
pF <- ggplot() +
  geom_curve(data = bridge,
             aes(x1, y1, xend = x2, yend = y2, colour = sign, linewidth = abs(weight)),
             curvature = 0.14, alpha = 0.68) +
  geom_point(data = bridge_nodes,
             aes(x, y, fill = node_class, size = consensus_degree),
             shape = 21, colour = ink, stroke = 0.45) +
  geom_text(data = bridge_nodes[bridge_nodes$x == 0, ],
            aes(x - 0.055, y, label = node), hjust = 1,
            family = font_family, fontface = "bold", size = 2.48, colour = ink) +
  geom_text(data = bridge_nodes[bridge_nodes$x == 1, ],
            aes(x + 0.055, y, label = node), hjust = 0,
            family = font_family, fontface = "bold", size = 2.48, colour = ink) +
  annotate("label", x = 0.50, y = max(c(reg_nodes$y, tar_nodes$y)) + 1.05,
           label = paste0(length(regs), "/", length(priority), " replicated activities directly intersected\n",
                          length(targs), " genes in the cross-cohort module"),
           family = font_family, fontface = "bold", size = 2.48, lineheight = 1.00,
           linewidth = 0, fill = "#F1F5F7", colour = ink) +
  scale_colour_manual(values = c("Positive prior" = teal, "Negative prior" = red), name = NULL) +
  scale_fill_manual(values = c("Priority activity" = gold, "Locked target" = orange,
                               "Expanded target" = cyan), name = NULL) +
  scale_linewidth_continuous(range = c(0.35, 1.35), guide = "none") +
  scale_size_continuous(range = c(4.0, 6.8), guide = "none") +
  scale_x_continuous(limits = c(-0.34, 1.36), breaks = c(0, 1),
                     labels = c("Priority activity", "Consensus module target"), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.35, max(c(reg_nodes$y, tar_nodes$y)) + 1.35), clip = "off") +
  theme_final(8.0) +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), axis.text.x = element_text(face = "bold", size = 7.5),
        legend.text = element_text(size = 7.0),
        legend.position = "bottom", legend.box = "vertical",
        panel.grid = element_blank(), plot.margin = margin(5, 7, 3, 7, "mm")) +
  guides(colour = guide_legend(nrow = 1, order = 1),
         fill = guide_legend(nrow = 1, order = 2, override.aes = list(size = 3.5)))
save_final(pF, "Fig6F_regulator_module_bridge_network_90mm", 108)
write.csv(bridge_nodes, file.path(data_dir, "Fig6F_bridge_nodes.csv"), row.names = FALSE)
write.csv(bridge, file.path(data_dir, "Fig6F_bridge_edges.csv"), row.names = FALSE)

notes <- c(
  "FIGURE 6 FINAL 90-mm EXPORT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "All panels: 90 mm wide, RGB, 600 dpi, LZW TIFF, Arial.",
  "QA PNG files were rendered at 96 dpi at the same physical dimensions.",
  paste0("Module: ", summary_tab$module_genes, " genes; ", summary_tab$consensus_edges,
         " consensus edges; validation density empirical p=",
         signif(summary_tab$validation_density_empirical_p, 5), "."),
  "Figure 6 reports cross-cohort co-expression reproducibility, not a formal network-preservation statistic.",
  "Panel B rho and p use all 19,503 gene pairs; 7,389 thresholded consensus edges build the network.",
  "Panel C null comprises 500 same-sized random gene sets.",
  "Regulator-target bridges are frozen-prior overlaps and do not establish causal regulation.",
  "Bulk coexpression cannot completely separate intrinsic regulation from stable glandular-cell abundance.",
  "",
  capture.output(sessionInfo())
)
writeLines(notes, file.path(out, "R_sessionInfo_and_export_notes.txt"))
