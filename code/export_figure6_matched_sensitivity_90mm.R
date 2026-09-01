options(stringsAsFactors = FALSE)
options(error = function() { traceback(2); q(status = 1) })

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ragg)
  library(grid)
})

root <- "."
out <- file.path(root, "deliverables", "Figure6_matched_permutation_sensitivity_90mm")
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
light <- "#DCE6EB"

theme_final <- function(base_size = 8.4) {
  theme_classic(base_family = font_family, base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = ink),
      axis.ticks = element_line(linewidth = 0.30, colour = ink),
      axis.text = element_text(size = base_size - 0.2, colour = ink),
      axis.title = element_text(size = base_size, face = "bold", colour = ink),
      plot.title = element_text(size = base_size + 0.4, face = "bold", colour = ink,
                                margin = margin(0, 0, 2.2, 0, "mm")),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#5F7180",
                                   margin = margin(0, 0, 2.5, 0, "mm")),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(3.0, 3.5, 3.0, 3.5, "mm")
    )
}

save_final <- function(p, stem, height_mm) {
  tf <- file.path(out, paste0(stem, ".tiff"))
  ragg::agg_tiff(tf, width = 90, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white")
  print(p)
  dev.off()

  pf <- file.path(out, paste0(stem, ".pdf"))
  cairo_pdf(pf, width = 90 / 25.4, height = height_mm / 25.4,
            family = font_family, bg = "white")
  print(p)
  dev.off()

  qf <- file.path(qa_dir, paste0(stem, "_QA96.png"))
  ragg::agg_png(qf, width = 90, height = height_mm, units = "mm", res = 96,
                background = "white")
  print(p)
  dev.off()
}

summary_tab <- read.csv(file.path(root, "results", "advanced_coexpression",
                                  "module_preservation_summary.csv"))
perm <- read.csv(file.path(root, "results", "advanced_coexpression",
                           "module_density_matched_permutation.csv"))
strata <- read.csv(file.path(root, "results", "advanced_coexpression",
                             "module_density_matched_strata_audit.csv"))
metrics <- read.csv(file.path(root, "results", "advanced_coexpression",
                              "module_density_matching_metrics.csv"))

observed <- summary_tab$validation_absolute_density[[1]]
empirical_p <- summary_tab$validation_density_matched_empirical_p[[1]]
seed <- summary_tab$permutation_seed[[1]]
n_perm <- nrow(perm)
null_q <- quantile(perm$matched_permuted_validation_density,
                   probs = c(0.025, 0.5, 0.975), names = FALSE)

write.csv(summary_tab, file.path(data_dir, "FigS18_summary_statistics.csv"), row.names = FALSE)
write.csv(perm, file.path(data_dir, "FigS18B_matched_density_permutation.csv"), row.names = FALSE)
write.csv(strata, file.path(data_dir, "FigS18A_matching_strata.csv"), row.names = FALSE)
write.csv(metrics, file.path(data_dir, "FigS18_matching_metrics.csv"), row.names = FALSE)

# A: the exact 5 x 5 matching structure. Bubble area is the number of locked
# module genes that every random draw must reproduce in that stratum.
strata$expression_bin <- as.integer(sub("_.*$", "", strata$stratum)) + 1L
strata$variance_bin <- as.integer(sub("^.*_", "", strata$stratum)) + 1L
grid_all <- expand.grid(expression_bin = 1:5, variance_bin = 1:5)
strata_plot <- merge(grid_all,
                     strata[, c("expression_bin", "variance_bin", "module_genes", "candidate_genes")],
                     by = c("expression_bin", "variance_bin"), all.x = TRUE)
strata_plot$module_genes[is.na(strata_plot$module_genes)] <- 0
strata_plot$candidate_genes[is.na(strata_plot$candidate_genes)] <- 0
strata_plot$label <- ifelse(strata_plot$module_genes > 0,
                            strata_plot$module_genes, "")

pA <- ggplot(strata_plot, aes(expression_bin, variance_bin)) +
  geom_point(aes(size = module_genes, fill = module_genes), shape = 21,
             colour = "white", stroke = 0.55, alpha = 0.97) +
  geom_text(aes(label = label), family = font_family, fontface = "bold",
            colour = "white", size = 2.45) +
  scale_size_area(max_size = 13, limits = c(0, max(strata_plot$module_genes)),
                  breaks = c(5, 20, 40, 60), name = "Module genes") +
  scale_fill_gradientn(colours = c(cyan, teal, blue, navy),
                       limits = c(0, max(strata_plot$module_genes)), guide = "none") +
  scale_x_continuous(breaks = 1:5, labels = c("Low", "2", "3", "4", "High"),
                     expand = expansion(mult = c(0.10, 0.10))) +
  scale_y_continuous(breaks = 1:5, labels = c("Low", "2", "3", "4", "High"),
                     expand = expansion(mult = c(0.12, 0.12))) +
  labs(title = "Expression- and variance-matched sampling",
       subtitle = "Every draw preserves the frozen module's\n5 x 5 stratum counts",
       x = "Mean expression rank", y = "Residual variance rank") +
  coord_fixed() +
  theme_final(8.4) +
  theme(
    panel.grid.major = element_line(colour = light, linewidth = 0.35),
    axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.height = unit(6.5, "mm"),
    legend.box.margin = margin(0, 0, 0, -1, "mm")
  ) +
  guides(size = guide_legend(override.aes = list(fill = teal, colour = "white", alpha = 1)))
save_final(pA, "FigS18A_matching_strata_90mm", 78)

# B: matched null distribution. A density silhouette and a deterministic
# subsample of random draws preserve distributional information at final size.
set.seed(seed)
perm$sampled <- FALSE
perm$sampled[sample(seq_len(n_perm), min(420L, n_perm))] <- TRUE
density_obj <- density(perm$matched_permuted_validation_density, n = 512)
density_df <- data.frame(x = density_obj$x,
                         y = density_obj$y / max(density_obj$y))
jitter_df <- perm[perm$sampled, , drop = FALSE]
jitter_df$y <- runif(nrow(jitter_df), -0.085, -0.020)

pB <- ggplot() +
  geom_ribbon(data = density_df, aes(x = x, ymin = 0, ymax = y),
              fill = alpha(teal, 0.72), colour = NA) +
  geom_line(data = density_df, aes(x = x, y = y), colour = navy, linewidth = 0.65) +
  geom_point(data = jitter_df,
             aes(x = matched_permuted_validation_density, y = y),
             shape = 16, size = 0.62, alpha = 0.28, colour = blue) +
  annotate("segment", x = null_q[1], xend = null_q[3], y = 0.08, yend = 0.08,
           colour = "white", linewidth = 1.20, lineend = "round") +
  annotate("point", x = null_q[2], y = 0.08, shape = 21, size = 2.5,
           fill = "white", colour = navy, stroke = 0.55) +
  geom_vline(xintercept = observed, colour = orange, linewidth = 1.05) +
  annotate("label", x = observed, y = 0.92,
           label = paste0("Observed = ", sprintf("%.3f", observed),
                          "\nEmpirical P = ", sprintf("%.6f", empirical_p)),
           hjust = 1.03, vjust = 1, family = font_family, fontface = "bold",
           size = 2.55, lineheight = 1.02, linewidth = 0,
           fill = alpha("white", 0.92), colour = orange) +
  annotate("label", x = 0.242, y = 0.66,
           label = paste0(comma(n_perm), " matched gene sets\n",
                          "Median ", sprintf("%.3f", null_q[2]),
                          "; 95% range ", sprintf("%.3f", null_q[1]),
                          "-", sprintf("%.3f", null_q[3])),
           hjust = 0, family = font_family, size = 2.45, lineheight = 1.05,
           fontface = "bold", linewidth = 0,
           fill = alpha("white", 0.90), colour = ink) +
  scale_x_continuous(labels = label_number(accuracy = 0.01),
                     limits = c(min(density_df$x) - 0.006, observed + 0.008),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.11, 1.02), breaks = NULL, expand = c(0, 0)) +
  labs(title = "The preserved module exceeds the matched null",
       subtitle = "Validation-cohort mean absolute gene-gene correlation",
       x = "Module density | GSE36830", y = NULL) +
  theme_final(8.4) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(3.0, 3.5, 3.5, 3.5, "mm"))
save_final(pB, "FigS18B_matched_density_permutation_90mm", 66)

session_lines <- c(
  paste0("R version: ", R.version.string),
  paste0("ggplot2: ", as.character(packageVersion("ggplot2"))),
  paste0("ragg: ", as.character(packageVersion("ragg"))),
  paste0("scales: ", as.character(packageVersion("scales"))),
  paste0("Seed: ", seed),
  paste0("Matched permutations: ", n_perm),
  "Matching: validation mean expression x within-group residual variance; 5 x 5 rank-quantile strata",
  "Sampling universe: eligible non-module genes with finite positive residual variance",
  "Sampling: without replacement within each stratum; exact module stratum counts preserved"
)
writeLines(session_lines, file.path(out, "R_session_and_analysis_parameters.txt"))

message("Exported Figure S18 matched-permutation sensitivity panels to: ", normalizePath(out))
