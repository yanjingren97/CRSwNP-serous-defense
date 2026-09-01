options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(scales)
})

project_dir <- "."
font_family <- "Arial"
ink <- "#263442"
grid_col <- "#DCE3E8"

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

rrho <- read.csv(file.path(
  project_dir,
  "deliverables", "Figure3_data_and_R_scripts_package", "01_CSV",
  "04_top_ranked_overlap_RRHO.csv"
))

p_rrho <- ggplot(
  rrho,
  aes(
    paired_top_ranked_percentage,
    external_top_ranked_percentage,
    fill = minus_log10_p
  )
) +
  geom_tile() +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_gradientn(
    colours = c("#F5F7F8", "#BFD6DE", "#4A9A9A", "#F0B85B", "#B94132"),
    name = expression(-log[10](italic(p))),
    limits = c(0, 120),
    oob = squish
  ) +
  scale_x_continuous(
    expand = c(0, 0), breaks = c(10, 25, 40, 50),
    labels = label_percent(scale = 1)
  ) +
  scale_y_continuous(
    expand = c(0, 0), breaks = c(10, 25, 40, 50),
    labels = label_percent(scale = 1)
  ) +
  labs(
    x = "Top-ranked genes: paired cohort",
    y = "Top-ranked genes: external cohort"
  ) +
  theme_sci(9) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    text = element_text(family = font_family),
    axis.text = element_text(family = font_family, size = 5.5),
    axis.title = element_text(family = font_family, size = 6.1, face = "bold"),
    strip.text = element_text(family = font_family, size = 6.4, face = "bold"),
    legend.title = element_text(family = font_family, size = 5.7),
    legend.text = element_text(family = font_family, size = 5.3),
    plot.margin = margin(4, 5, 4, 4),
    panel.spacing = grid::unit(2.0, "mm")
  )

output_dir <- Sys.getenv(
  "FIG3_CORRECTED_DIR",
  unset = file.path(project_dir, "deliverables", "Figure3_corrected_A-F")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "Fig3B_RRHO_landscape_corrected_lowercase_p.tiff")
qa_file <- file.path(
  output_dir, "Fig3B_RRHO_landscape_corrected_lowercase_p_QA.png"
)

agg_tiff(
  output_file,
  width = 90, height = 62, units = "mm", res = 600,
  compression = "lzw", background = "white"
)
print(p_rrho)
dev.off()

agg_png(
  qa_file,
  width = 90, height = 62, units = "mm", res = 300,
  background = "white"
)
print(p_rrho)
dev.off()

message(normalizePath(output_file, winslash = "/", mustWork = TRUE))
message(normalizePath(qa_file, winslash = "/", mustWork = TRUE))
