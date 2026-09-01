options(stringsAsFactors = FALSE)

root <- "."
input_file <- file.path(root, "results", "bulk_crs_validation", "GSE136825_scores.csv")
output_pairs <- file.path(root, "results", "bulk_crs_validation",
                          "GSE136825_paired_module_scores.csv")
output_stats <- file.path(root, "results", "bulk_crs_validation",
                          "GSE136825_paired_module_statistics.csv")

d <- read.csv(input_file, check.names = FALSE)

extract_pair <- function(x, title_prefix, group_name, value_name) {
  q <- x[x$group == group_name & grepl(paste0("^", title_prefix, " "), x$title), ]
  q$pair_id <- sub(paste0("^", title_prefix, " "), "", q$title)
  q <- q[, c("pair_id", "sample", "title", "score")]
  names(q)[2:4] <- paste0(value_name, c("_sample", "_title", "_score"))
  q
}

# The GEO sample titles encode the matched tissue pairs:
# PY <-> SCP and TT <-> SCT, with the same numeric patient identifier.
np <- rbind(
  extract_pair(d, "PY", "CRSwNP_NP", "polyp"),
  extract_pair(d, "TT", "CRSwNP_NP", "polyp")
)
it <- rbind(
  extract_pair(d, "SCP", "CRSwNP_IT", "inferior_turbinate"),
  extract_pair(d, "SCT", "CRSwNP_IT", "inferior_turbinate")
)

np$pair_series <- ifelse(grepl("^PY ", np$polyp_title), "PY-SCP", "TT-SCT")
it$pair_series <- ifelse(grepl("^SCP ", it$inferior_turbinate_title), "PY-SCP", "TT-SCT")

pairs <- merge(np, it, by = c("pair_series", "pair_id"), all = FALSE, sort = FALSE)
pairs$pair_numeric <- suppressWarnings(as.integer(pairs$pair_id))
pairs <- pairs[order(pairs$pair_series, pairs$pair_numeric), ]
pairs$delta_polyp_minus_it <- pairs$polyp_score - pairs$inferior_turbinate_score
pairs$direction <- ifelse(
  pairs$delta_polyp_minus_it < 0, "Lower in polyp",
  ifelse(pairs$delta_polyp_minus_it > 0, "Higher in polyp", "No change")
)

stopifnot(nrow(pairs) == 30L)
stopifnot(sum(pairs$delta_polyp_minus_it < 0) == 26L)

wt <- wilcox.test(
  pairs$polyp_score,
  pairs$inferior_turbinate_score,
  paired = TRUE,
  alternative = "two.sided",
  exact = TRUE,
  correct = FALSE,
  conf.int = FALSE
)

stats <- data.frame(
  dataset = "GSE136825",
  contrast = "CRSwNP polyp minus matched CRSwNP inferior turbinate",
  pairing_rule = "PY-SCP and TT-SCT title-prefix pairs with matching numeric identifier",
  n_pairs = nrow(pairs),
  n_lower_in_polyp = sum(pairs$delta_polyp_minus_it < 0),
  proportion_lower_in_polyp = mean(pairs$delta_polyp_minus_it < 0),
  median_polyp_score = median(pairs$polyp_score),
  median_inferior_turbinate_score = median(pairs$inferior_turbinate_score),
  median_paired_difference_polyp_minus_it = median(pairs$delta_polyp_minus_it),
  wilcoxon_signed_rank_statistic_V = unname(wt$statistic),
  wilcoxon_exact_two_sided_p = wt$p.value,
  wilcoxon_method = unname(wt$method),
  stringsAsFactors = FALSE
)

write.csv(pairs, output_pairs, row.names = FALSE, quote = TRUE)
write.csv(stats, output_stats, row.names = FALSE, quote = TRUE)

print(stats, digits = 15)
message(normalizePath(output_pairs, winslash = "/", mustWork = TRUE))
message(normalizePath(output_stats, winslash = "/", mustWork = TRUE))
