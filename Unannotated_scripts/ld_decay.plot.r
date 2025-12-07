#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

# --- Defaults ---
ld_path   <- ifelse(length(args) >= 1, args[1], "ld_decay.ld.gz")
bin_kb    <- ifelse(length(args) >= 2, as.numeric(args[2]), 5)    # bin width in kb
max_kb    <- ifelse(length(args) >= 3, as.numeric(args[3]), 500)  # analyze distances <= max_kb
subsample <- ifelse(length(args) >= 4, as.numeric(args[4]), 2e6)  # speed guard for huge files
r2_cutoff <- 0.1

message("Input: ", ld_path)
message("bin_kb=", bin_kb, " | max_kb=", max_kb, " | subsample=", subsample, " | r2_cutoff=", r2_cutoff)
message("Loading LD table...")

# --- Load LD table (PLINK --r2 gz) ---
# Expected columns (PLINK 1.9/2.0): CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2 (plus optional extras)
ld <- fread(ld_path, sep = "\t", header = TRUE, showProgress = TRUE)

message("Sanity checks and standardization...")
# Sanity checks and standardization
need <- c("BP_A", "BP_B", "R2")
if (!all(need %in% names(ld))) {
  stop("Input file must have columns: ", paste(need, collapse = ", "),
       ". Got: ", paste(names(ld), collapse = ", "))
}

message("Cleaning and deriving distance...")
# Clean and derive distance
ld <- ld[is.finite(R2)]
ld[, dist_bp := abs(BP_B - BP_A)]
ld <- ld[dist_bp > 0 & dist_bp <= max_kb * 1000]

# Optional subsample for speed/plotting
if (nrow(ld) > subsample) {
  set.seed(1L)
  ld <- ld[sample(.N, subsample)]
  message("Subsampled to ", subsample, " rows for speed.")
}

message("Binning...")
# --- Bin distances and compute median r^2 per bin ---
bin_size <- bin_kb * 1000
ld[, bin_start := (dist_bp %/% bin_size) * bin_size ]
agg <- ld[, .(
  median_r2 = median(R2, na.rm = TRUE),
  n_pairs   = .N
), by = bin_start][order(bin_start)]
agg[, bin_mid_bp := bin_start + bin_size / 2]
agg[, bin_mid_kb := bin_mid_bp / 1000]

message("Establishing threshold...")
# --- Find first distance where median r^2 <= cutoff ---
thr_row <- agg[median_r2 <= r2_cutoff][1]
threshold_kb <- if (nrow(thr_row) == 1) thr_row$bin_mid_kb else NA_real_

if (is.na(threshold_kb)) {
  message("Median r^2 did not drop below ", r2_cutoff, " within ", max_kb, " kb.")
} else {
  message(sprintf("Median r^2 drops below %.2f at ~%.1f kb.", r2_cutoff, threshold_kb))
}

message("Saving file...")
# --- Save summary table ---
fwrite(agg[, .(bin_mid_kb, median_r2, n_pairs)], "ld_decay_median_by_bin.tsv", sep = "\t")

message("Plotting...")
# --- Plot ---
# Light scatter of r^2 vs distance + median line + threshold markers
# (Downsample points again for plotting if still huge)
plot_pts <- ld
if (nrow(plot_pts) > 2e5) {
  set.seed(2L)
  plot_pts <- plot_pts[sample(.N, 2e5)]
}

p <- ggplot() +
  geom_point(
    data = plot_pts,
    aes(x = dist_bp/1000, y = R2),
    alpha = 0.08, size = 0.3
  ) +
  geom_line(
    data = agg,
    aes(x = bin_mid_kb, y = median_r2),
    linewidth = 1
  ) +
  geom_hline(yintercept = r2_cutoff, linetype = "dashed") +
  {
    if (!is.na(threshold_kb))
      geom_vline(xintercept = threshold_kb, linetype = "dotted")
  } +
  {
    if (!is.na(threshold_kb))
      annotate("text", x = threshold_kb, y = r2_cutoff + 0.03,
               label = sprintf("~%.1f kb @ r^2=%.2f", threshold_kb, r2_cutoff),
               hjust = -0.05, vjust = 0, size = 3.5)
  } +
  labs(
    x = "Pairwise distance (kb)",
    y = expression(r^2),
    title = "LD decay: r² vs genomic distance",
    subtitle = paste0("Binned median (", bin_kb, " kb) • up to ", max_kb, " kb")
  ) +
  theme_bw(base_size = 12)

ggsave("ld_decay_plot.pdf", p, width = 7, height = 5, device = cairo_pdf)
ggsave("ld_decay_plot.png", p, width = 7, height = 5, dpi = 200)


# --- Also write a tiny summary text file ---
sum_lines <- c(
  paste0("input\t", ld_path),
  paste0("bin_kb\t", bin_kb),
  paste0("max_kb\t", max_kb),
  paste0("r2_cutoff\t", r2_cutoff),
  paste0("threshold_kb\t", ifelse(is.na(threshold_kb), "NA", sprintf("%.1f", threshold_kb)))
)
writeLines(sum_lines, "ld_decay_summary.txt")

message("Done!")