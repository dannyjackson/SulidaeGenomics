#!/usr/bin/env Rscript

## -----------------------------------------------------------
## localFstats_manhattan.R
##
## Usage:
##   Rscript localFstats_manhattan.R \
##       <infile> <convfile> <cutoff> [metric_for_outliers] [outdir]
##
## Examples:
##   # Use 99.99th percentile of f_dM as cutoff (top 0.01%)
##   Rscript localFstats_manhattan.R \
##       BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__50_25.txt \
##       /xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt \
##       0.9999
##
##   # Use absolute cutoff of 0.5 on f_d
##   Rscript localFstats_manhattan.R \
##       somefile.txt chromconversion.txt 0.5 f_d
##
## -----------------------------------------------------------

suppressPackageStartupMessages({
  .libPaths("~/R/library_elgato")
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(rlang)
  library(stringr)
})

## -------------------------------
## Parse arguments
## -------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage:\n",
      "  Rscript localFstats_manhattan.R <infile> <convfile> <cutoff> [metric_for_outliers] [outdir]\n\n",
      "  cutoff in (0,1): quantile (e.g. 0.9999 = 99.99th percentile)\n",
      "  cutoff >= 1:     absolute value cutoff on metric\n",
      "  metric_for_outliers (optional): default = f_dM\n",
      "  outdir (optional): default = current working directory\n",
      sep = "")
  quit(status = 1)
}

infile   <- args[1]  # localFstats file
convfile <- args[2]  # chrom conversion file
cutoff_in <- as.numeric(args[3])

if (is.na(cutoff_in)) {
  stop("cutoff must be numeric (either quantile in (0,1) or absolute value).")
}

metric_for_outliers <- if (length(args) >= 4) args[4] else "f_dM"
outdir              <- if (length(args) >= 5) args[5] else getwd()

## -------------------------------
## Derived names from infile
## -------------------------------
# Expected pattern: <pop_name>_localFstats__<win>.txt
base_name <- basename(infile)
pop_name  <- sub("_localFstats__.*$", "", base_name)
win       <- sub("^.*_localFstats__", "", base_name)
win       <- sub("\\.txt$", "", win)

cat("Input file: ", infile, "\n", sep = "")
cat("Pop name:   ", pop_name, "\n", sep = "")
cat("Window tag: ", win, "\n", sep = "")
cat("Outdir:     ", outdir, "\n", sep = "")
cat("Outlier metric: ", metric_for_outliers, "\n", sep = "")
cat("Cutoff input:   ", cutoff_in, "\n\n", sep = "")

## -------------------------------
## Settings
## -------------------------------
metrics <- c("D", "f_d", "f_dM")    # which metrics to plot
metric_cutoffs <- c(D = 0, `f_d` = 0, `f_dM` = 0)  # horizontal line at 0

color1 <- "#444444"
color2 <- "#1f78b4"

# Desired chromosome order
chromo_levels <- c(1, "1A", 2:4, "4A", 5:33, "Z") |> as.character()

## -------------------------------
## Read data
## -------------------------------
cat("Reading localFstats data...\n")
data <- read.table(infile, header = TRUE, as.is = TRUE)

req_cols <- c("chr", "windowStart", "windowEnd")
if (!all(req_cols %in% names(data))) {
  stop("Input file is missing required columns: ",
       paste(setdiff(req_cols, names(data)), collapse = ", "))
}

data <- data %>%
  mutate(
    position = (windowStart + windowEnd) / 2,
    chr = as.character(chr)
  )

## -------------------------------
## Read chromosome conversion map
## Expect two columns: chromo, chr
## -------------------------------
cat("Reading chromosome conversion map...\n")
conv_guess <- read_csv(
  convfile,
  col_names = c("chromo", "chr"),
  show_col_types = FALSE
) %>%
  transmute(
    chr    = as.character(chr),
    chromo = as.character(chromo)
  )

## -------------------------------
## Join mapping + clean
## -------------------------------
cat("Applying chromosome renaming...\n")
data <- data %>%
  left_join(conv_guess, by = "chr")

if (any(is.na(data$chromo))) {
  missing_ids <- unique(data$chr[is.na(data$chromo)])
  warning("Some scaffolds did not map to a chromosome label. Showing first few:\n",
          paste(head(missing_ids, 10), collapse = ", "))
  # drop unmapped
  data <- data %>% filter(!is.na(chromo))
}

# Remove unwanted scaffolds explicitly (as in original)
data <- data %>%
  filter(!chr %in% c("CM062595.1", "CM062599.1", "CM062600.1", "CM062610.1")) %>%
  mutate(chromo = factor(as.character(chromo), levels = chromo_levels))

## -------------------------------
## Prepare cumulative BP coordinates
## -------------------------------
cat("Preparing data for plotting...\n")
plot_data <- data %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position), .groups = "drop") %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarise(center = mean(BPcum), .groups = "drop")

## -------------------------------
## Plot function
## -------------------------------
plot_metric <- function(metric) {
  stopifnot(metric %in% names(plot_data))
  metric_cutoff <- metric_cutoffs[[metric]]

  metric_dir <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", win))
  dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Generating plot for ", metric, "...\n", sep = "")
  p <- ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
    geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
    scale_color_manual(
      values = rep(c(color1, color2), length.out = length(unique(plot_data$chromo)))
    ) +
    scale_x_continuous(
      labels = axisdf$chromo,
      breaks = axisdf$center,
      guide  = guide_axis(n.dodge = 2)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Chromosome", y = metric) +
    geom_hline(yintercept = metric_cutoff) +
    theme_bw(base_size = 22) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  outfile <- file.path(metric_dir, paste0(pop_name, ".", metric, ".", win, ".sigline.png"))
  ggsave(filename = outfile, plot = p, width = 20, height = 5, units = "in")
  cat("Saved: ", outfile, "\n", sep = "")
}

## -------------------------------
## Plot all requested metrics
## -------------------------------
for (m in metrics) {
  if (m %in% names(plot_data)) {
    plot_metric(m)
  } else {
    warning("Metric '", m, "' not found in input; skipping.")
  }
}

## -------------------------------
## Outlier detection + highlight plot
## -------------------------------
if (!metric_for_outliers %in% names(plot_data)) {
  stop("Outlier metric '", metric_for_outliers, "' is not a column in the data.")
}

metric <- metric_for_outliers

cat("\nComputing outliers for metric: ", metric, "\n", sep = "")

metric_values <- plot_data[[metric]]

if (cutoff_in > 0 && cutoff_in < 1) {
  # Interpret as quantile
  cutoff_value <- as.numeric(quantile(metric_values, probs = cutoff_in, na.rm = TRUE))
  cutoff_type  <- paste0("quantile (p=", cutoff_in, ")")
} else {
  # Interpret as absolute value
  cutoff_value <- cutoff_in
  cutoff_type  <- "absolute"
}

if (is.na(cutoff_value)) {
  stop("Computed cutoff is NA (check metric values and cutoff).")
}

cat("Cutoff type:  ", cutoff_type, "\n", sep = "")
cat("Cutoff value: ", signif(cutoff_value, 6), "\n\n", sep = "")

plot_data <- plot_data %>%
  mutate(is_outlier = .data[[metric]] >= cutoff_value)

out_dir <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", win))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir,
                      paste0(pop_name, ".", metric, ".", win, ".", cutoff_in, ".outliers.tsv"))

plot_data %>%
  filter(is_outlier) %>%
  mutate(cutoff_used = cutoff_value,
         cutoff_type = cutoff_type) %>%
  select(
    chromo, chr, position, BPcum, windowStart, windowEnd,
    D, f_d, f_dM, d_f, cutoff_used, cutoff_type
  ) %>%
  arrange(chromo, position) %>%
  write_tsv(out_file)

cat("Outliers written to: ", out_file, "\n", sep = "")

## -------------------------------
## Plot with outliers highlighted
## -------------------------------
base_palette <- rep(c(color1, color2), length.out = length(levels(plot_data$chromo)))

p_out <- ggplot() +
  geom_point(
    data = subset(plot_data, !is_outlier),
    aes(x = BPcum, y = .data[[metric]], color = chromo),
    alpha = 0.8, size = 1, show.legend = FALSE
  ) +
  scale_color_manual(values = setNames(base_palette, levels(plot_data$chromo))) +
  geom_point(
    data = subset(plot_data, is_outlier),
    aes(x = BPcum, y = .data[[metric]]),
    color = "red3", alpha = 0.95, size = 1.6
  ) +
  scale_x_continuous(
    labels = axisdf$chromo,
    breaks = axisdf$center,
    guide  = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = cutoff_value, linetype = "dashed",
             linewidth = 0.5, color = "red3") +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

outfile_png <- file.path(out_dir,
                         paste0(pop_name, ".", metric, ".", win, ".outliers.sigline.png"))

ggsave(outfile_png, p_out, width = 20, height = 5, units = "in")

cat("Saved outlier plot: ", outfile_png, "\n")
cat("Done.\n")
