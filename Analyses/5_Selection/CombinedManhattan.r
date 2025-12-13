#!/usr/bin/env Rscript

# ============================================================
# Multi-format Manhattan plots (single multi-page PDF)
# Files: 3 DFM + 2 FST (handled automatically)
# Order preserved exactly as listed below
# ============================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(rlang)
  library(stringr)
  library(tidyr)
})

# -------------------------------
# File paths (order preserved)
# -------------------------------
MABO_NABO_DFM <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo/MABO_AtlCar_MABO_IndoPacific_NABO_localFstats__5000_200.txt"
MABO_NABO_FST <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fst.50000.Ztransformed.csv"
BFBO_PEBO_DFM <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/bfbo_pebo/BFBO_GofCA_BFBO_southern_PEBO_localFstats__5000_200.txt"
BFBO_PEBO_FST <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fst.50000.Ztransformed.csv"
BRBO_4taxa    <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/brbo_4taxa/BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__5000_200.txt"

convfile <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"
out_pdf  <- "combined_manhattans_outliers.pdf"

color1 <- "#666666"
color2 <- "#000000"
chromo_levels <- c(as.character(1:28), as.character(30:32))
drop_scaffolds <- c("CM062595.1","CM062599.1","CM062600.1","CM062610.1")

# -------------------------------
# Helpers
# -------------------------------
read_conv <- function(convfile) {
  read_csv(convfile, col_names = c("chromo","chr"), show_col_types = FALSE) %>%
    transmute(chr = as.character(chr),
              chromo = as.character(chromo))
}

read_dfm <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)
  req <- c("chr","windowStart","windowEnd")
  stopifnot(all(req %in% names(df)))
  df %>%
    mutate(position = (windowStart + windowEnd)/2,
           source_type = "DFM")
}

read_fst_ready <- function(path) {
  df <- read_tsv(
    path,
    col_types = cols(
      region = col_character(),
      chromo = col_character(),
      position = col_double(),
      Nsites = col_double(),
      fst = col_double(),
      z = col_double(),
      neg_log_pvalues_one_tailed = col_double()
    )
  )
  df %>%
    transmute(
      chromo = as.character(chromo),
      position = position,
      neglog10p = neg_log_pvalues_one_tailed,
      source_type = "FST_READY"
    )
}

prep_for_plot <- function(df, conv, chromo_levels) {
  if ("chromo" %in% names(df)) {
    df2 <- df %>%
      mutate(chromo = as.character(chromo))
  } else {
    df2 <- df %>%
      mutate(chr = as.character(chr)) %>%
      left_join(conv, by = "chr") %>%
      filter(!chr %in% drop_scaffolds)
    if (any(is.na(df2$chromo))) {
      df2 <- df2 %>% filter(!is.na(chromo))
    }
  }

  if (length(intersect(unique(df2$chromo), chromo_levels)) > 0) {
    df2$chromo <- factor(df2$chromo, levels = chromo_levels)
  } else {
    df2$chromo <- as.factor(df2$chromo)
  }

  lens <- df2 %>%
    group_by(chromo) %>%
    summarise(chr_len = max(position, na.rm = TRUE), .groups = "drop") %>%
    arrange(chromo) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(chromo, tot)

  plot_data <- lens %>%
    left_join(df2, by = "chromo") %>%
    arrange(chromo, position) %>%
    mutate(BPcum = position + tot)

  axisdf <- plot_data %>%
    group_by(chromo) %>%
    summarize(center = mean(BPcum), .groups = "drop")

  list(plot_data = plot_data, axisdf = axisdf)
}

plot_manhattan <- function(plot_data, axisdf, metric, title, ylab_override = NULL) {
  stopifnot(metric %in% names(plot_data))

  # determine cutoff and flag outliers
  cutoff <- quantile(plot_data[[metric]], probs = 0.999, na.rm = TRUE)
  plot_data <- plot_data %>% mutate(is_outlier = .data[[metric]] >= cutoff)

  pal <- rep(c(color1, color2), length.out = length(levels(plot_data$chromo)))
  names(pal) <- levels(plot_data$chromo)

  ggplot() +
    geom_point(
      data = subset(plot_data, !is_outlier),
      aes(x = BPcum, y = .data[[metric]], color = chromo),
      alpha = 0.8, size = 0.8, show.legend = FALSE
    ) +
    geom_point(
      data = subset(plot_data, is_outlier),
      aes(x = BPcum, y = .data[[metric]]),
      color = "red3", alpha = 0.95, size = 1.6
    ) +
    scale_color_manual(values = pal) +
    scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Chromosome",
      y = ifelse(is.null(ylab_override), metric, ylab_override),
      title = paste0(title, "  (top 0.001% in red)")

    ) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "red3", linewidth = 0.5) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

# -------------------------------
# Load conversion
# -------------------------------
conv <- read_conv(convfile)

jobs <- tibble::tibble(
  path   = c(MABO_NABO_DFM, MABO_NABO_FST, BFBO_PEBO_DFM, BFBO_PEBO_FST, BRBO_4taxa),
  type   = c("DFM",        "FST_READY",   "DFM",         "FST_READY",   "DFM"),
  metric = c("f_dM",       "neglog10p",   "f_dM",        "neglog10p",   "f_dM"),
  title  = basename(c(MABO_NABO_DFM, MABO_NABO_FST, BFBO_PEBO_DFM, BFBO_PEBO_FST, BRBO_4taxa))
)

# -------------------------------
# Generate multi-page PDF
# -------------------------------
pdf(out_pdf, width = 20, height = 5, onefile = TRUE)

for (i in seq_len(nrow(jobs))) {
  cat("Processing:", jobs$path[i], "(", jobs$type[i], ") -> metric:", jobs$metric[i], "\n")

  df_raw <- switch(
    jobs$type[i],
    "DFM"       = read_dfm(jobs$path[i]),
    "FST_READY" = read_fst_ready(jobs$path[i]),
    stop("Unknown type: ", jobs$type[i])
  )

  pp <- prep_for_plot(df_raw, conv, chromo_levels)

  ylab <- if (jobs$metric[i] == "neglog10p") expression(-log[10](p)) else jobs$metric[i]
  p  <- plot_manhattan(pp$plot_data, pp$axisdf, jobs$metric[i], jobs$title[i], ylab_override = ylab)
  print(p)
}

dev.off()
cat("Saved multi-page PDF with outliers highlighted:", out_pdf, "\n")
