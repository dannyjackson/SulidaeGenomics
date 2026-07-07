#!/usr/bin/env Rscript

# ============================================================
# Combined Manhattan plot (all datasets on same x-axis)
# Output: single-page PDF with shared chromosome axis
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(rlang)
  library(stringr)
  library(tidyr)
  library(patchwork)
})

# -------------------------------
# File paths
# -------------------------------
MABO_NABO_DFM <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo/MABO_AtlCar_MABO_IndoPacific_NABO_localFstats__5000_200.txt"
MABO_NABO_FST <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fst.50000.Ztransformed.csv"
BFBO_PEBO_DFM <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/bfbo_pebo/BFBO_GofCA_BFBO_southern_PEBO_localFstats__5000_200.txt"
BFBO_PEBO_FST <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fst.50000.Ztransformed.csv"
BRBO_4taxa    <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/brbo_4taxa/BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__5000_200.txt"

convfile <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"

out_pdf  <- "combined_manhattans_oneaxis.pdf"
out_png  <- "combined_manhattans_oneaxis.png"

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

# FST reader 

read_fst_ready <- function(path) {
  df <- read_tsv(
    path, col_types = cols( region = col_character(), 
    chromo = col_character(), position = col_double(), 
    Nsites = col_double(), 
    fst = col_double(), 
    z = col_double(), 
    neg_log_pvalues_one_tailed = col_double() ) 
  )
  df %>% 
    transmute(
      chromo = as.character(chromo),
      position = position,
      neglog10p = neg_log_pvalues_one_tailed,
      source_type = "FST_READY" ) }

prep_for_plot <- function(df, conv, chromo_levels) {
  if ("chromo" %in% names(df)) {
    df2 <- df
  } else {
    df2 <- df %>%
      mutate(chr = as.character(chr)) %>%
      left_join(conv, by = "chr") %>%
      filter(!chr %in% drop_scaffolds)
  }

  df2$chromo <- factor(df2$chromo, levels = chromo_levels)

  lens <- df2 %>%
    group_by(chromo) %>%
    summarise(chr_len = max(position, na.rm = TRUE), .groups = "drop") %>%
    arrange(chromo) %>%
    mutate(tot = cumsum(chr_len) - chr_len)

  plot_data <- df2 %>%
    left_join(lens, by = "chromo") %>%
    mutate(BPcum = position + tot)

  axisdf <- lens %>%
    mutate(center = tot + chr_len/2)

  list(plot_data = plot_data, axisdf = axisdf)
}

make_panel <- function(plot_data, axisdf, metric, title, ylab, show_x = FALSE) {
  cutoff <- quantile(plot_data[[metric]], 0.999, na.rm = TRUE)
  plot_data <- plot_data %>%
    mutate(is_outlier = .data[[metric]] >= cutoff)

  pal <- rep(c(color1, color2), length.out = length(levels(plot_data$chromo)))
  names(pal) <- levels(plot_data$chromo)

  p <- ggplot() +
    geom_point(data = subset(plot_data, !is_outlier),
               aes(x = BPcum, y = .data[[metric]], color = chromo),
               alpha = 0.8, size = 0.6) +
    geom_point(data = subset(plot_data, is_outlier),
               aes(x = BPcum, y = .data[[metric]]),
               color = "red3", size = 1.3) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "red3", linewidth = 0.4) +
    scale_color_manual(values = pal) +
    scale_x_continuous(labels = if (show_x) axisdf$chromo else NULL,
                       breaks = axisdf$center,
                       guide = guide_axis(n.dodge = 2)) +
    labs(x = NULL, y = ylab, title = title) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0, face = "plain"),  # left-aligned
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank()
    )
  p
}

standardize_for_axis <- function(df, conv, chromo_levels, drop_scaffolds) {
  if (!("chromo" %in% names(df))) {
    df <- df %>%
      mutate(chr = as.character(chr)) %>%
      left_join(conv, by = "chr") %>%
      filter(!chr %in% drop_scaffolds)
  }

  df %>%
    mutate(
      chromo = as.character(chromo),
      chromo = factor(chromo, levels = chromo_levels)
    ) %>%
    filter(!is.na(chromo))  # drop anything not in chromo_levels
}

build_global_axis <- function(all_dfs) {
  lens <- bind_rows(all_dfs) %>%
    group_by(chromo) %>%
    summarise(chr_len = max(position, na.rm = TRUE), .groups = "drop") %>%
    arrange(chromo) %>%
    mutate(tot = cumsum(chr_len) - chr_len)

  axisdf <- lens %>% mutate(center = tot + chr_len / 2)

  list(lens = lens, axisdf = axisdf, x_max = max(lens$tot + lens$chr_len, na.rm = TRUE))
}

add_bp_cum <- function(df, lens_global) {
  df %>%
    left_join(lens_global %>% select(chromo, tot), by = "chromo") %>%
    mutate(BPcum = position + tot)
}

# -------------------------------
# Combine data and plot
# -------------------------------
conv <- read_conv(convfile)

jobs <- tibble::tibble(
  path   = c(MABO_NABO_DFM, MABO_NABO_FST, BFBO_PEBO_DFM, BFBO_PEBO_FST, BRBO_4taxa),
  type   = c("DFM",        "FST_READY",   "DFM",         "FST_READY",   "DFM"),
  metric = c("f_dM",       "neglog10p",   "f_dM",        "neglog10p",   "f_dM")
)

# Custom titles (in the same order as 'jobs')
titles <- list(
  "Introgression between IndoPacific Masked and Nazca Boobies",
  expression(F[ST]*" Masked vs. Nazca Booby"),
  "Introgression between Blue-Footed and Peruvian Boobies",
  expression(F[ST]*" Blue-Footed vs. Peruvian Booby"),
  "Introgression between Atlantic/Caribbean Brown Boobies and the 4-taxa Clade"
)


# 1) Read + standardize all datasets first
raw_list <- vector("list", nrow(jobs))

for (i in seq_len(nrow(jobs))) {
  cat("Reading:", jobs$path[i], "\n")
  df_raw <- switch(
    jobs$type[i],
    "DFM"       = read_dfm(jobs$path[i]),
    "FST_READY" = read_fst_ready(jobs$path[i])
  )

  raw_list[[i]] <- standardize_for_axis(df_raw, conv, chromo_levels, drop_scaffolds)
}

# 2) Build ONE global axis from all datasets combined
axis_global <- build_global_axis(raw_list)
shared_axis <- axis_global$axisdf
lens_global <- axis_global$lens
x_max <- axis_global$x_max

# 3) Build plots using the SAME BPcum mapping + same x limits
plots <- list()
n_panels <- nrow(jobs)

for (i in seq_len(n_panels)) {
  df_i <- add_bp_cum(raw_list[[i]], lens_global)

  ylab <- if (jobs$metric[i] == "neglog10p") expression(-log[10](p)) else jobs$metric[i]
  show_x <- (i == n_panels)

  plots[[i]] <-
    make_panel(df_i, shared_axis, jobs$metric[i], titles[[i]], ylab, show_x = show_x) +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = shared_axis$center,
      labels = if (show_x) shared_axis$chromo else NULL,
      guide = guide_axis(n.dodge = 2)
    ) +
    labs(x = if (show_x) "Chromosome" else NULL)
}

combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.margin = margin(10, 10, 20, 10),
      plot.tag = element_text(face = "bold")
    )
  )

ggsave(out_png, combined_plot, width = 20, height = 16)
cat("Saved single-page combined Manhattan plot:", out_png, "\n")

ggsave(out_pdf, combined_plot, width = 20, height = 16)
cat("Saved single-page combined Manhattan plot:", out_pdf, "\n")


