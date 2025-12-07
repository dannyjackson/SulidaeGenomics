#!/usr/bin/env Rscript

# ============================================================
# Plot log-likelihood vs K and compute Evanno ΔK
# Works with:
#  - K, logL  (one run per K), or
#  - K, replicate, logL (multiple runs per K)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# -------------------------------
# Input file
# -------------------------------
infile <- "admix_likelihoods.tsv"

# Expect at least: K, logL
lik <- read_tsv(infile, show_col_types = FALSE)

if (!all(c("K", "logL") %in% names(lik))) {
  stop("Input file must have at least columns: K and logL")
}

# Make sure K is numeric and ordered
lik <- lik %>%
  mutate(K = as.numeric(K)) %>%
  arrange(K)

# ------------------------------------------------------------
# 1) If no replicate column: treat as one run per K
#    If replicate exists: compute mean & SD per K (Evanno style)
# ------------------------------------------------------------
if ("replicate" %in% names(lik)) {
  message("Detected 'replicate' column: using replicate-based Evanno method.")

  summary_lik <- lik %>%
    group_by(K) %>%
    summarise(
      mean_logL = mean(logL),
      sd_logL   = sd(logL),
      n         = n(),
      .groups   = "drop"
    ) %>%
    arrange(K)

  # For later calculations, rename to common columns
  summary_lik <- summary_lik %>%
    rename(logL = mean_logL)

} else {
  message("No 'replicate' column detected: treating as one run per K.")
  
  # Fake SD = 1 for ΔK denominator to avoid division by zero
  summary_lik <- lik %>%
    mutate(sd_logL = 1,
           n = 1) %>%
    arrange(K)
}

# ------------------------------------------------------------
# 2) Compute Evanno-style stats
#    L'(K) = L(K) - L(K-1)
#    L''(K) = L'(K+1) - L'(K)
#    ΔK = |L''(K)| / SD(L(K))
# ------------------------------------------------------------
evanno <- summary_lik %>%
  arrange(K) %>%
  mutate(
    L_prime       = logL - lag(logL),
    L_prime_next  = lead(logL) - logL,
    L_doubleprime = L_prime_next - L_prime,
    DeltaK        = abs(L_doubleprime) / sd_logL
  )

# Write out table for inspection
out_evanno <- "admix_evanno_stats.tsv"
write_tsv(evanno, out_evanno)
message("Wrote Evanno stats to: ", out_evanno)

# ------------------------------------------------------------
# 3) Plots
# ------------------------------------------------------------

## 3a. Log-likelihood vs K
p_logL <- ggplot(summary_lik, aes(x = K, y = logL)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "K (number of clusters)",
    y = "Log-likelihood",
    title = "NGSadmix log-likelihood vs K"
  )

ggsave(
  filename = "admix_logL_vs_K.pdf",
  plot     = p_logL,
  width    = 6,
  height   = 4
)

## 3b. ΔK vs K (Evanno)
# ΔK is defined for K between min+1 and max-1 (because of lag/lead)
evanno_for_plot <- evanno %>%
  filter(!is.na(DeltaK))

p_DeltaK <- ggplot(evanno_for_plot, aes(x = K, y = DeltaK)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw(base_size = 14) +
  labs(
    x = "K",
    y = expression(Delta*K),
    title = "Evanno ΔK vs K"
  )

ggsave(
  filename = "admix_DeltaK_vs_K.pdf",
  plot     = p_DeltaK,
  width    = 6,
  height   = 4
)

message("Saved plots")