# --- packages ---
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

WIN <- 50000
OUT <- "hets_CM062595.1"

# -------------------------
# Read sample names
# -------------------------
samples <- sub("\\.bam$", "",
               basename(readLines(paste0(OUT, ".samples.txt"))))

colnames_geno <- c("CHROM", "POS", samples)

geno <- read_tsv(
  paste0(OUT, ".genotypes.tsv"),
  col_names = colnames_geno,
  col_types = cols(
    CHROM = col_character(),
    POS   = col_double(),
    .default = col_character()
  ),
  progress = FALSE
)

# -------------------------
# Long format
# -------------------------
long <- geno |>
  pivot_longer(
    cols = all_of(samples),
    names_to = "individual",
    values_to = "GT"
  )

# -------------------------
# Heterozygote flag
# -------------------------
is_het <- function(gt) {
  ifelse(is.na(gt) | gt %in% c("./.", ".|."), NA,
         {
           m <- str_match(gt, "^([0-9.]+)[/|]([0-9.]+)$")
           ifelse(is.na(m[,1]), NA,
                  (m[,2] != m[,3]) & (m[,2] != ".") & (m[,3] != "."))
         })
}

# -------------------------
# Windowed heterozygosity
# -------------------------
hets <- long |>
  mutate(
    het = is_het(GT),
    win_start = floor((POS - 1) / WIN) * WIN,
    win_mid   = win_start + WIN / 2
  ) |>
  group_by(individual, CHROM, win_start, win_mid) |>
  summarise(
    n_sites  = sum(!is.na(het)),
    n_het    = sum(het, na.rm = TRUE),
    het_rate = n_het / n_sites,
    .groups = "drop"
  )

# -------------------------
# Join sex information
# -------------------------
sex <- read_csv(
  "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv",
  col_names = c("individual", "sex")
)

hets_sex <- hets |>
  left_join(sex, by = "individual") |>
  mutate(
    sex = str_squish(sex),
    sex = coalesce(sex, "Unknown")
  )

# -------------------------
# Per-window female:male ratio + test
# -------------------------
hets_ratio <- hets_sex |>
  group_by(CHROM, win_start, win_mid) |>
  summarise(
    n_sites = mean(n_sites, na.rm = TRUE),
    mean_female = mean(het_rate[sex == "Female"], na.rm = TRUE),
    mean_male   = mean(het_rate[sex == "Male"],   na.rm = TRUE),
    ratio = mean_female / mean_male,
    pval = tryCatch(
      t.test(het_rate[sex == "Female"],
             het_rate[sex == "Male"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) |>
  mutate(
    padj = p.adjust(pval, method = "BH")
  ) |>
  select(CHROM, win_start, win_mid, n_sites, ratio, pval, padj)

# -------------------------
# Collapse consecutive runs
# -------------------------
padj_cutoff <- 0.001
min_run <- 3

collapse_runs <- function(df, is_sig = TRUE) {
  df |>
    arrange(CHROM, win_start) |>
    group_by(CHROM) |>
    mutate(
      flag = if (is_sig) padj < padj_cutoff else padj >= padj_cutoff,
      run_id = cumsum(flag != lag(flag, default = first(flag)))
    ) |>
    filter(flag) |>
    group_by(CHROM, run_id) |>
    summarise(
      n_windows   = n(),
      start       = min(win_start),
      end         = max(win_start) + WIN,
      mid_start   = min(win_mid),
      mid_end     = max(win_mid),
      min_padj    = min(padj, na.rm = TRUE),
      min_pval    = min(pval, na.rm = TRUE),
      total_sites = sum(n_sites, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(n_windows >= min_run) |>
    arrange(CHROM, start)
}

sig_regions    <- collapse_runs(hets_ratio, is_sig = TRUE)
nonsig_regions <- collapse_runs(hets_ratio, is_sig = FALSE)

# -------------------------
# Write outputs
# -------------------------
write_tsv(sig_regions,    "regions.sig_runs.NR.tsv")
write_tsv(nonsig_regions, "regions.nonsig_runs.PAR.tsv")
