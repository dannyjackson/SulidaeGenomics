library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ---- inputs ----
merged_path <- "merged.snpden.tsv"
sex_path    <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"

# ---- read data ----
wide <- read_tsv(merged_path, show_col_types = FALSE)    # columns: CHROM, BIN_START, BFBO501, BFBO502, ...
sex <- read_csv(
  "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv",
  col_names = c("individual", "sex"),
  show_col_types = FALSE
)

# standardize sex table column names if needed
# (tries to find a column that matches the individual IDs, then rename to 'individual')
if (!"individual" %in% names(sex)) {
  cand <- intersect(names(sex), c("individual","ID","id","sample","Sample","SampleID","sample_id","sampleID"))
  if (length(cand) == 0) stop("Could not find an ID column in sex table; expected something like 'individual' or 'SampleID'.")
  sex <- sex %>% rename(individual = !!sym(cand[1]))
}
if (!"sex" %in% names(sex)) {
  cand <- intersect(names(sex), c("sex","Sex","SEX"))
  if (length(cand) == 0) stop("Could not find a 'sex' column in sex table.")
  sex <- sex %>% rename(sex = !!sym(cand[1]))
}

# Ensure ID formatting matches column names in merged file (e.g., drop any suffixes in sex file if present)
# Example: convert any ".idxstats" -> "" in sex IDs, trim spaces
sex <- sex %>%
  mutate(
    individual = str_replace(individual, "\\.idxstats$", ""),
    individual = str_trim(individual)
  )

# ---- pivot to long format ----
long <- wide %>%
  pivot_longer(
    cols = -c(CHROM, BIN_START),
    names_to = "individual",
    values_to = "SNP_COUNT"
  ) %>%
  # if your merged used "NA" strings or blanks, coerce to numeric
  mutate(SNP_COUNT = suppressWarnings(as.numeric(SNP_COUNT)))

# ---- join sex info ----
dat <- long %>%
  left_join(sex, by = "individual")

# Optional: set Unknown if missing
dat <- dat %>%
  mutate(sex = ifelse(is.na(sex) | sex == "", "Unknown", sex))

# ---- plotting ----
p <- ggplot(dat, aes(x = BIN_START, y = SNP_COUNT,
                     group = individual, color = sex)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_color_manual(values = c(Female = "#E26D5A", Male = "#4F7CAC", Unknown = "grey50")) +
  labs(x = "Position (bp)", y = "SNPs per 100 kb", color = "Sex") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "null"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))

print(p)

# ---- save (optional) ----
ggsave("snp_density_by_sex.faceted.chr29.png", p, width = 12, height = 3, dpi = 300)

