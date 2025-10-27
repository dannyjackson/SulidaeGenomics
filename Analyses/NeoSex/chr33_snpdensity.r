cd /xdisk/mcnew/dannyjackson/sulidae/analyses/snpdensity/chr33
# visualize the density of snps
VCF="/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062599.1.qualitysort_filtered_mind2.vcf"

# 10 kb non-overlapping windows
while read -r IND; do
vcftools --gzvcf $VCF --SNPdensity 10000 --indv ${IND} --out snpden.${IND} --non-ref-ac 1 --mac 1 --max-mac 2
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

# merge into one big matrix
#!/usr/bin/env bash

# inputs
LIST="/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt"
OUT="merged.snpden.tsv"

# read individuals from list
mapfile -t INDS < "$LIST"

# collect existing snpden files
files=()
for IND in "${INDS[@]}"; do
  f="snpden.${IND}.snpden"
  if [[ -s "$f" ]]; then
    files+=("$f")
  else
    echo "Warning: missing or empty file: $f" >&2
  fi
done

if [[ ${#files[@]} -eq 0 ]]; then
  echo "No snpden files found. Exiting." >&2
  exit 1
fi

# merge with AWK, then sort by CHROM and BIN_START
awk -v OFS='\t' -v header="$(printf "%s " "${INDS[@]}")" '
BEGIN {
  n = split(header, H, " ");
  # print header line
  printf "CHROM\tBIN_START";
  for (i = 1; i <= n; i++) printf "\t%s", H[i];
  print "";
}
FNR == 1 { next }  # skip per-file header
{
  # deduce IND from filename: snpden.<IND>.snpden
  fn = FILENAME;
  sub(/.*snpden\./, "", fn);
  sub(/\.snpden.*/, "", fn);

  key = $1 OFS $2;      # CHROM \t BIN_START
  keys[key] = 1;
  chrom[key] = $1;
  bin[key]   = $2;

  # third column is SNP_COUNT
  cnt[key, fn] = $3;
}
END {
  for (k in keys) {
    printf "%s\t%s", chrom[k], bin[k];
    for (i = 1; i <= n; i++) {
      ind = H[i];
      v = cnt[k, ind];
      if (v == "") v = 0;  # fill missing with 0
      printf "\t%s", v;
    }
    print "";
  }
}
' "${files[@]}" | sort -k1,1 -k2,2n > "$OUT"

echo "Wrote: $OUT"

# plot it
# ---- packages ----
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
  labs(x = "Position (bp)", y = "SNPs per 10 kb", color = "Sex") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "null"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))

print(p)

# ---- save (optional) ----
ggsave("snp_density_by_sex.faceted.chr33.png", p, width = 4, height = 3, dpi = 300)

