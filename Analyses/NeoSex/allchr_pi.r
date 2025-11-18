cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered

mkdir hetplot_allchr
cd hetplot_allchr

module load vcftools samtools htslib
# compare heterozygosity of all chromosomes
IND=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

while read -r CHROM; do
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/"${CHROM}".qualitysort_filtered_mind2.vcf        # bgzipped + indexed VCF
# gunzip "${VCF}.gz"
# sed -i 's/VCFv4\.3/VCFv4\.1/' ${VCF}
# bgzip ${VCF}
# tabix "${VCF}.gz"
vcftools --gzvcf "$VCF.gz" --chr "${CHROM}" \
         --keep ${IND} \
         --site-pi \
         --out "${CHROM}".pi
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt

while read -r CHROM; do
awk '{sum += $3} END {print "Mean π =", sum/NR}' "${CHROM}".pi.sites.pi > "${CHROM}".pi.avg
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt


#!/usr/bin/env bash
set -euo pipefail

INDFILE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
CONTIGS=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
VCF_DIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms


while read -r IND; do
while read -r CHROM; do
  VCF="${VCF_DIR}/${CHROM}.qualitysort_filtered_mind2.vcf.gz"

  # per-site π for this chromosome, restricted to IND list (population π)
  vcftools --gzvcf "${VCF}" --chr "${CHROM}" \
           --indv "${IND}" \
           --site-pi \
           --out "${CHROM}.${IND}.pi"
done < "${CONTIGS}"
done < "${INDFILE}"

#!/usr/bin/env bash
set -euo pipefail

# Directory containing the *.pi.sites.pi files
PI_DIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered/hetplot_allchr/pifiles
OUT=mean_pi_matrix.tsv

# Write header
{
  printf "CHROM"
  while read -r IND; do printf "\t%s" "$IND"; done < ${INDFILE}
  printf "\n"
} > "$OUT"

# For each chromosome, compute mean π per individual
while read -r CHROM; do
  printf "%s" "$CHROM" >> "$OUT"
  while read -r IND; do
    file="${PI_DIR}/${CHROM}.${IND}.pi.sites.pi"
    if [[ -f "$file" ]]; then
      mean=$(awk 'NR>1 && $3 ~ /^[0-9.eE+-]+$/ {sum+=$3; n++} END{if(n>0) printf("%.10g", sum/n); else printf("NA")}' "$file")
    else
      mean="NA"
    fi
    printf "\t%s" "$mean" >> "$OUT"
  done < ${INDFILE}
  printf "\n" >> "$OUT"
done < ${CONTIGS}

echo "Wrote mean π matrix to: $OUT"


# module purge 
module load micromamba
# micromamba create -n r_puma r-base r-tidyr -c conda-forge

micromamba create -n r_puma \
  -c https://conda.anaconda.org/conda-forge \
  -c https://conda.anaconda.org/bioconda \
  r-base r-tidyr

micromamba activate r_puma


# Load packages
library(tidyverse)

# --- Input files ---
pi_file   <- "mean_pi_matrix.tsv"
sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"
conv_file <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"

# --- Read data ---
pi  <- read_tsv(pi_file, na = c("", "NA", ".", "NaN"))
sex <- read_csv(sex_file, col_names = c("individual", "sex"))
conv <- read_csv(conv_file, col_names = c("chrom_num", "CHROM"))

# --- Reshape & join ---
pi_sex <- pi %>%
  pivot_longer(-CHROM, names_to = "individual", values_to = "pi", values_drop_na = TRUE) %>%
  left_join(sex, by = "individual") %>%
  left_join(conv, by = "CHROM") %>%       # attach chromosome number
  mutate(sex = replace_na(sex, "Unknown")) %>%
  filter(is.finite(pi))

# Ensure chrom is character
pi_sex <- pi_sex %>% mutate(chrom_num = as.character(chrom_num))


# Compute desired factor levels once
lev_base   <- str_sort(unique(pi_sex$chrom_num), numeric = TRUE)
specials   <- c("Z","W","X","MT","M","Un","Unplaced")
lev_order  <- c(setdiff(lev_base, specials),
                specials[specials %in% lev_base])   # push specials to the end (if present)


# Apply levels
pi_sex <- pi_sex %>%
  mutate(
    chrom_num = factor(chrom_num, levels = lev_order),
    sex = factor(sex, levels = c("Female", "Male", "Unknown"))
  )
  

# --- Plot: chromosomes on x-axis, violins by sex ---
pd <- position_dodge(width = 0.8)

gg <- ggplot(pi_sex, aes(x = chrom_num, y = pi, fill = sex)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Heterozygosity (", pi, ")"))
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  )


# --- Save ---
ggsave("heterozygosity_by_chromosome.sex.pdf", gg, width = 12, height = 3, dpi = 300)


pi_sex <- pi_sex %>%
  mutate(species = str_extract(individual, "^[A-Z]{4}"))


keep_chroms <- c("1", "2", "3", "29", "33", "Z")

pi_species <- pi_sex %>%
  filter(chrom_num %in% keep_chroms) %>%
  mutate(chrom_num = factor(chrom_num, levels = keep_chroms))

gg <- ggplot(pi_species, aes(x = chrom_num, y = pi, fill = species)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(BFBO="#4FE98C", PEBO="#2894FF", RFBO="#F3447D", BRBO="#000000", MABO="#E8CA3D", NABO="#B32904")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Heterozygosity (", pi, ")"))
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  )


# --- Save ---
ggsave("heterozygosity_by_chromosome.species.pdf", gg, width = 6, height = 3, dpi = 300)

