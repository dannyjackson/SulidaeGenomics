# WG heterozygosity

# compare heterozygosity across chr29
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered

cd hetplot_allchr


VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.v42.gz        # bgzipped + indexed VCF

### Per individual heterozygosity
# --- inputs you edit ---
OUT=hets_wg                  # prefix for intermediates

# index if needed
tabix -p vcf "$VCF"

# keep biallelic SNPs on that chromosome (clean het calling)
bcftools view -v snps -m2 -M2 -O z -o ${OUT}.snps.vcf.gz "$VCF"
tabix -p vcf ${OUT}.snps.vcf.gz

# sample order (matches columns in the query below)
bcftools query -l ${OUT}.snps.vcf.gz > ${OUT}.samples.txt

# CHROM POS then one GT per sample (e.g., 0/0, 0/1, 1|0, 1/1, ./.)
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ${OUT}.snps.vcf.gz > ${OUT}.genotypes.tsv

Rscript allchr_heterozygosity.r

# Rscript all_heterozygosity.1.r
#!/usr/bin/env bash
#SBATCH --job-name=all_heterozygosity.2
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=168G
#SBATCH --time=3:00:00
#SBATCH --output=slurm_output/all_heterozygosity.2.out
#SBATCH --mail-type=ALL
# sbatch all_heterozygosity.2.sh

module load micromamba
source ~/.bashrc
micromamba activate r_ocelote

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered/hetplot_allchr/
Rscript /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered/hetplot_allchr/WG_heterozygosity.2.r

# Rscript all_heterozygosity.3.r

sed -i 's/\"//g' 50kb_heterozygosity.labelled.csv

## add sex labels and chromosome numbers
R
library(readr)
library(dplyr)

# Read in the data
het <- read_csv("50kb_heterozygosity.labelled.csv")
chrommap <- read_csv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt",
                     col_names = c("chrom_num", "CHROM"))
sexinfo <- read_csv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv",
                    col_names = c("individual", "sex"))

# Merge in chromosome number and sex
merged <- het %>%
  left_join(chrommap, by = "CHROM") %>%
  left_join(sexinfo, by = "individual") %>%
  select(individual, species, sex, chrom_num, CHROM, win_start, win_mid, n_sites, n_het, het_rate)

# Save
write_csv(merged, "50kb_heterozygosity.labelled_with_chrnum_sex.csv")



########################################
## compute H F:M ratio
########################################


library(readr)
library(dplyr)

# Inputs
het_file   <- "50kb_heterozygosity.labelled_with_chrnum_sex.csv"
sex_file   <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"

# Read
het   <- read_csv(het_file)                                   # has: individual, CHROM, win_start, win_mid, het_rate, ...
sexdf <- read_csv(sex_file, col_names = c("individual","sex")) # "Male"/"Female"


# Assumes `het` has: chrom_num, win_start, win_mid, het_rate, sex, species
# sex values like "Female"/"Male"; species values like "RFBO","BRBO","MABO","NABO","PEBO","BFBO"

ratio_win <- het %>%
  group_by(chrom_num, win_start, win_mid) %>%
  summarise(
    # Baselines across all individuals in the window
    mean_female_all = mean(het_rate[sex == "Female"], na.rm = TRUE),
    mean_male_all   = mean(het_rate[sex == "Male"],   na.rm = TRUE),
    n_female_all    = sum(sex == "Female"),
    n_male_all      = sum(sex == "Male"),

    # Species-specific female means and counts
    mean_female_RFBO = mean(het_rate[sex == "Female" & species == "RFBO"], na.rm = TRUE),
    n_female_RFBO    = sum(sex == "Female" & species == "RFBO"),

    mean_female_BRBO = mean(het_rate[sex == "Female" & species == "BRBO"], na.rm = TRUE),
    n_female_BRBO    = sum(sex == "Female" & species == "BRBO"),

    mean_female_MABO = mean(het_rate[sex == "Female" & species == "MABO"], na.rm = TRUE),
    n_female_MABO    = sum(sex == "Female" & species == "MABO"),

    mean_female_NABO = mean(het_rate[sex == "Female" & species == "NABO"], na.rm = TRUE),
    n_female_NABO    = sum(sex == "Female" & species == "NABO"),

    mean_female_PEBO = mean(het_rate[sex == "Female" & species == "PEBO"], na.rm = TRUE),
    n_female_PEBO    = sum(sex == "Female" & species == "PEBO"),

    mean_female_BFBO = mean(het_rate[sex == "Female" & species == "BFBO"], na.rm = TRUE),
    n_female_BFBO    = sum(sex == "Female" & species == "BFBO"),
    .groups = "drop"
  ) %>%
  mutate(
    # Overall female:male ratio
    ratio_F_over_M = ifelse(
      n_female_all == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_all / mean_male_all
    ),

    # Species-specific female : all-males ratios
    ratio_RFBO = ifelse(
      n_female_RFBO == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_RFBO / mean_male_all
    ),
    ratio_BRBO = ifelse(
      n_female_BRBO == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_BRBO / mean_male_all
    ),
    ratio_MABO = ifelse(
      n_female_MABO == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_MABO / mean_male_all
    ),
    ratio_NABO = ifelse(
      n_female_NABO == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_NABO / mean_male_all
    ),
    ratio_PEBO = ifelse(
      n_female_PEBO == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_PEBO / mean_male_all
    ),
    ratio_BFBO = ifelse(
      n_female_BFBO == 0 | n_male_all == 0 | is.na(mean_male_all) | mean_male_all == 0,
      NA_real_, mean_female_BFBO / mean_male_all
    )
  ) %>%
  arrange(chrom_num, win_start, win_mid)


# Save
out_file <- "window_het_ratio_F_over_M.csv"
write_csv(ratio_win, out_file)
message("Wrote: ", out_file)





###################################################
# Make manhattan plot
###################################################
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read in the ratio file
ratio <- read_csv("window_het_ratio_F_over_M.csv")
ratio <- ratio %>%
  rename(CHROM = chrom_num)
# Optional: drop windows with missing ratios
ratio <- ratio %>% filter(!is.na(ratio_F_over_M))

# Order chromosomes numerically
ratio <- ratio %>%
  mutate(
    CHROM = factor(CHROM, levels = unique(CHROM[order(as.numeric(gsub("\\D", "", CHROM)))]))
  )

# Compute chromosome offsets for cumulative plotting
chr_offsets <- ratio %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(win_mid, na.rm = TRUE)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

ratio <- ratio %>%
  left_join(chr_offsets, by = "CHROM") %>%
  arrange(CHROM, win_mid) %>%
  mutate(BPcum = win_mid + tot)

# X-axis labels: chromosome midpoints
axis_df <- ratio %>%
  group_by(CHROM) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)



# Extract chromosome numbers
# Assign chromosome categories for coloring
ratio <- ratio %>%
  mutate(
    chrom_num = as.numeric(gsub("\\D", "", CHROM)),   # extract numeric part
    color_group = case_when(
      grepl("Z", CHROM, ignore.case = TRUE) ~ "Z",
      chrom_num == 29 ~ "29",
      TRUE ~ "other"
    )
  )


# Manhattan plot
p <- ggplot(ratio, aes(x = BPcum, y = ratio_F_over_M, color = color_group)) +
  geom_point(size = 0.1, alpha = 0.8) +
  scale_color_manual(values = c("other" = "gray70", "29" = "#2894FF", "Z" = "#D73027"), guide = "none") +
  scale_x_continuous(
    label = axis_df$CHROM,
    breaks = axis_df$center,
    expand = c(0.01, 0.01)
  ) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "Chromosome",
    y = "H F:M"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )


ggsave("window_het_ratio_FoverM_manhattan_highlightZ29.pdf",
       p, width = 10, height = 2, dpi = 300)



# Read in the ratio file
ratio <- read_csv("window_het_ratio_F_over_M.csv")
ratio <- ratio %>%
  rename(CHROM = chrom_num)
# Optional: drop windows with missing ratios
ratio <- ratio %>% filter(!is.na(ratio_F_over_M))

ratio_filt <- ratio %>%
  filter(CHROM %in% c("28", "29", "30", "31", "32", "33"))

# Order chromosomes numerically
ratio_filt <- ratio_filt %>%
  mutate(
    CHROM = factor(CHROM, levels = unique(CHROM[order(as.numeric(gsub("\\D", "", CHROM)))]))
  )

# Compute chromosome offsets for cumulative plotting
chr_offsets <- ratio_filt %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(win_mid, na.rm = TRUE)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

ratio_filt <- ratio_filt %>%
  left_join(chr_offsets, by = "CHROM") %>%
  arrange(CHROM, win_mid) %>%
  mutate(BPcum = win_mid + tot)

# X-axis labels: chromosome midpoints
axis_df <- ratio_filt %>%
  group_by(CHROM) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)


# Assign chromosome categories for coloring
ratio_filt <- ratio_filt %>%
  mutate(
    chrom_num = as.numeric(gsub("\\D", "", CHROM)),   # extract numeric part
    color_group = case_when(
      grepl("Z", CHROM, ignore.case = TRUE) ~ "Z",
      chrom_num == 29 ~ "29",
      TRUE ~ "other"
    )
  )

# Manhattan plot

p <- ggplot(ratio_filt, aes(x = BPcum, y = ratio_F_over_M, color = color_group)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("other" = "gray70", "29" = "#2894FF", "Z" = "#D73027"), guide = "none") +
  scale_x_continuous(
    label = axis_df$CHROM,
    breaks = axis_df$center,
    expand = c(0.01, 0.01)
  ) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.4) +
  labs(
    x = NULL,
    y = NULL,
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )
  

ggsave("window_het_ratio_FoverM_manhattan_autosome_vs_Z29.pdf",
       p, width = 5, height = 2, dpi = 300)


# do the same for Z

ratio_filt <- ratio %>%
  filter(CHROM %in% c("Z"))

# Order chromosomes numerically
ratio_filt <- ratio_filt %>%
  mutate(
    CHROM = factor(CHROM, levels = unique(CHROM[order(as.numeric(gsub("\\D", "", CHROM)))]))
  )

# Compute chromosome offsets for cumulative plotting
chr_offsets <- ratio_filt %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(win_mid, na.rm = TRUE)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

ratio_filt <- ratio_filt %>%
  left_join(chr_offsets, by = "CHROM") %>%
  arrange(CHROM, win_mid) %>%
  mutate(BPcum = win_mid + tot)

# X-axis labels: chromosome midpoints
axis_df <- ratio_filt %>%
  group_by(CHROM) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)


# Assign chromosome categories for coloring
ratio_filt <- ratio_filt %>%
  mutate(
    chrom_num = as.numeric(gsub("\\D", "", CHROM)),   # extract numeric part
    color_group = case_when(
      grepl("Z", CHROM, ignore.case = TRUE) ~ "Z",
      chrom_num == 29 ~ "29",
      TRUE ~ "other"
    )
  )

# Manhattan plot

p <- ggplot(ratio_filt, aes(x = BPcum, y = ratio_F_over_M, color = color_group)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c("other" = "gray70", "29" = "#2894FF", "Z" = "#D73027"), guide = "none") +
  scale_x_continuous(
    label = axis_df$CHROM,
    breaks = axis_df$center,
    expand = c(0.01, 0.01)
  ) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.4) +
  labs(
    x = NULL,
    y = NULL,
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )
  

ggsave("window_het_ratio_FoverM_manhattan_Z.pdf",
       p, width = 5, height = 2, dpi = 300)
















# plot as lines 
p <- ggplot(ratio_filt, aes(x = BPcum, y = ratio_F_over_M, color = color_group, group = chrom_num)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  scale_color_manual(
    values = c("other" = "gray70", "29" = "#2894FF", "Z" = "#D73027"),
    guide = "none"
  ) +
  scale_x_continuous(
    labels = axis_df$CHROM,
    breaks = axis_df$center,
    expand = c(0.01, 0.01)
  ) +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "Chromosome",
    y = expression(Female / Male ~ heterozygosity),
    title = "Female-to-Male Heterozygosity Ratio (50 kb windows)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ggsave("window_het_ratio_FoverM_lines_highlightZ29.png",
       p, width = 10, height = 2, dpi = 300)


### Plot with species specific ratios


# 1) Read color codes for species
cc <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/colorcodes.tsv", show_col_types = FALSE)
stopifnot(all(c("Species", "Code") %in% names(cc)))
cols <- setNames(cc$Code, cc$Species)

# 2) Pivot species ratio columns long
ratio_long <- ratio_filt %>%
  select(chrom_num, BPcum,
         ratio_RFBO, ratio_BRBO, ratio_MABO, ratio_NABO, ratio_PEBO, ratio_BFBO, ratio_F_over_M) %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "Species",
    values_to = "ratio"
  ) %>%
  mutate(Species = toupper(sub("^ratio_", "", Species)))

# Keep only colors for present species
cols_use <- cols[intersect(names(cols), unique(ratio_long$Species))]

# 3) Plot: line-based Manhattan-style plot per species

# 3) Plot: lines per species + bold aggregate line on top
p_species_lines <- ggplot() +
  # species lines first
  geom_line(
    data = ratio_long %>% filter(Species != "F_OVER_M"),
    aes(x = BPcum, y = ratio, color = Species, group = interaction(Species, chrom_num)),
    linewidth = 0.3, alpha = 0.5
  ) +
  scale_color_manual(values = cols_use, name = "Species") +

  # aggregate Female/Male ratio line overtop with chr-specific coloring
  # overall Female/Male line: chr 28 & 30 (light gray)
  geom_line(
    data = ratio_filt %>%
      filter(as.character(chrom_num) %in% c("28", "chr28", "30", "chr30")),
    aes(x = BPcum, y = ratio_F_over_M, group = chrom_num),
    color = "gray25", linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +
  # overall Female/Male line: chr 29 (dark gray)
  geom_line(
    data = ratio_filt %>%
      filter(as.character(chrom_num) %in% c("29", "chr29")),
    aes(x = BPcum, y = ratio_F_over_M, group = chrom_num),
    color = "lightgray", linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +

  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(
    labels = axis_df$CHROM,
    breaks = axis_df$center,
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Chromosome",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

ggsave("species_ratio_lines_manhattan.png",
       p_species_lines, width = 10, height = 2, dpi = 300)


p_species_lines <- ggplot() +

  # overall Female/Male line: chr 28 & 30 (light gray)
  geom_line(
    data = ratio_filt %>%
      filter(as.character(chrom_num) %in% c("28", "chr28", "30", "chr30")),
    aes(x = BPcum, y = ratio_F_over_M, group = chrom_num),
    color = "gray25", linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +
  # overall Female/Male line: chr 29 (dark gray)
  geom_line(
    data = ratio_filt %>%
      filter(as.character(chrom_num) %in% c("29", "chr29")),
    aes(x = BPcum, y = ratio_F_over_M, group = chrom_num),
    color = "lightgray", linewidth = 1, alpha = 0.9, inherit.aes = FALSE
  ) +

  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(
    labels = axis_df$CHROM,
    breaks = axis_df$center,
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Chromosome",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none"
  )

ggsave("species_ratio_lines_chr29.png",
       p_species_lines, width = 10, height = 3, dpi = 300)
