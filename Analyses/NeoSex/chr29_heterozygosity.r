# compare heterozygosity across chr29
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered

mkdir hetplot_chr29
cd hetplot_chr29

# inputs
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/snpdensity/CM062595.1.qualitysort_filtered_mind2.vcf        # bgzipped + indexed VCF
# VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf
CHR=CM062595.1       # chromosome/scaffold name
WIN=1000000            # window size
STEP=1000000           # step size (non-overlapping)

# make sure VCF is indexed
# bgzip "$VCF"
# tabix -p vcf "$VCF.gz"

# compute π in windows
vcftools --gzvcf "$VCF.gz" --chr "$CHR" \
         --window-pi $WIN --window-pi-step $STEP \
         --out ${CHR}.pi_${WIN}

# Output: ${CHR}.pi_${WIN}.windowed.pi
# columns: CHROM  BIN_START  BIN_END  N_SNPs  PI

library(readr); library(dplyr); library(ggplot2)

win <- read_tsv("CM062567.1.pi_1000000.windowed.pi",
                show_col_types = FALSE,
                col_types = cols(
                  CHROM = col_character(),
                  BIN_START = col_double(),
                  BIN_END = col_double(),
                  N_VARIANTS = col_double(),
                  PI = col_double()
                ))

ggplot(win, aes(x = (BIN_START + BIN_END)/2, y = PI)) +
  geom_line(linewidth = 0.4) +
  labs(title = "Windowed nucleotide diversity (π)",
       x = "Position (bp)", y = "π per 50 kb") +
  theme_bw()


ggsave("pi.CM062595.1.allindv.png", width = 12, height = 24, dpi = 300)


### Per individual heterozygosity
# --- inputs you edit ---
CHR=CM062595.1                    # chromosome/scaffold to plot (change as needed)
OUT=hets_${CHR}                   # prefix for intermediates
VCF="$VCF.gz"
# index if needed
tabix -p vcf "$VCF"

# keep biallelic SNPs on that chromosome (clean het calling)
bcftools view -r "$CHR" -v snps -m2 -M2 -O z -o ${OUT}.snps.vcf.gz "$VCF"
tabix -p vcf ${OUT}.snps.vcf.gz

# sample order (matches columns in the query below)
bcftools query -l ${OUT}.snps.vcf.gz > ${OUT}.samples.txt

# CHROM POS then one GT per sample (e.g., 0/0, 0/1, 1|0, 1/1, ./.)
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ${OUT}.snps.vcf.gz > ${OUT}.genotypes.tsv

R
# --- packages ---
library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2)

WIN <- 50000
OUT <- "hets_CM062595.1"   # must match the bash OUT above (or set dynamically)

samples <- readLines(paste0(OUT, ".samples.txt"))
colnames_geno <- c("CHROM","POS", samples)

# Read genotype table: POS numeric; GT columns character
geno <- read_tsv(paste0(OUT, ".genotypes.tsv"),
                 col_names = colnames_geno,
                 col_types = cols(
                   CHROM = col_character(),
                   POS   = col_double(),
                   .default = col_character()
                 ),
                 progress = FALSE)

# Long format
long <- geno |>
  pivot_longer(cols = all_of(samples),
               names_to = "individual",
               values_to = "GT")

# Heterozygote flag (accepts / or |; treats ./. as NA)
is_het <- function(gt) {
  ifelse(is.na(gt) | gt %in% c("./.", ".|."), NA,
         {
           m <- str_match(gt, "^([0-9.]+)[/|]([0-9.]+)$")
           ifelse(is.na(m[,1]), NA,
                  (m[,2] != m[,3]) & (m[,2] != ".") & (m[,3] != "."))
         })
}

# Bin positions to 50 kb windows
hets <- long |>
  mutate(
    het = is_het(GT),
    win_start = floor((POS - 1) / WIN) * WIN,
    win_mid   = win_start + WIN/2
  ) |>
  group_by(individual, CHROM, win_start, win_mid) |>
  summarise(
    n_sites = sum(!is.na(het)),
    n_het   = sum(het, na.rm = TRUE),
    het_rate = n_het / n_sites,
    .groups = "drop"
  )

# --- Species assignment ---
# If IDs start with species code (e.g., "BFBO501"), use substring:
hets <- hets |>
  mutate(species = substr(individual, 1, 4))
# If you have a metadata table instead, do something like:
# meta <- read_tsv("samples_metadata.tsv")  # columns: individual, species
# hets <- hets |> left_join(meta, by="individual")

# Optional: filter out low-data windows (helps noisy tails)
min_sites <- 100                         # tweak to your data density
hets_f <- hets |>
  filter(n_sites >= min_sites)

# --- Build a per-species palette mapping: individual -> color ---
# Use a strong qualitative palette; repeat per species (that's fine)
pal_for <- function(n) grDevices::hcl.colors(n, palette = "Dynamic")

inds_by_sp <- hets_f |>
  distinct(species, individual) |>
  arrange(species, individual) |>
  group_by(species) |>
  summarise(ind_list = list(individual), .groups = "drop")

# one palette *per* species; same colors can recur in a different species facet
pal_tbl <- inds_by_sp |>
  rowwise() |>
  mutate(cols = list(setNames(pal_for(length(ind_list)), ind_list))) |>
  ungroup()

pal_individual <- unlist(pal_tbl$cols)  # named vector: names = individual, values = hex colors

sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"

sex <- read_csv(sex_file, col_names = c("individual", "sex"))

# --- Join with sex and chr ---
hets_sex <- hets_f %>%
  left_join(sex, by = "individual") %>%                    # add sex
  mutate(
    sex = str_squish(sex),                                 # trim any stray spaces
    sex = coalesce(sex, "Unknown")                         # replace NA with "Unknown"
  )


# --- Plot ---
p <- ggplot(hets_sex,
            aes(x = win_mid, y = het_rate, color = sex, group = individual)) +
  geom_line(linewidth = 0.2, alpha = 0.95) +
  scale_color_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    title = paste0("Chromosome 29 (50 kb bins)"),
    x = "Position (bp)", y = "Heterozygosity", color = "Individual"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "right",
    legend.key.height = unit(0.35, "lines"),
    legend.text = element_text(size = 7)
  )


ggsave(paste0(OUT, ".hets_50kb.pdf"), p, width = 12, height = 3, dpi = 300)

# plot colored as individual

# --- Plot ---
p <- ggplot(hets_sex,
            aes(x = win_mid, y = het_rate, color = individual, group = individual)) +
  geom_line(linewidth = 0.2, alpha = 0.95) +
  labs(
    title = paste0("Chromosome 29 (50 kb bins)"),
    x = "Position (bp)", y = "Heterozygosity", color = "Individual"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "right",
    legend.key.height = unit(0.35, "lines"),
    legend.text = element_text(size = 7)
  )


ggsave(paste0(OUT, ".hets_50kb.individual.pdf"), p, width = 12, height = 3, dpi = 300)

library(patchwork)

plots_facets <- hets_f %>%
  split(.$species) %>%
  lapply(function(df) {
    inds <- sort(unique(df$individual))
    pal  <- setNames(pal_for(length(inds)), inds)

    ggplot(df, aes(x = win_mid, y = het_rate,
                   color = individual, group = individual)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      scale_color_manual(values = pal) +
      labs(title = unique(df$species), x = "Position (bp)", y = NULL, color = NULL) +
      theme_bw(base_size = 10) +
      theme(legend.position = "right",
            legend.key.height = unit(0.4, "lines"),
            plot.margin = margin(2, 5, 2, 5))
  })

wrap_plots(plots_facets, ncol = 1)
ggsave(paste0(OUT, ".hets_50kb.faceted.pdf"), width = 12, height = 8, dpi = 300)
