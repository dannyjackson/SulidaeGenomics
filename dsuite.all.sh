# Dsuite all
# filter vcf to just autosomes

#!/usr/bin/env bash
#SBATCH --job-name=filtervcf
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/filtervcf%A_%a.out
#SBATCH --mail-type=ALL

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz
VCF_OUT=/xdisk/mcnew/dannyjackson/sulidae/datafiles/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.autosomes.vcf.gz

module load bcftools


bcftools view -t ^CM062595.1,CM062600.1,CM062610.1,CM062599.1 -Oz -o "${VCF_OUT}" "${VCF}"

# Dsuite
# Save this as SETS.txt\
BFBO501 BFBO_southern
BFBO502 BFBO_southern
BFBO503 BFBO_southern
BFBO504 BFBO_GofCA
BFBO505 BFBO_southern
BRBO201 BRBO_Pacific
BRBO202 BRBO_Pacific
BRBO203 BRBO_AtlCar
BRBO205 BRBO_AtlCar
MABO302 MABO_IndoPacific
MABO304 MABO_AtlCar
MABO305 MABO_AtlCar
MABO306 MABO_IndoPacific
NABO402 NABO
NABO403 NABO
NABO404 NABO
NABO405 NABO
NABO406 NABO
PEBO601 PEBO
PEBO603 PEBO
PEBO604 PEBO
PEBO605 PEBO
PEBO606 PEBO
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup


#!/usr/bin/env bash
#SBATCH --job-name=all_dsuite
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/all_dsuite.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch all_dsuite.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite
# all
mkdir -p all
cd all

VCF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t all.nwk

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all
~/programs/Dsuite/Build/Dsuite Fbranch all.nwk SETS_tree.txt > all_Fbranch.txt

python3 ~/programs/Dsuite/utils/dtools.py all_Fbranch.txt all.nwk --tree-label-size 5 --ladderize --use_distances


echo -e 'MABO_AtlCar\tMABO_IndoPacific\tNABO' > test_trios.txt
echo -e 'BFBO_GofCA\tBFBO_southern\tPEBO' >> test_trios.txt

# Dinvestigate MABO NABO and BFBO PEBO
~/programs/Dsuite/Build/Dsuite Dinvestigate -w 50,25 $VCF SETS.txt test_trios.txt

# plot MABO NABO

# Load packages
.libPaths("~/R/library_elgato")

library(dplyr)
library(readr)
library(ggplot2)
library(rlang)
library(stringr)

# read in the results with 50 SNP windows and a step of 25 SNPs
bigStep <- read.table("MABO_AtlCar_MABO_IndoPacific_NABO_localFstats__50_25.txt",as.is=T,header=T)

# install.packages("stringr")




# -------------------------------
# Inputs (edit if needed)
# -------------------------------
infile   <- "MABO_AtlCar_MABO_IndoPacific_NABO_localFstats__50_25.txt"
convfile <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"

# Derive pop_name + window label from infile (robust to extra dots/underscores)
# Expected pattern: <pop_name>_localFstats__<win>.txt
pop_name <- sub("_localFstats__.*$", "", basename(infile))
win      <- sub("^.*_localFstats__", "", basename(infile)) |> sub("\\.txt$", "", x = _)

# Output base directory; you can change this to your project root
outdir <- getwd()

# Which metrics to plot
metrics <- c("D", "f_d", "f_dM")

# Cutoffs (edit to taste)
metric_cutoffs <- c(D = 0, `f_d` = 0, `f_dM` = 0)

# Alternating colors
color1 <- "#444444"
color2 <- "#1f78b4"

# Desired chromosome order
chromo_levels <- c(1, "1A", 2:4, "4A", 5:29, "Z") |> as.character()

# -------------------------------
# Read data
# -------------------------------
cat("Reading bigStep...\n")
data <- read.table(infile, header = TRUE, as.is = TRUE)

# Ensure required columns exist
req_cols <- c("chr", "windowStart", "windowEnd")
stopifnot(all(req_cols %in% names(data)))

# Use midpoints as 'position'
data <- data %>%
  mutate(position = (windowStart + windowEnd) / 2)

# -------------------------------
# Read chromosome conversion map
# Expect two columns: original scaffold ID and desired chromo label.
# It can be headered or headerless; handle both.
# -------------------------------
cat("Reading chromosome conversion map...\n")
conv_guess <- read_csv(convfile, col_names = c("chromo", "chr"), show_col_types = FALSE)

# Coerce to character
conv_guess <- conv_guess %>%
  transmute(chr = as.character(chr),
            chromo = as.character(chromo))

# -------------------------------
# Join mapping
# -------------------------------
cat("Applying chromosome renaming...\n")
data <- data %>%
  mutate(chr = as.character(chr)) %>%
  left_join(conv_guess, by = "chr")

# Sanity check
if (any(is.na(data$chromo))) {
  missing_ids <- unique(data$chr[is.na(data$chromo)])
  warning("Some scaffolds did not map to a chromosome label. Showing first few:\n",
          paste(head(missing_ids, 10), collapse = ", "))
  # If desired, drop unmapped
  data <- data %>% filter(!is.na(chromo))
}

data <- data %>%
  mutate(chromo = as.character(chromo)) %>%
  # remove unwanted scaffolds explicitly
  filter(!chr %in% c("CM062595.1","CM062599.1","CM062600.1","CM062610.1")) 

data$chromo <- factor(data$chromo, levels = c(1, "1A", 2:4, "4A", 5:33, "Z"))
# -------------------------------
# Prepare cumulative BP coordinates (your template)
# -------------------------------
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
  summarize(center = mean(BPcum), .groups = "drop")

# -------------------------------
# Plot function (matches your style)
# -------------------------------
plot_metric <- function(metric) {
  stopifnot(metric %in% names(plot_data))
  metric_cutoff <- metric_cutoffs[[metric]]
  metric_dir <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", win))
  dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Generating plot for", metric, "...\n")
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
  cat("Saved:", outfile, "\n")
}

# -------------------------------
# Run for each requested metric
# -------------------------------
for (m in metrics) plot_metric(m)

cat("Script completed successfully!\n")



# -------------------------------
# Filter by top 0.001 f_dM
# -------------------------------
library(readr)

# --- choose metric + cutoff ---
metric <- "f_dM"  # <- change if needed
# global 99.99th percentile across all chromosomes
fdm_cut <- quantile(plot_data[[metric]], probs = 0.9999, na.rm = TRUE)

# flag outliers
plot_data <- plot_data %>%
  mutate(is_outlier = .data[[metric]] >= fdm_cut)

# write outliers to file
out_dir  <- file.path(outdir, "analyses", metric, paste0(pop_name, "/", win))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, paste0(pop_name, ".", metric, ".", win, ".top0.01pct.outliers.tsv"))
plot_data %>%
  filter(is_outlier) %>%
  # keep useful columns; add cutoff used
  mutate(cutoff_used = fdm_cut) %>%
  select(chromo, chr, position, BPcum, windowStart, windowEnd, D, f_d, f_dM, d_f, cutoff_used) %>%
  arrange(chromo, position) %>%
  write_tsv(out_file)

message("Outliers written to: ", out_file)

# --- plotting: base layer for normal points (alternating chrom colors) ---
base_palette <- rep(c(color1, color2), length.out = length(levels(plot_data$chromo)))

p <- ggplot() +
  # non-outliers in alternating chrom colors
  geom_point(
    data = subset(plot_data, !is_outlier),
    aes(x = BPcum, y = .data[[metric]], color = chromo),
    alpha = 0.8, size = 1, show.legend = FALSE
  ) +
  scale_color_manual(values = setNames(base_palette, levels(plot_data$chromo))) +

  # outliers in red on top
  geom_point(
    data = subset(plot_data, is_outlier),
    aes(x = BPcum, y = .data[[metric]]),
    color = "red3", alpha = 0.95, size = 1.6
  ) +

  # x-axis ticks/labels from axisdf
  scale_x_continuous(
    labels = axisdf$chromo,
    breaks = axisdf$center,
    guide  = guide_axis(n.dodge = 2)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = fdm_cut, linetype = "dashed", linewidth = 0.5, color = "red3") +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

outfile_png <- file.path(out_dir, paste0(pop_name, ".", metric, ".", win, ".top0.01pct.sigline.png"))
ggsave(outfile_png, p, width = 20, height = 5, units = "in")

message("Saved plot: ", outfile_png)
message("f_dM 99.99th percentile cutoff: ", signif(fdm_cut, 4))




chromo	chr	position	BPcum	windowStart	windowEnd	D	f_d	f_dM	d_f	cutoff_used
1	CM062567.1	64794200	64794200	64793825	64794575	0.103548	0.911685	0.911685	0.053628	0.7324406396000691
1	CM062567.1	73516758.5	73516758.5	73516464	73517053	0.107023	0.992556	0.948429	0.058201	0.7324406396000691
2	CM062568.1	110824328	280963631	110824025	110824631	0.283058	0.827295	0.738943	0.164862	0.7324406396000691
2	CM062568.1	110824612.5	280963915.5	110824373	110824852	0.44246	0.86838	0.8324	0.288486	0.7324406396000691
2	CM062568.1	110824931	280964234	110824668	110825194	0.240404	0.810627	0.753642	0.137891	0.7324406396000691



GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
############################################
# Are there genes within these signals? 
############################################
# first signal of introgression
awk '$1 == "CM062567.1" && $4 >= 64793825 && $5 <= 64794575' \
  $GFF | grep 'ID=gene'

# second 
awk '$1 == "CM062567.1" && $4 >= 73516464 && $5 <= 73517053' \
  $GFF | grep 'ID=gene'

# third
awk '$1 == "CM062567.1" && $4 >= 110824025 && $5 <= 110824631' \
  $GFF | grep 'ID=gene'

# fourth
awk '$1 == "CM062567.1" && $4 >= 110824373 && $5 <= 110824852' \
  $GFF | grep 'ID=gene'

# fifth
awk '$1 == "CM062567.1" && $4 >= 110824668 && $5 <= 110825194' \
  $GFF | grep 'ID=gene'

# No.
############################################

############################################
# Are there genes nearby these signals? 10kb either
############################################

# first signal of introgression
echo 'one'
awk '$1 == "CM062567.1" && $4 >= 64593825 && $5 <= 64994575' \
  $GFF | grep 'ID=gene'
# CM062567.1      Liftoff gene    64684546        64808864        .       -       .       ID=gene-CTDP1;Dbxref=GeneID:104050901;Name=CTDP1;description=CTD phosphatase subunit 1;gbkey=Gene;gene=CTDP1;gene_biotype=protein_coding;coverage=0.993;sequence_ID=0.951;valid_ORFs=3;extra_copy_number=0;copy_num_ID=gene-CTDP1_0

echo 'two'
# second 
awk '$1 == "CM062567.1" && $4 >= 73416464 && $5 <= 73617053' \
  $GFF | grep 'ID=gene'
# CM062567.1      Liftoff gene    73531157        73586898        .       -       .       ID=gene-DSP;Dbxref=GeneID:104049457;Name=DSP;description=desmoplakin;gbkey=Gene;gene=DSP;gene_biotype=protein_coding;coverage=0.998;sequence_ID=0.973;valid_ORFs=4;extra_copy_number=0;copy_num_ID=gene-DSP_0

echo 'three'
# third
awk '$1 == "CM062567.1" && $4 >= 110724025 && $5 <= 110924631' \
  $GFF | grep 'ID=gene'
# CM062567.1      Liftoff gene    110820478       110901791       .       -       .       ID=gene-CEP192;Dbxref=GeneID:104047969;Name=CEP192;description=centrosomal protein 192;gbkey=Gene;gene=CEP192;gene_biotype=protein_coding;coverage=0.998;sequence_ID=0.941;valid_ORFs=0;extra_copy_number=0;copy_num_ID=gene-CEP192_0
# CM062567.1      Liftoff gene    110904891       110918376       .       -       .       ID=gene-SEH1L;Dbxref=GeneID:104040201;Name=SEH1L;description=SEH1 like nucleoporin;gbkey=Gene;gene=SEH1L;gene_biotype=protein_coding;coverage=0.993;sequence_ID=0.965;valid_ORFs=2;extra_copy_number=0;copy_num_ID=gene-SEH1L_0

echo 'four'
# fourth
awk '$1 == "CM062567.1" && $4 >= 110724373 && $5 <= 110924852' \
  $GFF | grep 'ID=gene'
#CM062567.1      Liftoff gene    110820478       110901791       .       -       .       ID=gene-CEP192;Dbxref=GeneID:104047969;Name=CEP192;description=centrosomal protein 192;gbkey=Gene;gene=CEP192;gene_biotype=protein_coding;coverage=0.998;sequence_ID=0.941;valid_ORFs=0;extra_copy_number=0;copy_num_ID=gene-CEP192_0
#CM062567.1      Liftoff gene    110904891       110918376       .       -       .       ID=gene-SEH1L;Dbxref=GeneID:104040201;Name=SEH1L;description=SEH1 like nucleoporin;gbkey=Gene;gene=SEH1L;gene_biotype=protein_coding;coverage=0.993;sequence_ID=0.965;valid_ORFs=2;extra_copy_number=0;copy_num_ID=gene-SEH1L_0

echo 'five'
# fifth
awk '$1 == "CM062567.1" && $4 >= 110724668 && $5 <= 110925194' \
  $GFF | grep 'ID=gene'
#CM062567.1      Liftoff gene    110820478       110901791       .       -       .       ID=gene-CEP192;Dbxref=GeneID:104047969;Name=CEP192;description=centrosomal protein 192;gbkey=Gene;gene=CEP192;gene_biotype=protein_coding;coverage=0.998;sequence_ID=0.941;valid_ORFs=0;extra_copy_number=0;copy_num_ID=gene-CEP192_0
#CM062567.1      Liftoff gene    110904891       110918376       .       -       .       ID=gene-SEH1L;Dbxref=GeneID:104040201;Name=SEH1L;description=SEH1 like nucleoporin;gbkey=Gene;gene=SEH1L;gene_biotype=protein_coding;coverage=0.993;sequence_ID=0.965;valid_ORFs=2;extra_copy_number=0;copy_num_ID=gene-SEH1L_0


CTDP1
DSP
CEP192
SEH1L