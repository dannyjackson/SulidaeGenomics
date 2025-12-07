

# do it with brbo vs 4 taxa
BFBO501 4taxa
BFBO502 4taxa
BFBO503 4taxa
BFBO504 4taxa
BFBO505 4taxa
BRBO201 BRBO_Pacific
BRBO202 BRBO_Pacific
BRBO203 BRBO_AtlCar
BRBO205 BRBO_AtlCar
MABO302 4taxa
MABO304 4taxa
MABO305 4taxa
MABO306 4taxa
NABO402 4taxa
NABO403 4taxa
NABO404 4taxa
NABO405 4taxa
NABO406 4taxa
PEBO601 4taxa
PEBO603 4taxa
PEBO604 4taxa
PEBO605 4taxa
PEBO606 4taxa
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup

echo "(((BRBO_AtlCar,BRBO_Pacific),4taxa),Outgroup);" > BRBO_4taxa.nwk

#!/usr/bin/env bash
#SBATCH --job-name=BRBO4taxa_dsuite
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBO4taxa_dsuite.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO4taxa_dsuite.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite
# all
mkdir -p BRBO_4taxa/
cd BRBO_4taxa/

VCF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t BRBO_4taxa.nwk

~/programs/Dsuite/Build/Dsuite Fbranch BRBO_4taxa.nwk SETS_tree.txt > BRBO_4taxa_Fbranch.txt


# Redo with likelihoods
#!/usr/bin/env bash
#SBATCH --job-name=BRBO4taxa_dsuite
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBO4taxa_dsuite.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO4taxa_dsuite.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO_4taxa_likelihoods

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t BRBO_4taxa.nwk --ABBAclustering -g
~/programs/Dsuite/Build/Dsuite Fbranch BRBO_4taxa.nwk SETS_tree.txt > BRBO_4taxa_Fbranch.txt


module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py BRBO_4taxa_Fbranch.txt BRBO_4taxa.nwk --tree-label-size 3 --ladderize --use_distances


#!/usr/bin/env bash
#SBATCH --job-name=BRBO4taxa_dinvestigate
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --output=slurm_output/BRBO4taxa_dinvestigate.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO4taxa_dinvestigate.sh

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

echo -e 'BRBO_Pacific\tBRBO_AtlCar\t4taxa' > test_trios.txt
~/programs/Dsuite/Build/Dsuite Dinvestigate -g -w 50,25 $VCF SETS.txt test_trios.txt

gunzip /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf

ruby get_fixed_site_gts.rb $VCF BRBO_4taxa.fixed.txt BRBO201,BRBO202 BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,MABO302,MABO304,MABO305,MABO306,NABO402,NABO403,NABO404,NABO405,NABO406,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_4taxa.fixed.txt BRBO_4taxa.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_BFBO.fixed.txt BRBO201,BRBO202 BFBO501,BFBO502,BFBO503,BFBO504,BFBO505 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_BFBO.fixed.txt BRBO_BFBO.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_MABO.fixed.txt BRBO201,BRBO202 MABO302,MABO304,MABO305,MABO306 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_MABO.fixed.txt BRBO_MABO.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_NABO.fixed.txt BRBO201,BRBO202 NABO402,NABO403,NABO404,NABO405,NABO406 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_NABO.fixed.txt BRBO_NABO.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_PEBO.fixed.txt BRBO201,BRBO202 PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_PEBO.fixed.txt BRBO_PEBO.fixed.svg 1.0 1000





# Load packages
.libPaths("~/R/library_elgato")

library(ggplot2)
library(dplyr)


# Read data
bigStep <- read.table(
  "BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__50_25.txt",
  header = TRUE, as.is = TRUE
)


# install.packages("stringr")


library(dplyr)
library(readr)
library(ggplot2)
library(rlang)
library(stringr)


# -------------------------------
# Inputs (edit if needed)
# -------------------------------
infile   <- "BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__50_25.txt"
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
4	CM062570.1	76150455.5	499231318.5	76150228	76150683	0.14719	0.575061	0.55087	0.067784	0.5003591676001057
13	CM062579.1	26413271	921799107	26412558	26413984	0.249313	0.625807	0.607195	0.105941	0.5003591676001057
22	CM062588.1	930854.5	1065245245	920676	941033	0.524932	0.620804	0.519498	0.253805	0.5003591676001057
22	CM062588.1	940591.5	1065254982	939221	941962	0.475157	0.626907	0.522118	0.229384	0.5003591676001057
23	CM062589.1	292582	1073795662.5	292358	292806	0.81585	0.533434	0.516381	0.587031	0.5003591676001057


GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
############################################
# Are there genes within these signals? 
############################################
# first signal of introgression
awk '$1 == "CM062570.1" && $4 >= 76150228 && $5 <= 76150683' \
  $GFF | grep 'ID=gene'

# second signal of introgression
awk '$1 == "CM062579.1" && $4 >= 26412558 && $5 <= 26413984' \
  $GFF | grep 'ID=gene'

# third signal of introgression
awk '$1 == "CM062588.1" && $4 >= 920676 && $5 <= 941033' \
  $GFF | grep 'ID=gene'
# CM062588.1      Liftoff gene    932580  932874  .       +       .       ID=gene-LOC135317588;Dbxref=GeneID:135317588;Name=LOC135317588;description=uncharacterized LOC135317588;gbkey=Gene;gene=LOC135317588;gene_biotype=lncRNA;coverage=0.540;sequence_ID=0.470;extra_copy_number=0;copy_num_ID=gene-LOC135317588_0;low_identity=True

# fourth signal of introgression
awk '$1 == "CM062588.1" && $4 >= 939221 && $5 <= 941962' \
  $GFF | grep 'ID=gene'

# fifth signal of introgression
awk '$1 == "CM062589.1" && $4 >= 292358 && $5 <= 292806' \
  $GFF | grep 'ID=gene'