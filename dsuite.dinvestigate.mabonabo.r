# plot MABO NABO

# Load packages

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

data$chromo <- factor(data$chromo, levels = c(1:28, 30:32, "Z"))
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
fdm_cut <- quantile(plot_data[[metric]], probs = 0.999, na.rm = TRUE)

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



GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
############################################
# Are there genes within these signals? 
############################################

WIN=MABO_AtlCar_MABO_IndoPacific_NABO.f_dM.50_25.top0.01pct.outliers.tsv

awk -v FS='\t' -v OFS='\t' 'NR>1 {
  chr=$2; ws=$5; we=$6;
  print chr, ws-1, we, "win_"NR-1 "|" $1 "|" $3   # name keeps track of source row
}' "$WIN" > windows.bed

awk -v FS='\t' -v OFS='\t' '$3=="gene" {
  id=""; name="";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^ID=/)   id=substr(a[i],4);
    if(a[i] ~ /^Name=/) name=substr(a[i],6);
  }
  if(name=="") name=id;
  print $1, $4-1, $5, name
}' "$GFF" > genes.bed

bedtools intersect -a genes.bed -b windows.bed -wa -wb > gene_window_overlaps.tsv

awk '{print $1,$2,$4}' gene_window_overlaps.tsv | sort | uniq

CM062567.1 39078937 WIPF3
CM062567.1 39106091 SCRN1
CM062567.1 64684597 CTDP1
CM062568.1 123724941 BBX
CM062568.1 31314558 FNDC3A
CM062568.1 58753424 HS6ST3
CM062569.1 49528678 AIG1
CM062569.1 54326766 TAAR1
CM062569.1 86566520 EYS
CM062570.1 1878992 NPFFR2
CM062571.1 58902249 TSPAN12
CM062572.1 63963401 CACNA1H
CM062573.1 12446351 bMorBas2_egapxtmp_005019
CM062575.1 28217304 PTPRT
CM062581.1 22878774 HCN4
CM062585.1 397695 THUMPD2

# Two genes are putatively overrepresented and associated with GO term monoatomic ion transport; HCN4 and CACNA1H
# This could relate to salt tolerance, as the salt levels differ between AtlCar (higher) and IndoPacific (lower)

############################################

############################################
# Are there genes nearby these signals? 10kb either
############################################


WIN=MABO_AtlCar_MABO_IndoPacific_NABO.f_dM.50_25.top0.01pct.outliers.tsv

awk -v FS='\t' -v OFS='\t' 'NR>1 {
  chr=$2; ws=$5-10000; we=$6+10000;
  print chr, ws-1, we, "win_"NR-1 "|" $1 "|" $3   # name keeps track of source row
}' "$WIN" > windows.10kb.bed

awk -v FS='\t' -v OFS='\t' '$3=="gene" {
  id=""; name="";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^ID=/)   id=substr(a[i],4);
    if(a[i] ~ /^Name=/) name=substr(a[i],6);
  }
  if(name=="") name=id;
  print $1, $4-1, $5, name
}' "$GFF" > genes.10kb.bed

bedtools intersect -a genes.10kb.bed -b windows.10kb.bed -wa -wb > gene_window_overlaps.10kb.tsv

awk '{print $1, $2, $4}' gene_window_overlaps.10kb.tsv | sort | uniq

# Chromosome 1
CM062567.1 39078937 WIPF3
CM062567.1 39106091 SCRN1
CM062567.1 64684597 CTDP1

# Chromosome 2
CM062568.1 31314558 FNDC3A
CM062568.1 58753424 HS6ST3
CM062568.1 123724941 BBX

# Chromosome 3
CM062569.1 49528678 AIG1
CM062569.1 54326766 TAAR1
CM062569.1 54339056 bMorBas2_egapxtmp_000806 # trace amine-associated receptor 2-like; TAAR2
CM062569.1 86566520 EYS

# Chromosome 4
CM062570.1 1878992 NPFFR2

# Chromosome 5
CM062571.1 29862725 IFT56
CM062571.1 58902249 TSPAN12
# Chromosome 6
CM062572.1 63963401 CACNA1H
# Chromosome 7
CM062573.1 12446351 bMorBas2_egapxtmp_005019 # ankyrin repeat and KH domain-containing protein 1 isoform 12-T12; ANKHD1
# Chromosome 9
CM062575.1 28217304 PTPRT
# Chromosome 12
CM062578.1 25668204 GPR148 # 0.8297822
# Chromosome 15
CM062581.1 22878774 HCN4
# Chromosome 19
CM062585.1 397695 THUMPD2
CM062585.1 411988 GPCPD1

# GPR148 -- detection of smell!