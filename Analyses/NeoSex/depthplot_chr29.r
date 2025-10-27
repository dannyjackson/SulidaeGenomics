#!/usr/bin/env bash
set -euo pipefail

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth_filtered/depthplot_chr29

CHR="CM062595.1"   # chromosome to analyze
WIN=50000          # window size in bp (10k–100k typical)
INPAT="*_${CHR}_depth.txt"
OUTDIR="depth_summaries"
mkdir -p "$OUTDIR"


# If your files are gzipped, set READER to "zcat"
READER="cat"

# Process each matching file one-by-one
find . -maxdepth 1 -type f -name "$INPAT" -print0 |
while IFS= read -r -d '' F; do
  IND="$(basename "$F" | cut -d'_' -f1)"
  OUT="${OUTDIR}/${IND}.${CHR}.bins${WIN}.tsv"
  echo "Binning $F  ->  $OUT"

  $READER "$F" | awk -v W="$WIN" -v IND="$IND" -v CHR="$CHR" -v OFS="\t" '
    BEGIN { FS = OFS = "\t"; sum = cnt = 0; cur = -1 }
    $1 == CHR {
      bin = int($2 / W)
      if (cur < 0) { cur = bin }                 # initialize
      if (bin != cur) {
        if (cnt > 0) print IND, CHR, cur*W, (cur+1)*W - 1, sum / cnt
        cur = bin; sum = cnt = 0
      }
      sum += $3; cnt++
    }
    END {
      if (cnt > 0) print IND, CHR, cur*W, (cur+1)*W - 1, sum / cnt
    }' > "$OUT"
done

# Combine all individuals
COMB="${OUTDIR}/binned_depth.${CHR}.WIN${WIN}.tsv"
printf "individual\tchrom\tstart\tend\tmean_depth\n" > "$COMB"
cat "${OUTDIR}"/*."${CHR}".bins${WIN}.tsv >> "$COMB"
echo "Wrote: $COMB"

######################################
# Keep ~0.2% of rows (tune P between 0.001 and 0.01 depending on file size)
#!/usr/bin/env bash


# Loop over files serially
find . -maxdepth 1 -type f -name "$INPAT" -print0 |
while IFS= read -r -d '' F; do
  IND="$(basename "$F" | cut -d'_' -f1)"
  OUT="${OUTDIR}/${IND}.${CHR}.sample.tsv"
  echo "Sampling $F -> $OUT (p=${P})"

  $READER "$F" | awk -v p="$P" -v IND="$IND" -v CHR="$CHR" -v OFS="\t" '
    BEGIN { FS = OFS = "\t"; srand() }          # for reproducibility, use: srand(12345)
    ($1 == CHR) && (rand() < p) { print IND, $1, $2, $3 }
  ' > "$OUT"
done

# Combine into one tidy table
{
  printf "individual\tchrom\tpos\tdepth\n"
  cat "$OUTDIR"/*."${CHR}".sample.tsv
} > "$OUTDIR/sampled_depth.${CHR}.p${P}.tsv"

echo "Wrote: $OUTDIR/sampled_depth.${CHR}.p${P}.tsv"

######################################
# PLOT
######################################

library(data.table)
library(ggplot2)

# Inputs from steps 1–2
bins_file <- "depth_summaries/binned_depth.CM062595.1.WIN50000.tsv"
samp_file <- "depth_samples/sampled_depth.CM062595.1.p0.002.tsv"  # optional

bins <- fread(bins_file)
bins[, bin_mid := (start + end) / 2]

# OPTIONAL: normalize per individual (e.g., divide by that individual’s median bin depth)
bins[, med_ind := median(mean_depth), by = individual]
bins[, depth_norm := mean_depth / med_ind]

# If you have a 2-column mapping file "individual<TAB>sex", join it here:
# sexmap <- fread("sex_map.tsv", col.names = c("individual","sex"))
# bins <- sexmap[bins, on="individual"]

# Base plot: per-individual smoothed mean depth across the chromosome

p1 <- ggplot(bins, aes(x = bin_mid, y = depth_norm, group = individual)) +
  geom_line(alpha = 0.8, linewidth = 0.3) +
  facet_wrap(~ individual, ncol = 4, scales = "free_y") +
  labs(x = "Position on CM062595.1 (bp)",
       y = "Normalized mean depth (per individual median = 1)",
       title = "Depth profiles by individual (binned, fast)") +
  theme_bw(base_size = 10)


print(p1)
ggsave("depth_profiles.CM062595.1.binned_norm.png", p1, width = 12, height = 8, dpi = 300)


bins[, species := sub("([A-Za-z]+)[0-9]+", "\\1", individual)]

# Example assuming `bins` has columns: individual, chrom, start, end, bin_mid, depth_norm, species

p1 <- ggplot(bins, aes(x = bin_mid, y = depth_norm, color = species, group = individual)) +
  geom_line(alpha = 0.6, linewidth = 0.4) +
  labs(
    x = "Position on CM062595.1 (bp)",
    y = "Normalized mean depth (per individual median = 1)",
    title = "Depth profiles across CM062595.1 by species",
    color = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90")
  )

print(p1)
ggsave("depth_profiles.CM062595.1.all_species.png", p1, width = 12, height = 6, dpi = 300)



########
# facet by species
########

# install.packages("ggnewscale")
library(ggnewscale)


library(data.table)
library(ggplot2)
library(patchwork)

# If needed
# bins[, species := sub("([A-Za-z]+)[0-9]+", "\\1", individual)]

# Distinct palette per species (good separation even with many IDs)
pal_for <- function(n) grDevices::hcl.colors(n, palette = "Dynamic")

make_species_plot <- function(sp, data = bins) {
  d <- data[species == sp]
  # stable order for colors/legend
  inds <- sort(unique(d$individual))
  pal  <- setNames(pal_for(length(inds)), inds)

  ggplot(d, aes(x = bin_mid, y = depth_norm,
                color = individual, group = individual)) +
    geom_line(alpha = 0.85, linewidth = 0.5) +
    scale_color_manual(values = pal, guide = guide_legend(ncol = 1, override.aes = list(linewidth = 1))) +
    labs(title = sp, x = "Position (bp)", y = "Normalized depth", color = "Individual") +
    theme_bw(base_size = 10) +
    theme(
      # put a small legend inside the plotting area (top-right)
      legend.position      = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = scales::alpha("white", 0.75), colour = "grey80"),
      legend.key.height    = unit(9, "pt"),
      legend.key.width     = unit(16, "pt"),
      legend.text          = element_text(size = 7),
      legend.title         = element_text(size = 8),
      plot.title           = element_text(hjust = 0.5, face = "bold")
    )
}

plots <- lapply(sort(unique(bins$species)), make_species_plot)

p_all <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Depth profiles by species (distinct colors + per-panel legends)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("depth_profiles.CM062595.1.faceted_species_indiv_legend_inpanel.png",
       p_all, width = 12, height = 24, dpi = 300)


###########################
# plot without outliers
###########################

library(scales)   # for squish

# Distinct palette per species (good separation even with many IDs)
pal_for <- function(n) grDevices::hcl.colors(n, palette = "Dynamic")

make_species_plot <- function(sp, data = bins, cover = 0.98) {
  d <- data[species == sp]
  inds <- sort(unique(d$individual))
  pal  <- setNames(pal_for(length(inds)), inds)

  # Robust per-species y-limits: keep central `cover` fraction (e.g., 98%)
  alpha <- 1 - cover
  lims  <- quantile(d$depth_norm, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  ggplot(d, aes(x = bin_mid, y = depth_norm, color = individual, group = individual)) +
    geom_line(alpha = 0.85, linewidth = 0.25) +
    scale_color_manual(values = pal,
                       guide = guide_legend(ncol = 1, override.aes = list(linewidth = 1))) +
    scale_y_continuous(limits = lims, oob = squish) +  # << key line
    labs(title = sp, x = "Position (bp)", y = "Normalized depth", color = "Individual") +
    theme_bw(base_size = 10) +
    theme(
      legend.position      = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = scales::alpha("white", 0.75), colour = "grey80"),
      legend.key.height    = unit(9, "pt"),
      legend.key.width     = unit(16, "pt"),
      legend.text          = element_text(size = 7),
      legend.title         = element_text(size = 8),
      plot.title           = element_text(hjust = 0.5, face = "bold")
    )
}

plots <- lapply(sort(unique(bins$species)), make_species_plot)  # default cover = 0.98

p_all <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Depth profiles by species (central 98% shown; outliers squished)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("depth_profiles.CM062595.1.faceted_species_indiv_legend_inpanel.trim98.png",
       p_all, width = 12, height = 24, dpi = 300)

