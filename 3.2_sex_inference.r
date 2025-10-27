# Sex Inference

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats

mkdir sexing_by_depth
cd sexing_by_depth

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/

BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/sortedbamfiles

#!/bin/bash
#SBATCH --account=mcnew
#SBATCH --job-name=depthstats
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/depthstats.%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
# NOTE: set the array at submission time (see below)

set -euo pipefail

# ---- Required env vars (export before sbatch or edit here) ----
PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
OUTDIR=${PROJDIR}/datafiles/bamstats/sexing_by_depth
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
LIST_FILE=${PROJDIR}/raw_sequences/filenames_SRA.txt

# Pick the Nth (array index) line from the list
bird="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_FILE" | tr -d '\r')"

# Skip empty/comment lines defensively
if [[ -z "${bird}" || "${bird}" =~ ^[[:space:]]*# ]]; then
  echo "Empty or comment line at index ${SLURM_ARRAY_TASK_ID}; nothing to do."
  exit 0
fi

module load samtools || true

bam="${BAMDIR}/${bird}.sorted.bam"

if [[ ! -s "$bam" ]]; then
  echo "ERROR: BAM not found or empty: $bam" >&2
  exit 2
fi

# Write sample name as first line, then append depths
out="${OUTDIR}/${bird}_depthstats.txt"
echo "${bird}" > "$out"
samtools depth "${BAMDIR}/${bird}.sorted.bam" >> "$out"



PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/

N=$(grep -v -E '^[[:space:]]*(#|$)' "${PROJDIR}/raw_sequences/filenames_SRA.txt" | wc -l)
sbatch --array=1-"$N" depth.sh



PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
OUTDIR=${PROJDIR}/datafiles/bamstats/sexing_by_depth
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
LIST_FILE=${PROJDIR}/raw_sequences/filenames_SRA.txt

ls ${BAMDIR}/*.sorted.bam | wc -l
ls ${BAMDIR}/*.sorted.bam.bai | wc -l
while read -r bird; do
samtools index ${BAMDIR}/"$bird".sorted.bam 
done <  ${PROJDIR}/raw_sequences/filenames_SRA.txt 

# get mean depth
# Inputs
LIST="${PROJDIR}/raw_sequences/filenames_SRA.txt"
# Build the file paths once, then process them all in a single awk
sed "s|^|${OUTDIR}/|; s|$|_depthstats.txt|" "$LIST" \
| xargs -r awk 'BEGINFILE{sum=0;n=0} {sum+=$3;n++} ENDFILE{printf "%s\t%.6f\n", FILENAME, (n? sum/n:0)}'

# Remove these individuals, <5x mean depth across genome
SRR19149587	MABO301
SRR19149582	NABO401
SRR19149570	BFBO506
SRR19149568	PEBO602
SRR19149557 BRBO204 # really high stdev across chromosomes

# Identify sex of individuals based on bamfiles pre-filtering
PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
OUTDIR=${PROJDIR}/datafiles/bamstats/sexing_by_depth
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
while read -r bird; do
samtools idxstats ${BAMDIR}/"$bird".sorted.bam | grep 'CM' > ${OUTDIR}/"$bird".idxstats
done <  ${PROJDIR}/raw_sequences/filenames_SRA.mindepth5x.txt 


mv SRR19149590.idxstats RFBO101.idxstats
mv SRR19149589.idxstats RFBO102.idxstats
mv SRR19149578.idxstats RFBO103.idxstats
mv SRR19149567.idxstats RFBO104.idxstats
mv SRR19149562.idxstats RFBO105.idxstats
mv SRR19149561.idxstats RFBO106.idxstats
mv SRR19149560.idxstats BRBO201.idxstats
mv SRR19149559.idxstats BRBO202.idxstats
mv SRR19149558.idxstats BRBO203.idxstats
mv SRR19149557.idxstats BRBO204.idxstats
mv SRR19149588.idxstats BRBO205.idxstats
mv SRR19149587.idxstats MABO301.idxstats
mv SRR19149586.idxstats MABO302.idxstats
mv SRR19149585.idxstats MABO304.idxstats
mv SRR19149584.idxstats MABO305.idxstats
mv SRR19149583.idxstats MABO306.idxstats
mv SRR19149582.idxstats NABO401.idxstats
mv SRR19149581.idxstats NABO402.idxstats
mv SRR19149580.idxstats NABO403.idxstats
mv SRR19149579.idxstats NABO404.idxstats
mv SRR19149577.idxstats NABO405.idxstats
mv SRR19149576.idxstats NABO406.idxstats
mv SRR19149575.idxstats BFBO501.idxstats
mv SRR19149574.idxstats BFBO502.idxstats
mv SRR19149573.idxstats BFBO503.idxstats
mv SRR19149572.idxstats BFBO504.idxstats
mv SRR19149571.idxstats BFBO505.idxstats
mv SRR19149570.idxstats BFBO506.idxstats
mv SRR19149569.idxstats PEBO601.idxstats
mv SRR19149568.idxstats PEBO602.idxstats
mv SRR19149566.idxstats PEBO603.idxstats
mv SRR19149565.idxstats PEBO604.idxstats
mv SRR19149564.idxstats PEBO605.idxstats
mv SRR19149563.idxstats PEBO606.idxstats


ls "$PWD"/*idxstats > ${PROJDIR}/referencelists/idxindex.txt

R
library(tidyverse)

# list of idxstats files
file_list <- readLines("/xdisk/mcnew/dannyjackson/sulidae/referencelists/idxindex.txt")

# read and combine into one dataframe
all_data <- file_list %>%
  set_names(basename(.)) %>%   # name each dataset by filename
  map_dfr(~ read.table(.x, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
          .id = "individual")

# add column names
colnames(all_data) <- c("individual", "chrom", "length", "mapped", "unmapped")


# NORMALIZE DEPTH
# compute per-chrom depth
all_data <- all_data %>%
  mutate(depth = mapped / length)

# calculate normalization factors for each individual
norm_factors <- all_data %>%
  filter(chrom %in% c("CM062567.1", "CM062568.1", "CM062569.1")) %>%
  group_by(individual) %>%
  summarize(norm_factor = mean(depth, na.rm = TRUE))

# join and normalize
all_data <- all_data %>%
  left_join(norm_factors, by = "individual") %>%
  mutate(norm_depth = depth / norm_factor)

# for ordering chromosomes
all_data <- all_data %>%
  mutate(chrom = factor(chrom, levels = unique(chrom)))

keep_chroms <- c("CM062567.1", "CM062568.1", "CM062569.1", "CM062600.1")

all_data <- all_data %>%
  filter(chrom %in% keep_chroms) %>%
  mutate(chrom = factor(chrom, levels = keep_chroms))

# plot
ggplot(all_data, aes(x = chrom, y = norm_depth, color = individual)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Chromosome", y = "Normalized depth")

ggsave("normalized_depth_plot.png", width = 12, height = 6)


# assign sex
sex_calls <- all_data %>%
  filter(chrom == "CM062600.1") %>%
  select(individual, chrom, norm_depth) %>%
  mutate(sex = ifelse(norm_depth > 0.9, "male", "female"))

# add sex back to all data
all_data <- all_data %>%
  left_join(sex_calls %>% select(individual, sex), by = "individual")

# plot chromosomes colored by sex
ggplot(all_data, aes(x = chrom, y = norm_depth, color = sex)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  scale_color_manual(values = c("male" = "blue", "female" = "red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Chromosome", y = "Normalized depth", color = "Sex")

ggsave("normalized_depth_by_sex.png", width = 12, height = 6)

# clean up names (remove .idxstats)
sex_calls <- all_data %>%
  filter(chrom == "CM062600.1") %>%
  select(individual, chrom, norm_depth) %>%
  mutate(individual = sub("\\.idxstats$", "", individual),
         sex = ifelse(norm_depth > 0.9, "male", "female"))

# list individuals
males <- sex_calls %>% filter(sex == "male") %>% pull(individual)
females <- sex_calls %>% filter(sex == "female") %>% pull(individual)

males
females

> males
 [1] "BFBO501" "BFBO502" "BFBO503" "BFBO505" "BRBO203" "MABO302" "MABO305"
 [8] "NABO404" "PEBO606" "RFBO101" "RFBO103" "RFBO104"
> females
 [1] "BFBO504" "BRBO201" "BRBO202" "BRBO205" "MABO304" "MABO306" "NABO402"
 [8] "NABO403" "NABO405" "NABO406" "PEBO601" "PEBO603" "PEBO604" "PEBO605"
[15] "RFBO102" "RFBO105" "RFBO106"
> 

mv SRR19149590_depthstats.txt RFBO101_depthstats.txt
mv SRR19149589_depthstats.txt RFBO102_depthstats.txt
mv SRR19149578_depthstats.txt RFBO103_depthstats.txt
mv SRR19149567_depthstats.txt RFBO104_depthstats.txt
mv SRR19149562_depthstats.txt RFBO105_depthstats.txt
mv SRR19149561_depthstats.txt RFBO106_depthstats.txt
mv SRR19149560_depthstats.txt BRBO201_depthstats.txt
mv SRR19149559_depthstats.txt BRBO202_depthstats.txt
mv SRR19149558_depthstats.txt BRBO203_depthstats.txt
mv SRR19149557_depthstats.txt BRBO204_depthstats.txt
mv SRR19149588_depthstats.txt BRBO205_depthstats.txt
mv SRR19149587_depthstats.txt MABO301_depthstats.txt
mv SRR19149586_depthstats.txt MABO302_depthstats.txt
mv SRR19149585_depthstats.txt MABO304_depthstats.txt
mv SRR19149584_depthstats.txt MABO305_depthstats.txt
mv SRR19149583_depthstats.txt MABO306_depthstats.txt
mv SRR19149582_depthstats.txt NABO401_depthstats.txt
mv SRR19149581_depthstats.txt NABO402_depthstats.txt
mv SRR19149580_depthstats.txt NABO403_depthstats.txt
mv SRR19149579_depthstats.txt NABO404_depthstats.txt
mv SRR19149577_depthstats.txt NABO405_depthstats.txt
mv SRR19149576_depthstats.txt NABO406_depthstats.txt
mv SRR19149575_depthstats.txt BFBO501_depthstats.txt
mv SRR19149574_depthstats.txt BFBO502_depthstats.txt
mv SRR19149573_depthstats.txt BFBO503_depthstats.txt
mv SRR19149572_depthstats.txt BFBO504_depthstats.txt
mv SRR19149571_depthstats.txt BFBO505_depthstats.txt
mv SRR19149570_depthstats.txt BFBO506_depthstats.txt
mv SRR19149569_depthstats.txt PEBO601_depthstats.txt
mv SRR19149568_depthstats.txt PEBO602_depthstats.txt
mv SRR19149566_depthstats.txt PEBO603_depthstats.txt
mv SRR19149565_depthstats.txt PEBO604_depthstats.txt
mv SRR19149564_depthstats.txt PEBO605_depthstats.txt
mv SRR19149563_depthstats.txt PEBO606_depthstats.txt

# Plot Z scores along Z chromosome for each individual


#!/usr/bin/env bash
#SBATCH --account=mcnew
#SBATCH --job-name=plotZ
#SBATCH --partition=standard
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm_output/plotZ.%A_%a.out
#SBATCH --mail-type=ALL

set -euo pipefail
module load R
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/female_depth_Zchrom.txt                  # each line is a full path to an input file
CHROM="CM062600.1"
INFILE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")"
BASENAME="$(basename "$INFILE" | sed 's/\.[^.]*$//')"
OUTPNG="plots/${BASENAME}_${CHROM}_Z.png"
OUTFILE="outfiles/${BASENAME}_${CHROM}_Z.tsv"

if [[ -z "${INFILE}" ]]; then
  echo "ERROR: Empty INFILE for task ${SLURM_ARRAY_TASK_ID}"; exit 2
fi
if [[ ! -s "${INFILE}" ]]; then
  echo "ERROR: File not found or empty: ${INFILE}"; exit 2
fi

mkdir -p plots slurm_output outfiles

Rscript analyze_depth_Z.r "${INFILE}" "${OUTPNG}" "${CHROM}" "${OUTFILE}"


srun Rscript - "${INFILE}" "${OUTPNG}" "${CHROM}" "${OUTFILE}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outpng <- args[2]; chrom <- args[3]; outfile <- args[4]

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(readr)
})

# Files have NO header and 3 cols: chrom, position, depth (tab/space)
dt <- fread(infile,
            header = FALSE,
            col.names = c("chrom","position","depth"),
            sep = "\t",    # if some are space-separated, fread will auto-handle; this is safe
            data.table = TRUE,
            nThread = getDTthreads(),
            showProgress = FALSE)

# Sanity
stopifnot(all(c("chrom","position","depth") %in% names(dt)))
if (!("position" %in% names(dt)) || nrow(dt) == 0) stop("Empty or malformed file: ", infile)

# Compute Z from depth within file (center/scale)
# (If you already have a z column elsewhere, skip this and read it instead.)
mu <- mean(dt$depth, na.rm = TRUE)
sdv <- sd(dt$depth,  na.rm = TRUE)
if (is.na(sdv) || sdv == 0) {
  # Fall back to mean-centered if variance is zero
  dt[, z := depth - mu]
} else {
  dt[, z := (depth - mu) / sdv]
}

# Optional: light binning for huge base-level series (smoother + memory friendly)
if (nrow(dt) > 5e6) {
  dt[, bin := position %/% 1000L]
  dt <- dt[, .(position = as.integer(1000L * bin),
               z = mean(z, na.rm = TRUE)), by = bin][, bin := NULL]
}

png(outpng, width = 2200, height = 900, res = 150)
print(
  ggplot(dt, aes(position, z)) +
    geom_point( alpha = 0.8, size = 1) +
    labs(x = "Position", y = "Z", title = chrom) +
    theme_minimal(base_size = 12)
)
dev.off()

write_delim(dt, outfile, quote = "none", delim = "\t")

RS








LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/female_depth_Zchrom.txt
N=$(wc -l < $LIST)

# submit with literal array size
sbatch --array=1-"$N" submit.plot_z.sh
# save as list_files.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BFBO504.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BRBO201.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BRBO202.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BRBO205.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/MABO304.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/MABO306.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/NABO402.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/NABO403.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/NABO405.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/NABO406.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/PEBO601.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/PEBO603.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/PEBO604.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/PEBO605.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/RFBO102.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/RFBO105.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/RFBO106.Zdepth.txt

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/male_depth_Zchrom.txt
N=$(wc -l < $LIST)

# submit with literal array size
sbatch --array=1-"$N" submit.plot_z.sh
# males

/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BFBO501.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BFBO502.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BFBO503.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BFBO505.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/BRBO203.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/MABO302.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/MABO305.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/NABO404.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/PEBO606.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/RFBO101.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/RFBO103.Zdepth.txt
/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/RFBO104.Zdepth.txt



# Compute average scaled depth for 1kb windows for males and females to compare

#!/bin/bash
set -euo pipefail

MALE_DIR="/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/outfiles/males"
FEMALE_DIR="/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth/outfiles/females"
OUTFILE="male_female_means.tsv"

male_files=(${MALE_DIR}/*.tsv)
female_files=(${FEMALE_DIR}/*.tsv)

# Compute mean male values
paste "${male_files[@]}" | \
  awk -v OFS="\t" '
  NR==1 { next }
  {
    pos=$1
    sum=0; count=0
    for(i=2;i<=NF;i+=2){
      sum+=$i; count++
    }
    mean_male=sum/count
    print pos, mean_male
  }' > tmp.males

# Compute mean female values
paste "${female_files[@]}" | \
  awk -v OFS="\t" '
  NR==1 { next }
  {
    pos=$1
    sum=0; count=0
    for(i=2;i<=NF;i+=2){
      sum+=$i; count++
    }
    mean_female=sum/count
    print pos, mean_female
  }' > tmp.females

# Join and add ratio, handling divide-by-zero
join -1 1 -2 1 <(sort -k1,1 tmp.males) <(sort -k1,1 tmp.females) \
  | awk -v OFS="\t" 'BEGIN{print "position","mean_male","mean_female","ratio_male_female","ratio_female_male"}
                     {
                       if ($3==0) ratio_MF="NA";
                       else ratio_MF=($2/$3);
                       if ($2==0) ratio_FM="NA";
                       else ratio_FM=($3/$2);
                       print $1,$2,$3,ratio_MF,ratio_FM
                     }' \
  > "$OUTFILE"

rm tmp.males tmp.females


# Plot 
Rscript - "male_female_means.tsv" "male_female_means.png" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outpng <- args[2]

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(readr)
})

# Files have NO header and 3 cols: chrom, position, depth (tab/space)
dt <- fread(infile,
            header = TRUE,
            sep = "\t",    # if some are space-separated, fread will auto-handle; this is safe
            data.table = TRUE,
            nThread = getDTthreads(),
            showProgress = FALSE)

png(outpng, width = 2200, height = 900, res = 150)
print(
  ggplot(dt, aes(position, ratio_male_female)) +
    geom_point( alpha = 0.8, size = 1) +
    labs(x = "Position", y = "Male:Female avg normalized depth") +
    theme_minimal(base_size = 12)
)
dev.off()
RS


# Plot the inverse

# Plot 
Rscript - "male_female_means.tsv" "female_male_means.png" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outpng <- args[2]

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(readr)
})

# Files have NO header and 3 cols: chrom, position, depth (tab/space)
dt <- fread(infile,
            header = TRUE,
            sep = "\t",    # if some are space-separated, fread will auto-handle; this is safe
            data.table = TRUE,
            nThread = getDTthreads(),
            showProgress = FALSE)

png(outpng, width = 2200, height = 900, res = 150)
print(
  ggplot(dt, aes(position, ratio_female_male)) +
    geom_point( alpha = 0.8, size = 1) +
    labs(x = "Position", y = "Female:Male avg normalized depth") +
    theme_minimal(base_size = 12)
)
dev.off()
RS

# plot just avg female

Rscript - "male_female_means.tsv" "female_means.png" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outpng <- args[2]

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(readr)
})

# Files have NO header and 3 cols: chrom, position, depth (tab/space)
dt <- fread(infile,
            header = TRUE,
            sep = "\t",    # if some are space-separated, fread will auto-handle; this is safe
            data.table = TRUE,
            nThread = getDTthreads(),
            showProgress = FALSE)

png(outpng, width = 2200, height = 900, res = 150)
print(
  ggplot(dt, aes(position, mean_female)) +
    geom_point( alpha = 0.8, size = 1) +
    labs(x = "Position", y = "Female avg normalized depth") +
    theme_minimal(base_size = 12)
)
dev.off()
RS


# plot just avg male

Rscript - "male_female_means.tsv" "male_means.png" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outpng <- args[2]

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(readr)
})

# Files have NO header and 3 cols: chrom, position, depth (tab/space)
dt <- fread(infile,
            header = TRUE,
            sep = "\t",    # if some are space-separated, fread will auto-handle; this is safe
            data.table = TRUE,
            nThread = getDTthreads(),
            showProgress = FALSE)

png(outpng, width = 2200, height = 900, res = 150)
print(
  ggplot(dt, aes(position, mean_male)) +
    geom_point( alpha = 0.8, size = 1) +
    labs(x = "Position", y = "Male avg normalized depth") +
    theme_minimal(base_size = 12)
)
dev.off()
RS

awk 'NR==FNR {map[$1]=$2; next} {for (i in map) gsub(i,map[i])}1' /xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt depthstats.parallel.tsv > depthstats.parallel.boobynames.tsv
