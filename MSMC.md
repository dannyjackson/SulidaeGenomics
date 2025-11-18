# run MSMC markdown

msmc2 -I 0,1 -o indivA_allchr out/chr*.multihetsep.txt

INDMASK=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/ind_mask.BFBO501.CM062567.1.bed.gz
REFMASK=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062567.1.mask.150.50.bed.gz
VCF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/phasedvcfs/BFBO501.CM062567.1.whatshap.samtools.vcf.gz
OUTPUT=/xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/BFBO501.CM062567.1.input.txt

~/programs/msmc-tools/generate_multihetsep.py \
    --mask=${INDMASK} \
    --mask=${REFMASK} \
    ${VCF} \
    > ${OUTPUT}

#!/usr/bin/env bash
#SBATCH --job-name=MSMC_input
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/MSMC_input.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-29 MSMC_input.sh /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

set -euo pipefail

VCFDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/phasedvcfs
MASKDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
IND="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

mkdir -p "$OUTDIR"

shopt -s nullglob
for vcf in "${VCFDIR}"/${IND}*.whatshap.samtools.vcf.gz; do
    base=$(basename "$vcf")   # e.g. BFBO501.CM062567.1.whatshap.samtools.vcf.gz
    chrom=${base#*.}          # CM062567.1.whatshap.samtools.vcf.gz
    chrom=${chrom%%.*}        # everything before next dot → CM062567.1

    INDMASK=${MASKDIR}/ind_mask.${IND}.${chrom}.1.bed.gz
    REFMASK=${MASKDIR}/GCA_031468815_${chrom}.1.mask.150.50.bed.gz
    OUTPUT=${OUTDIR}/${IND}.${chrom}.1.input.txt

    echo "Processing $IND $chrom"

    ~/programs/msmc-tools/generate_multihetsep.py \
        --mask="$INDMASK" \
        --mask="$REFMASK" \
        "$vcf" \
        > "$OUTPUT"
done



~/programs/msmc_2.0.0_linux64bit -t 16 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT
# msmc -i 20 -t 8 -p 1*2+25*1+1*2 -o indivA_allchr out/chr*.multihetsep.txt
# msmc2 -I 0,1 -o indivA_allchr out/chr*.multihetsep.tx

#!/usr/bin/env bash
#SBATCH --job-name=MSMC_run
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=6
#SBATCH --mem=200G
#SBATCH --time=120:00:00
#SBATCH --output=slurm_output/MSMC_run.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-29 MSMC_run.autosomes.sh /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

INDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/autosomes/
IND="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

~/programs/msmc_2.0.0_linux64bit -I 0,1 -o ${OUTDIR}/${IND}_allchr ${INDIR}/${IND}*.input.txt


# Plot all sample MSMC results together
Rscript plotMSMC.allsamples.R

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

# Parameters
mu <- 2.3e-09
gen <- 7

# Directory of input files
indir <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/autosomes/"
files <- list.files(indir, pattern="*.final.txt", full.names=TRUE)

# Extract data
all_data <- lapply(files, function(f) {
  dat <- read.csv(f, sep="\t")
  
  species <- str_extract(basename(f), "^[A-Z]+")  # grabs BFBO, PEBO, etc.
  indiv <- str_extract(basename(f), "^[A-Z]+[0-9]+")

  time <- (dat$left_time_boundary/mu*gen)
  pop.size <- (1/dat$lambda)/(2*mu)
  
  data.frame(time=time, pop.size=pop.size, species=species, individual=indiv)
})

all_data <- bind_rows(all_data)


# Load color codes
colorcodes <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/colorcodes.tsv")

# Make a named vector for scale_color_manual()
species_colors <- setNames(colorcodes$Code, colorcodes$Species)

# Plot
current_date <- format(Sys.Date(), "%m%d%y")
outfile <- paste0("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots/MSMC_all_species_", current_date, ".autosomes.pdf")

pdf(outfile, width = 16, height = 6)  # higher resolution and larger canvas

ggplot(all_data, aes(x = time, y = pop.size, color = species, group = individual)) +
  geom_step(linewidth = 1.4, alpha = 0.75) +  # thicker lines
  scale_x_log10(limits = c(1e4, 1e7)) +
  ylim(0, 500000) +
  scale_color_manual(values = species_colors) +  # apply your color scheme
  labs(
    x = "Years before present (log scale)",
    y = "Effective Population Size",
    title = "MSMC across individuals and species"
  ) +
  theme_minimal(base_size = 22) +  # larger base font size
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.6),
    legend.position = "null"
  )

dev.off()



# Make a facet plot with one panel per species
suppressMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# Output dir + date stamp
outdir  <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date(), "%m%d%y")

# Unique species
species_list <- sort(unique(all_data$species))

# Save one PNG per species
for (sp in species_list) {
  df <- filter(all_data, species == sp)

  if (nrow(df) == 0) next

  p_sp <- ggplot(df, aes(x = time, y = pop.size, color = individual, group = individual)) +
  geom_step(linewidth = 1.4, alpha = 0.75) +  # thicker lines
  scale_x_log10(limits = c(1e4, 1e7)) +
  ylim(0, 200000) +
  labs(
    x = "Years before present (log scale)",
    y = "Effective Population Size",
    title = paste0("MSMC per individual — ", sp),
    color = "Individual"
  ) +
  theme_minimal(base_size = 22) +  # larger base font size
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.6),
  )

  outfile <- file.path(outdir, paste0("MSMC_", sp, "_", stamp, ".pdf"))
  ggsave(outfile, p_sp, width = 10, height = 7, dpi = 300)
  message("Saved: ", outfile)
}

Rscript plotMSMC.allsamples.facet.R



# Can we make a single plot for multiple individuals?

# make list of files of all individuals

#!/usr/bin/env bash
#SBATCH --job-name=MSMC_multihetsep
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/MSMC_multihetsep.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-6 MSMC_multihetsep.sh
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt
IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
MASK_GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062567.1.mask.150.50.bed.gz
VCFLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/vcfprocessing/${IND}_phasedvcfs.txt
OUTFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species/${IND}.all.txt

mkdir -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species/

ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/phasedvcfs/renamedvcfs/${IND}*vcf.gz > ${VCFLIST}

~/programs/msmc-tools/generate_multihetsep.py --mask=$MASK_GENOME `cat ${VCFLIST}` > ${OUTFILE}





# different attempt

#!/usr/bin/env bash
#SBATCH --job-name=MSMC_multihetsep_chrom
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/MSMC_multihetsep_chrom.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-6 MSMC_multihetsep_chrom.sh

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt
IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
MASK_GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062567.1.mask.150.50.bed.gz

for SCAFFOLD in `cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.txt`
    do 
    VCFLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/${IND}_${SCAFFOLD}_phasedvcfs.txt
    OUTFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/${IND}_${SCAFFOLD}.txt
    MASKLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/${IND}_${SCAFFOLD}_masks.txt
    
    mkdir -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/

    ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/phasedvcfs/renamedvcfs/${IND}*${SCAFFOLD}*vcf.gz > ${VCFLIST}
    ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/*${IND}*${SCAFFOLD}.bed.gz > ${MASKLIST}

    ~/programs/msmc-tools/generate_multihetsep.py `cat ${MASKLIST} | sed 's/^/--mask=/'` \
        --mask=$MASK_GENOME `cat ${VCFLIST}` > ${OUTFILE}
done

# check status 
ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/BFBO* | wc -l


#!/usr/bin/env bash
#SBATCH --job-name=MSMC_run_species
#SBATCH --partition=standard
#SBATCH --constraint=hi_mem
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=32
#SBATCH --mem=3008G
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/MSMC_run_species.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-6 MSMC_run_species.sh
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt
IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
MSMC_OUTPUT=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/${IND}_allchr
NR_IND=` grep ${IND} /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt | wc -l `
MSMC_INPUT_FILE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/${IND}_msmc_input.txt
ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/${IND}*.txt  > ${MSMC_INPUT_FILE}
MSMC_INPUT=`cat ${MSMC_INPUT_FILE}`
n=$(expr ${NR_IND} + ${NR_IND} - 2)
INDEX=$(for num in `seq 0 ${n}`; do echo -n "${num},"; done; echo $(expr ${NR_IND} + ${NR_IND} - 1))

~/programs/msmc_2.0.0_linux64bit -t 32 -p 1*2+15*1+1*2 -i 100 -o /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/${IND}_allchr -I `echo $INDEX` ${MSMC_INPUT}


# Rerun RFBO 
#!/usr/bin/env bash
#SBATCH --job-name=MSMC_run_RFBO
#SBATCH --partition=standard
#SBATCH --constraint=hi_mem
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=94
#SBATCH --mem=3008G
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/MSMC_run_RFBO.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch MSMC_run_RFBO.sh
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt
IND=RFBO
MSMC_OUTPUT=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/${IND}_allchr_solorun
NR_IND=` grep ${IND} /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt | wc -l `
MSMC_INPUT_FILE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/${IND}_msmc_input.txt
ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/${IND}*.txt  > ${MSMC_INPUT_FILE}
MSMC_INPUT=`cat ${MSMC_INPUT_FILE}`
n=$(expr ${NR_IND} + ${NR_IND} - 2)
INDEX=$(for num in `seq 0 ${n}`; do echo -n "${num},"; done; echo $(expr ${NR_IND} + ${NR_IND} - 1))

~/programs/msmc_2.0.0_linux64bit -t 94 -p 1*2+15*1+1*2 -i 100 -o /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/${IND}_allchr -I `echo $INDEX` ${MSMC_INPUT}




# Plot all sample MSMC results together
Rscript plotMSMC.allsamples.R

library(ggplot2)
library(dplyr)
library(stringr)

# Parameters
mu <- 2.3e-09
gen <- 7

# Directory of input files
indir <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/"
files <- list.files(indir, pattern="*.final.txt", full.names=TRUE)

# Extract data
all_data <- lapply(files, function(f) {
  dat <- read.csv(f, sep="\t")
  
  species <- str_extract(basename(f), "^[A-Z]+")  # grabs BFBO, PEBO, etc.
  indiv <- str_extract(basename(f), "^[A-Z]+[0-9]+")

  time <- (dat$left_time_boundary/mu*gen)
  pop.size <- (1/dat$lambda)/(2*mu)
  
  data.frame(time=time, pop.size=pop.size, species=species, individual=indiv)
})

all_data <- bind_rows(all_data)

# Plot
current_date <- format(Sys.Date(), "%m%d%y")
outfile <- paste0("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots/MSMC_all_fulldata_", current_date, ".png")

png(outfile, width=800, height=600)

ggplot(all_data, aes(x=time, y=pop.size, color=species)) +
  geom_step(alpha=0.7) +
  scale_x_log10(limits = c(1e5, 1e7)) +   # restrict x axis to 10^4–10^7
  ylim(0, 500000) +
  labs(
    x = "Years before present (log scale)",
    y = "Effective Population Size",
    title = "MSMC across individuals and species"
  ) +
  theme_minimal(base_size=14)

dev.off()

outfile <- paste0("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots/MSMC_all_fulldata_zoomed_", current_date, ".png")

png(outfile, width=800, height=600)

ggplot(all_data, aes(x=time, y=pop.size, color=species)) +
  geom_step(alpha=0.7) +
  scale_x_log10(limits = c(1e4, 1e7)) +   # restrict x axis to 10^4–10^7
  ylim(0, 100000) +
  labs(
    x = "Years before present (log scale)",
    y = "Effective Population Size",
    title = "MSMC across individuals and species"
  ) +
  theme_minimal(base_size=14)

dev.off()

# Make a facet plot with one panel per species
suppressMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# Output dir + date stamp
outdir  <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date(), "%m%d%y")

# Unique species
species_list <- sort(unique(all_data$species))

# Save one PNG per species
for (sp in species_list) {
  df <- filter(all_data, species == sp)

  if (nrow(df) == 0) next

  p_sp <- ggplot(df, aes(x = time, y = pop.size, color = individual, group = individual)) +
    geom_step(alpha = 0.9, linewidth = 0.7) +
    scale_x_log10(
      limits = c(1e4, 1e6),
      labels = label_number(scale_cut = cut_si(" "))
    ) +
    coord_cartesian(ylim = c(0, 5e5)) +
    labs(
      x = "Years before present (log scale)",
      y = "Effective population size (Ne)",
      title = paste0("MSMC per individual — ", sp),
      color = "Individual"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  outfile <- file.path(outdir, paste0("MSMC_fulldata_", sp, "_", stamp, ".png"))
  ggsave(outfile, p_sp, width = 10, height = 7, dpi = 300)
  message("Saved: ", outfile)
}

Rscript plotMSMC.allsamples.facet.R


ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/BFBO* | wc -l


# Estimate divergence between species
# BLUE FOOTED AND PERUVIAN
BFBO PEBO
MABO NABO
BFBO PEBO MABO NABO
BFBO PEBO MABO NABO BRBO
BFBO PEBO MABO NABO BRBO RFBO
BRBO RFBO

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/divergencesets.txt
ALL=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

# 1) Get the Nth line (strip CRs/extra whitespace)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r' | xargs)

# 2) Split into an array of species codes (e.g., BFBO PEBO ...)
read -r -a codes <<< "$line"

# 3) Make a label like BFBO_PEBO_MABO
IND=$(IFS=_; printf '%s' "${codes[*]}")

# 4) Build a regex that matches any of the prefixes at line start: ^(BFBO|PEBO|...)
regex="^($(printf '%s|' "${codes[@]}" | sed 's/|$//'))"

# 5) Count individuals whose sample IDs start with any of those prefixes
#    Works whether ALL has "BFBO501" or "BFBO501 <tabs> other stuff"
NR_IND=$(grep -E -c "$regex" "$ALL")

# 6) (Optional) write out the matching sample IDs (first column) for downstream tools
grep -E "$regex" "$ALL" | awk '{print $1}' > "${IND}.samples"

# 7) Your MSMC output stem using the combined label
MSMC_OUTPUT=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/${IND}_allchr

echo "SET:         $line"
echo "IND:         $IND"
echo "NR_IND:      $NR_IND"
echo "SAMPLE LIST: ${IND}.samples"
echo "MSMC_OUTPUT: $MSMC_OUTPUT"

MSMC_INPUT_FILE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/${IND}_msmc_input.txt
ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/msmc_input/species_scaffold/${IND}*.txt  > ${MSMC_INPUT_FILE}
MSMC_INPUT=`cat ${MSMC_INPUT_FILE}`
n=$(expr ${NR_IND} + ${NR_IND} - 2)
INDEX=$(for num in `seq 0 ${n}`; do echo -n "${num},"; done; echo $(expr ${NR_IND} + ${NR_IND} - 1))

~/programs/msmc_2.0.0_linux64bit -t 32 -p 1*2+15*1+1*2 -i 100 -o /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/${IND}_allchr -I `echo $INDEX` ${MSMC_INPUT}

FILES=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species
~/programs/msmc_2.0.0_linux64bit -I `echo $INDEX` -o BFBO_PEBO --skipAmbiguous 1 \
      $FILES/BFBO_chr*.txt $FILES/BRBO_chr*.txt