# align and sort 
## download Northern Gannett (Morus bassanus) reference genome 
module load sratoolkit/3.1.1

~/programs/datasets download genome accession GCA_031468815.1 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip 


#!/bin/bash
#SBATCH --job-name=align
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12            # threads per sample
#SBATCH --mem=64G                     # ~5G/core is plenty for bwa-mem2
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/align.%A_%a.out

set -euo pipefail

module load bwa
module load samtools
module list

LIST=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
TRIM=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq
OUT=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$OUT" "$SCRATCH" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
echo "[$(date)] aligning $SAMPLE with ${SLURM_CPUS_PER_TASK} threads"

# Use node-local scratch for sort temp; then move result
TMPDIR="$SCRATCH"

bwa mem -t ${SLURM_CPUS_PER_TASK} /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${SAMPLE}_trimmed_1P.fq.gz \
/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${SAMPLE}_trimmed_2P.fq.gz | samtools sort -@ $((SLURM_CPUS_PER_TASK-2)) -m 3G -T "$SCRATCH/${SAMPLE}.tmp" -o "$OUT/${SAMPLE}.bam" -

samtools index -@ 2 "$OUT/${SAMPLE}.bam"
echo "[$(date)] done $SAMPLE"


sbatch --array=1-34%10 align.sh


#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/sort.%A_%a.out


set -euo pipefail

module load bwa
module load picard
module load samtools
module load parallel
module list

LIST=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt
OUT=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$OUT" "$SCRATCH" slurm_output


SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
echo "[$(date)] aligning $SAMPLE with ${SLURM_CPUS_PER_TASK} threads"

# Use node-local scratch for sort temp; then move result
TMPDIR="$SCRATCH"

samtools sort -@ 2 "$OUT/${SAMPLE}.bam" -o "$OUT/${SAMPLE}.sorted.bam"
samtools index -@ 2 "$OUT/${SAMPLE}.sorted.bam"

echo "[$(date)] done $SAMPLE"

sbatch --array=1-34%10 sort.sh



#!/bin/bash
#SBATCH --job-name=addreadgroups
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/addreadgroups.%A_%a.out

set -euo pipefail

module load picard
module load samtools

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
IN=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
OUT=/xdisk/mcnew/dannyjackson/sulidae/bamfiles
SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$OUT" "$SCRATCH" slurm_output


SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
echo "[$(date)] aligning $SAMPLE with ${SLURM_CPUS_PER_TASK} threads"

# Use node-local scratch for sort temp; then move result
TMPDIR="$SCRATCH"

# picard AddOrReplaceReadGroups \
# I=${IN}/${SAMPLE}.sorted.bam \
# O=${IN}/${SAMPLE}.sorted_RGadded.bam \
# RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${SAMPLE}

picard MarkDuplicates \
I=${IN}/${SAMPLE}.sorted_RGadded.bam \
O=${OUT}/${SAMPLE}.sorted_RGadded_dupmarked.bam \
M=${OUT}/metrics/${SAMPLE}.duplicate.metrics.txt

samtools index ${OUT}/${SAMPLE}.sorted_RGadded_dupmarked.bam

# SUBMIT: sbatch --array=1-34%10 readgroups.sh
# SUBMIT: sbatch --array=1-34%10 markdups.sh

# rename sorted bams using species codes and sample numbers

#!/usr/bin/env bash

conv="/xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt"

while read -r old new; do
  for f in ${old}*; do
    if [[ -e "$f" ]]; then
      newname="${f/$old/$new}"
      mv "$f" "$newname"
    fi
  done
done < "$conv"
