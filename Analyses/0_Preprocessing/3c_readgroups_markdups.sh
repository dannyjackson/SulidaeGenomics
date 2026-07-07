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

picard AddOrReplaceReadGroups \
I=${IN}/${SAMPLE}.sorted.bam \
O=${IN}/${SAMPLE}.sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${SAMPLE}

picard MarkDuplicates \
I=${IN}/${SAMPLE}.sorted_RGadded.bam \
O=${OUT}/${SAMPLE}.sorted_RGadded_dupmarked.bam \
M=${OUT}/metrics/${SAMPLE}.duplicate.metrics.txt

samtools index ${OUT}/${SAMPLE}.sorted_RGadded_dupmarked.bam

