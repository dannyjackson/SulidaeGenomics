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
