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