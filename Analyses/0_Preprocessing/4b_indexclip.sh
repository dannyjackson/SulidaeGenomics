#!/bin/bash
#SBATCH --job-name=indexclip
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/indexclip.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
SORTBAM=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/
OUTDIR=${SORTBAM}/indelmaps

mkdir -p "$OUTDIR" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')

BAMFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.clip.bam

module load samtools
samtools index ${BAMFILE} -@ ${SLURM_CPUS_PER_TASK}