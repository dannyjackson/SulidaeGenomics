#!/bin/bash
#SBATCH --job-name=clipoverlap
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/clipoverlap.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
SORTBAM=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/
SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$OUTDIR" "$SCRATCH" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
INFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.bam
OUTFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.clip.bam
echo "[$(date)] aligning $SAMPLE with ${SLURM_CPUS_PER_TASK} threads"

~/programs/bamUtil-master/bin/bam clipOverlap \
--in ${INFILE} \
--out ${OUTFILE} \
--stats --params --poolSize 10000000

echo ${SAMPLE} >> /xdisk/mcnew/dannyjackson/sulidae/clippingdone.txt