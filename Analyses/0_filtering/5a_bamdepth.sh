#!/bin/bash
#SBATCH --job-name=bamstats
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/bamstats.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
BAMFILE=${BAMDIR}/${SAMPLE}.final.bam
OUTFILE=${OUTDIR}/${SAMPLE}_depthstats.txt
STATSFILE=${OUTDIR}/alldepthstats.txt
module load samtools
samtools depth "${BAMDIR}/${SAMPLE}.final.bam" >> "$OUTFILE"

echo $SAMPLE >> ${STATSFILE}

# average and standard deviation
awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' ${OUTFILE} >> ${STATSFILE}

while read -r bird;
do 
  echo $bird >> depthstats.txt

  # average and standard deviation
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$bird"_depthstats.txt >> depthstats.txt

done < ${LIST}