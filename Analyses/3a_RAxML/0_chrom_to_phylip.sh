#!/usr/bin/env bash
#SBATCH --job-name=chromtophylip
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/chromtophylip.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-35 chrom_to_phylip.sh 

set -euo pipefail

INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

module load micromamba
eval "$(micromamba shell hook bash)"

micromamba activate r_ocelote

echo 'vcf2phylip for $CHROM'

python ~/programs/vcf2phylip.py3 \
    -i "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.vcf \
    -o "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.phy \
    -r 

echo 'ascbias for $CHROM'

python ~/programs/ascbias.py -p "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.phy -o "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.noINV.phy

echo 'DONE'