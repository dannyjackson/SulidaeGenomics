#!/usr/bin/env bash
#SBATCH --job-name=splitvcfchrom
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/splitvcfchrom.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-35 2_split_vcf_chrom.sh 

set -euo pipefail

module load bcftools

INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

bcftools view -r "${CHROM}" "${INDIR}"/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz > "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.vcf
