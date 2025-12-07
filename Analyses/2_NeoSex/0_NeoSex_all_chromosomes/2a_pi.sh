#!/usr/bin/env bash
#SBATCH --job-name=pi
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/pi.%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/piplot_allchr

module load vcftools samtools htslib

INDFILE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
CONTIGS=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
VCF_DIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms

while read -r IND; do
while read -r CHROM; do
  VCF="${VCF_DIR}/${CHROM}.qualitysort_filtered_mind2.vcf.gz"

  # per-site π for this chromosome, restricted to IND list (population π)
  vcftools --gzvcf "${VCF}" --chr "${CHROM}" \
           --indv "${IND}" \
           --site-pi \
           --out "${CHROM}.${IND}.pi"
done < "${CONTIGS}"
done < "${INDFILE}"