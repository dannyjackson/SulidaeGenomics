#!/usr/bin/env bash
#SBATCH --job-name=filtervcf
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/filtervcf%A_%a.out
#SBATCH --mail-type=ALL

module load bcftools

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/MABO

bcftools view \
  -r $(paste -sd, /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt) \
  -s MABO302,MABO304,MABO305,MABO306,RFBO101,RFBO102,RFBO103,RFBO104,RFBO105,RFBO106 \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered.recode.vcf.gz \
  -Oz -o RFBO_MABO.autosomes.vcf.gz