#!/usr/bin/env bash
#SBATCH --job-name=all_dsuite
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --output=slurm_output/all_dsuite.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch all_dsuite.sh

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods

module load bcftools htslib

bgzip /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered.recode.vcf
tabix -p vcf /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered.recode.vcf.gz

echo 'filtering vcf to autosomes only'
bcftools view \
  -r $(paste -sd, /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt) \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered.recode.vcf.gz \
  -Oz -o Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

echo 'performing dsuite analysis'

VCF=Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t all.nwk --ABBAclustering -g