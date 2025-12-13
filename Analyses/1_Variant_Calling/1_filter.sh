#!/usr/bin/env bash
#SBATCH --job-name=filter
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/filter.out
#SBATCH --mail-type=ALL

module load vcftools
module load bcftools/1.19
module load plink

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs

grep -v 'JAVH' Sula_MorusBassanus_snps_multiallelic.vcf > Sula_MorusBassanus_snps_multiallelic.contigs.vcf 

bcftools reheader -s /xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt \
    -o Sula_MorusBassanus_snps_multiallelic.contigs.renamed.vcf \
    Sula_MorusBassanus_snps_multiallelic.contigs.vcf 

bcftools view -i 'QUAL>30' Sula_MorusBassanus_snps_multiallelic.contigs.renamed.vcf  > Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf

vcftools --vcf Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf \
    --min-meanDP 4 --max-meanDP 20 --remove-indels --recode \
    --out Sula_MorusBassanus.qualitysort_filtered

plink --vcf Sula_MorusBassanus.qualitysort_filtered.recode.vcf \
    --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid \
    --out Sula_MorusBassanus.qualitysort_filtered_mind2

bgzip Sula_MorusBassanus.qualitysort_filtered_mind2.vcf
bcftools index Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz
