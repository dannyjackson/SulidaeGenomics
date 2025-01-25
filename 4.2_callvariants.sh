# call variants
#!/bin/bash
module load bcftools/1.19

ref="/xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna"
bamdir="/xdisk/mcnew/dannyjackson/sulidae/indelrealignment/"
ID="sula"

cd /xdisk/mcnew/dannyjackson/sulidae/genotype_calls

bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.final.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


sbatch --account=mcnew \
--job-name=callvariants \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.callvariants.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=240:00:00 \
callvariants.sh
# Submitted batch job 12032173



#filter by quality
#!/bin/bash

module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module spider samtools/1.19.2

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/genotype_calls

bcftools view -i 'QUAL>100' sula_snps_multiallelic.vcf > sula_qualitysort.vcf

#filters by depth and removes indels
vcftools --vcf sula_qualitysort.vcf --min-meanDP 2 --max-meanDP 8 --remove-indels --recode --out sula_filtered


sbatch --account=mcnew \
--job-name=filtervcf \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.filtervcf.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=40:00:00 \
filtervcf.sh

# Submitted batch job 3503159



#!/bin/bash

module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module spider samtools/1.19.2

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/genotype_calls

plink --vcf sula_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out sula_filtered_mind2

cp sula_filtered_mind2.vcf sula_filtered_cleaned_zip.vcf

module load samtools 
bgzip sula_filtered_cleaned_zip.vcf

bcftools index sula_filtered_cleaned_zip.vcf.gz


sbatch --account=mcnew \
--job-name=filtervcf \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.filtervcf.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=40:00:00 \
filtervcf_2.sh

Submitted batch job 12039087