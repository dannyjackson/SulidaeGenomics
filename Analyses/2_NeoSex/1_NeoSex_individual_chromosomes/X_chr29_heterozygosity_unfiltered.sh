# chrZ heterozygosity unfiltered

# Make vcf








#!/usr/bin/env bash
#SBATCH --job-name=call
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/call.%A_%a.out
#SBATCH --mail-type=ALL

BASE=/xdisk/mcnew/dannyjackson/sulidae
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
LIST=$BASE/referencelists/allsamplecodes.txt
BAMDIR=$BASE/datafiles/bamfiles/
OUTDIR=$BASE/datafiles/bamstats/sexing_by_heterozygosity_filtered/unfiltered_bams/

module load bcftools/1.19
module load vcftools
module load plink

ID="Sula_MorusBassanus_unfilteredbams.CM062600.1"
bcftools mpileup -Ou -f "$REF" -r CM062600.1 -a FORMAT/AD,DP,INFO/AD,SP "$BAMDIR"/*.bam | bcftools call -mv -V indels > $OUTDIR/vcfs/"$ID".vcf

