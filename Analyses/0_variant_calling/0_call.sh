#!/usr/bin/env bash
#SBATCH --job-name=call
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/call.%A_%a.out
#SBATCH --mail-type=ALL

BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
LIST=$BASE/referencelists/allsamplecodes.txt
BAMDIR=$BASE/datafiles/finalbams
OUTDIR=$BASE/analyses/raxml_wg/

module load bcftools/1.19
module load vcftools
module load plink

ID="Sula_MorusBassanus_2"
bcftools mpileup -Ou -f "$REF" -a FORMAT/AD,DP,INFO/AD,SP "$BAMDIR"/*final.bam | bcftools call -mv -V indels > $OUTDIR/vcfs/"$ID"_snps_multiallelic.vcf