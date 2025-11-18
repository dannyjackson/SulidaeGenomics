# chrZ heterozygosity unfiltered

# Make vcf

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
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
LIST=$BASE/referencelists/allsamplecodes.txt
BAMDIR=$BASE/datafiles/bamfiles/
OUTDIR=$BASE/datafiles/bamstats/sexing_by_heterozygosity_filtered/unfiltered_bams/

module load bcftools/1.19
module load vcftools
module load plink

ID="Sula_MorusBassanus_unfilteredbams.CM062595.1"
bcftools mpileup -Ou -f "$REF" -r CM062595.1 -a FORMAT/AD,DP,INFO/AD,SP "$BAMDIR"/*.bam | bcftools call -mv -V indels > $OUTDIR/vcfs/"$ID".vcf


cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered/unfiltered_bams/chr29
# inputs
VCF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered/unfiltered_bams/vcfs/Sula_MorusBassanus_unfilteredbams.CM062595.1.vcf        # bgzipped + indexed VCF
# VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf
CHR=CM062595.1       # chromosome/scaffold name
WIN=50000            # window size
STEP=50000           # step size (non-overlapping)

# make sure VCF is indexed
# bgzip "$VCF"
# tabix -p vcf "$VCF.gz"

# compute π in windows
vcftools --vcf "$VCF" --chr "$CHR" \
         --window-pi $WIN --window-pi-step $STEP \
         --out ${CHR}.pi_${WIN}



library(readr); library(dplyr); library(ggplot2)

win <- read_tsv("CM062595.1.pi_50000.windowed.pi",
                show_col_types = FALSE,
                col_types = cols(
                  CHROM = col_character(),
                  BIN_START = col_double(),
                  BIN_END = col_double(),
                  N_VARIANTS = col_double(),
                  PI = col_double()
                ))

ggplot(win, aes(x = (BIN_START + BIN_END)/2, y = PI)) +
  geom_line(linewidth = 0.4) +
  labs(title = "Windowed nucleotide diversity (π)",
       x = "Position (bp)", y = "π per 50 kb") +
  theme_bw()


ggsave("pi.CM062595.1.allindv.png", width = 12, height = 3, dpi = 300)


### Per individual heterozygosity
# --- inputs you edit ---
CHR=CM062595.1                    # chromosome/scaffold to plot (change as needed)
OUT=hets_${CHR}                   # prefix for intermediates
VCF="$VCF.gz"
# index if needed
tabix -p vcf "$VCF"

# keep biallelic SNPs on that chromosome (clean het calling)
bcftools view -r "$CHR" -v snps -m2 -M2 -O z -o ${OUT}.snps.vcf.gz "$VCF"
tabix -p vcf ${OUT}.snps.vcf.gz

# sample order (matches columns in the query below)
bcftools query -l ${OUT}.snps.vcf.gz > ${OUT}.samples.txt

# CHROM POS then one GT per sample (e.g., 0/0, 0/1, 1|0, 1/1, ./.)
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ${OUT}.snps.vcf.gz > ${OUT}.genotypes.tsv


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

