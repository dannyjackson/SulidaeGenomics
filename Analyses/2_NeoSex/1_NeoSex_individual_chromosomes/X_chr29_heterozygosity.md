# compare heterozygosity across chr29

Set up environment and directory structure.

```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/

mkdir hetplot_chr29
cd hetplot_chr29

module load micromamba vcftools
source ~/.bashrc
micromamba activate r_ocelote
```

```
# inputs
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/snpdensity/CM062595.1.qualitysort_filtered_mind2.vcf        # bgzipped + indexed VCF
# VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf
CHR=CM062595.1       # chromosome/scaffold name
WIN=1000000            # window size
STEP=1000000           # step size (non-overlapping)

# make sure VCF is indexed
# bgzip "$VCF"
# tabix -p vcf "$VCF.gz"

# compute π in windows
vcftools --gzvcf "$VCF.gz" --chr "$CHR" \
         --window-pi $WIN --window-pi-step $STEP \
         --out ${CHR}.pi_${WIN}

# Output: ${CHR}.pi_${WIN}.windowed.pi
# columns: CHROM  BIN_START  BIN_END  N_SNPs  PI

Rscript chr29_heterozygosity.1.r


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


# explore genes in highly differentiated regions
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff
```
Plot all plots together
```
Rscript chr29_heterozygosity.3.r
```
