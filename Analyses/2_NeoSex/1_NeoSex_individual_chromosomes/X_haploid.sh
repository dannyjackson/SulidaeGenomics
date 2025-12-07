# Call Z as haploid
# 1) Make a samples file with sex codes (2 columns, space/tab-separated)
awk -F',' 'BEGIN{OFS="\t"} {code=($2=="Male"?"M":"F"); print $1, code}' sexID.csv > samples.sex

# map back to og file names
awk 'NR==FNR {map[$2]=$1; next} ($1 in map) {print map[$1], $2}' filenameconversion.txt samples.sex > samples.renamed.sex

# 2) Get the contig length and build a valid ploidy file
# make sure REF.fasta.fai exists


BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

LEN29=$(awk -v c=CM062595.1 '$1==c{print $2}' "${REF}.fai")
LENZ=$(awk -v c=CM062600.1 '$1==c{print $2}' "${REF}.fai")

cat > ploidy.txt <<EOF
# CHROM      FROM  TO         SEX PLOIDY
*            *     *          M   2
*            *     *          F   2
CM062595.1   1     $LEN29       M   2
CM062595.1   1     $LEN29       F   1
CM062600.1   1     $LENZ       M   2
CM062600.1   1     $LENZ       F   1
EOF

# Make ped file
cd /xdisk/mcnew/dannyjackson/sulidae/referencelists
R

library(readr)
library(dplyr)

# Read your CSV
sex <- read_csv("sexID.csv", col_names = c("IID", "Sex"))

# Convert to PED format
ped <- sex %>%
  mutate(
    FID = IID,
    PID = 0,
    MID = 0,
    SEX = case_when(
      Sex == "M" ~ 1,
      Sex == "F" ~ 2,
      TRUE ~ 0
    ),
    PHENOTYPE = -9
  ) %>%
  select(FID, IID, PID, MID, SEX, PHENOTYPE)

# Write to .ped file (tab-separated)
write.table(ped, "sexID.ped", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
LIST=$BASE/referencelists/allsamplecodes.txt
BAMDIR=$BASE/datafiles/finalbams
OUTDIR=$BASE/analyses/raxml_wg/
SEX=$BASE/referencelists/samples.renamed.sex
PLOIDY=$BASE/referencelists/ploidy.txt

# call for the neosex
bcftools mpileup -Ou \
  -f "$REF" \
  -a FORMAT/DP,FORMAT/AD,SP \
  -r CM062595.1 \
  "$BAMDIR"/*final.bam \
| bcftools call -mv -V indels \
  --ploidy-file ${PLOIDY} \
  -S ${SEX} \
  -Oz -o Z_haploid.vcf.gz

bcftools index -t 29_haploid.vcf.gz


bcftools reheader -s /xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt \
    -o 29_haploid.vcf.renamed.vcf.gz \
    29_haploid.vcf.gz

bcftools view -i 'QUAL>30' 29_haploid.vcf.renamed.vcf.gz  > 29_haploid.qualitysort.vcf.gz


vcftools --gzvcf 29_haploid.qualitysort.vcf.gz \
    --min-meanDP 2 --max-meanDP 20 --remove-indels --recode \
    --out 29_haploid.qualitysort_filtered


plink --vcf 29_haploid.qualitysort_filtered.recode.vcf \
    --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid \
    --out 29_haploid.qualitysort_filtered_mind2

bgzip 29_haploid.qualitysort_filtered_mind2.vcf
bcftools index 29_haploid.qualitysort_filtered_mind2.vcf.gz


# call for Z

bcftools mpileup -Ou \
  -f "$REF" \
  -a FORMAT/DP,FORMAT/AD,SP \
  -r CM062600.1 \
  "$BAMDIR"/*final.bam \
| bcftools call -mv -V indels \
  --ploidy-file ${PLOIDY} \
  -S ${SEX} \
  -Oz -o Z_haploid.vcf.gz

bcftools index -t Z_haploid.vcf.gz


bcftools reheader -s /xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt \
    -o Z_haploid.vcf.renamed.vcf.gz \
    Z_haploid.vcf.gz

bcftools view -i 'QUAL>30' Z_haploid.vcf.renamed.vcf.gz  > Z_haploid.qualitysort.vcf.gz


vcftools --gzvcf Z_haploid.qualitysort.vcf.gz \
    --min-meanDP 2 --max-meanDP 20 --remove-indels --recode \
    --out Z_haploid.qualitysort_filtered


plink --vcf Z_haploid.qualitysort_filtered.recode.vcf \
    --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid \
    --out Z_haploid.qualitysort_filtered_mind2

bgzip Z_haploid.qualitysort_filtered_mind2.vcf
bcftools index Z_haploid.qualitysort_filtered_mind2.vcf.gz
