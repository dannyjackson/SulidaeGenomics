# compare heterozygosity across chr29
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered

mkdir hetplot_chr29
cd hetplot_chr29

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


GFF=chr29.PC.genes.gff

Boundaries:
Diff: 
1100000-21500000
2.75-4.1

awk '$4 >= 1 && $5 <= 1100000' $GFF > region0.PC.diff.gff

awk '$4 >= 1100000 && $5 <= 2150000' $GFF > region1.PC.diff.gff
awk '$4 >= 2750000 && $5 <= 4100000' $GFF > region5.PC.diff.gff


Uncertain:
2.2-2.25
2.15-2.2
awk '$4 >= 2200000 && $5 <= 2250000' $GFF > region2.PC.diff.gff
awk '$4 >= 2150000 && $5 <= 2200000' $GFF > region3.PC.diff.gff


Not diff:
2.25-2.75
4.1-4.7

awk '$4 >= 2250000 && $5 <= 2750000' $GFF > region4.PC.diff.gff
awk '$4 >= 4100000 && $5 <= 4700000' $GFF > region6.PC.diff.gff

# pull out genes
for f in region*.PC.diff.gff; do
  awk -F'\t' '
    $3=="gene" {
      name="NA"; desc="NA";
      if (match($9,/;?Name=([^;]+)/,n))  name=n[1];
      if (match($9,/;?description=([^;]+)/,d)) desc=d[1];
      print name "\t" desc
    }' "$f" > "${f%.gff}.PC.names_desc.tsv"
done


awk -v FS='\t' -v OFS='\t' '
$3=="gene"{
  name="NA"; desc="NA";

  # split attributes on semicolons
  n = split($9, a, ";");
  for (i=1; i<=n; i++) {
    key = a[i]; sub(/=.*/, "", key);          # text before first =
    val = a[i]; sub(/^[^=]*=/, "", val);      # text after first =
    if (key=="Name")        name=val;
    else if (key=="description") desc=val;
  }

  print $1, $4, $5, name, desc;
}' "$GFF" > genecoords.PC.tsv

grep -Ei 'PMVK|FDPS|FLAD1' genecoords.PC.tsv > isoprenoidgenes.PC.tsv
grep 'keratin' genecoords.PC.tsv > keratingenes.PC.tsv


GFF=chr29.MB.genes.gff


Boundaries:
Diff: 
1100000-21500000
2.75-4.1

awk '$4 >= 1 && $5 <= 1100000' $GFF > region0.MB.diff.gff

awk '$4 >= 1100000 && $5 <= 2150000' $GFF > region1.MB.diff.gff
awk '$4 >= 2750000 && $5 <= 4100000' $GFF > region5.MB.diff.gff


Uncertain:
2.2-2.25
2.15-2.2
awk '$4 >= 2200000 && $5 <= 2250000' $GFF > region2.MB.diff.gff
awk '$4 >= 2150000 && $5 <= 2200000' $GFF > region3.MB.diff.gff


Not diff:
2.25-2.75
4.1-4.7

awk '$4 >= 2250000 && $5 <= 2750000' $GFF > region4.MB.diff.gff
awk '$4 >= 4100000 && $5 <= 4700000' $GFF > region6.MB.diff.gff


# pull out genes
for f in region*.MB.diff.gff; do
  awk -F'\t' '
    $3=="gene" {
      name="NA"; desc="NA";
      if (match($9,/;?Name=([^;]+)/,n))  name=n[1];
      if (match($9,/;?description=([^;]+)/,d)) desc=d[1];
      print name "\t" desc
    }' "$f" > "${f%.gff}.MB.names_desc.tsv"
done


awk -v FS='\t' -v OFS='\t' '
$3=="gene"{
  name="NA"; desc="NA";

  # split attributes on semicolons
  n = split($9, a, ";");
  for (i=1; i<=n; i++) {
    key = a[i]; sub(/=.*/, "", key);          # text before first =
    val = a[i]; sub(/^[^=]*=/, "", val);      # text after first =
    if (key=="Name")        name=val;
    else if (key=="description") desc=val;
  }

  print $1, $4, $5, name, desc;
}' "$GFF" > genecoords.MB.tsv

grep -Ei 'PMVK|FDPS|FLAD1' genecoords.MB.tsv > isoprenoidgenes.MB.tsv
grep 'keratin' genecoords.MB.tsv > keratingenes.MB.tsv



Region 0:
Region 1:
  PMVK,FDPS
  Isoprenoid metabolism

Region 2: NA
Region 3: NA
Region 4:
  Neurological stuff, lots of significant GO terms

Region 5:
  No GO terms survive FDR but the list has a shitload of keratin-like LOC unidentified genes
Region 6:
  None survive FDR, contains a variety of terms with p<0.05

Rscript chr29_heterozygosity.2.r

# Plot with genes
# chr 29 genes

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff
grep 'CM062595.1' $GFF | grep 'ID\=gene' > chr29.PC.genes.gff


GFF_RE=/xdisk/mcnew/dannyjackson/sulidae/analyses/repeats/vertebrata/ChrsOfInterest.fa.out.gff
grep 'CM062595.1' $GFF_RE > chr29.RE.gff

GFF_MB=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
grep 'CM062595.1' $GFF_MB | grep 'ID\=gene' > chr29.MB.genes.gff


Rscript chr29_heterozygosity.3.r