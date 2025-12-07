
## 1. Chromosome 29
### Pi
Set up working directory and environment.
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference

mkdir -p piplot_chr29
cd piplot_chr29

module load micromamba vcftools samtools htslib
source ~/.bashrc
micromamba activate r_ocelote

```
Define and prepare input file.
```
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062595.1.qualitysort_filtered_mind2.vcf
# VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf
CHR=CM062595.1       # chromosome/scaffold name
WIN=50000            # window size
STEP=50000           # step size (non-overlapping)

# make sure VCF is indexed
# bgzip "$VCF"
tabix -p vcf "$VCF.gz"
```
Compute windowed pi stats
```
# compute π in windows
vcftools --gzvcf "$VCF.gz" --chr "$CHR" \
         --window-pi $WIN --window-pi-step $STEP \
         --out ${CHR}.pi_${WIN}
```
Generate a plot of the output
```
Rscript 1a_chr29_pi.r
```
### SNP Density
Set up working directory
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference
mkdir -p snpdensity_chr29
cd snpdensity_chr29
```
Visualize the density of snps in 100kb non-overlapping windows
```
VCF="/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062595.1.qualitysort_filtered_mind2.vcf.gz"

while read -r IND; do
vcftools --gzvcf $VCF --SNPdensity 100000 --indv ${IND} --out snpden.${IND} --non-ref-ac 1 --mac 1 --max-mac 2
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
```
Merge individual outputs into one big matrix.
```
LIST="/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt"
OUT="merged.snpden.tsv"

# read individuals from list
mapfile -t INDS < "$LIST"

# collect existing snpden files
files=()
for IND in "${INDS[@]}"; do
  f="snpden.${IND}.snpden"
  if [[ -s "$f" ]]; then
    files+=("$f")
  else
    echo "Warning: missing or empty file: $f" >&2
  fi
done

if [[ ${#files[@]} -eq 0 ]]; then
  echo "No snpden files found. Exiting." >&2
  exit 1
fi

# merge with AWK, then sort by CHROM and BIN_START
awk -v OFS='\t' -v header="$(printf "%s " "${INDS[@]}")" '
BEGIN {
  n = split(header, H, " ");
  # print header line
  printf "CHROM\tBIN_START";
  for (i = 1; i <= n; i++) printf "\t%s", H[i];
  print "";
}
FNR == 1 { next }  # skip per-file header
{
  # deduce IND from filename: snpden.<IND>.snpden
  fn = FILENAME;
  sub(/.*snpden\./, "", fn);
  sub(/\.snpden.*/, "", fn);

  key = $1 OFS $2;      # CHROM \t BIN_START
  keys[key] = 1;
  chrom[key] = $1;
  bin[key]   = $2;

  # third column is SNP_COUNT
  cnt[key, fn] = $3;
}
END {
  for (k in keys) {
    printf "%s\t%s", chrom[k], bin[k];
    for (i = 1; i <= n; i++) {
      ind = H[i];
      v = cnt[k, ind];
      if (v == "") v = 0;  # fill missing with 0
      printf "\t%s", v;
    }
    print "";
  }
}
' "${files[@]}" | sort -k1,1 -k2,2n > "$OUT"

echo "Wrote: $OUT"
```
Plot the output
```
Rscript 1b_chr29_snpdensity.r
```
### Heterozygosity
Set up working directory and environment
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference

module load vcftools bcftools samtools htslib micromamba
source ~/.bashrc
micromamba activate r_ocelote

mkdir -p hetplot_chr29

cd hetplot_chr29

CHR=CM062595.1                    # chromosome/scaffold to plot (change as needed)
OUT=hets_${CHR}                   # prefix for intermediates
VCF="/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062595.1.qualitysort_filtered_mind2.vcf.gz"
# index if needed
tabix -p vcf "$VCF"
```
Filter to just keep biallelic SNPs on the chromosome of interest for clean heterozygosity calling.
```
# keep biallelic SNPs on that chromosome (clean het calling)
bcftools view -r "$CHR" -v snps -m2 -M2 -O z -o ${OUT}.snps.vcf.gz "$VCF"
tabix -p vcf ${OUT}.snps.vcf.gz
```
Convert to a genotype matrix while retaining sample names in a reference file.
```
# sample order (matches columns in the query below)
bcftools query -l ${OUT}.snps.vcf.gz > ${OUT}.samples.txt

# CHROM POS then one GT per sample (e.g., 0/0, 0/1, 1|0, 1/1, ./.)
bcftools query -f '%CHROM\t%POS[\t%GT]\n' ${OUT}.snps.vcf.gz > ${OUT}.genotypes.tsv
```

### Repeat elements
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference

mkdir -p repeat_elements
cd repeat_elements

module load micromamba

# micromamba create -n repeats -c conda-forge -c bioconda repeatmasker rmblast trf

micromamba activate repeats
micromamba update repeatmasker

echo 'CM062595.1' > ChrsOfInterest.txt
echo 'CM062599.1' >> ChrsOfInterest.txt
echo 'CM062600.1' >> ChrsOfInterest.txt

samtools faidx -r ChrsOfInterest.txt $REF > ChrsOfInterest.fa

mkdir -p vertebrata

cd vertebrata

BASE=/xdisk/mcnew/dannyjackson/sulidae
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

RepeatMasker -pa 16 -e ncbi -species "Vertebrata" -a -gff ../ChrsOfInterest.fa
```
### Final plot
Subset the GFFs to this chromosome
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/hetplot_chr29

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff
grep 'CM062595.1' $GFF | grep 'ID\=gene' > chr29.PC.genes.gff

GFF_RE=/xdisk/mcnew/dannyjackson/sulidae/analyses/repeats/vertebrata/ChrsOfInterest.fa.out.gff
grep 'CM062595.1' $GFF_RE > chr29.RE.gff

GFF_MB=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
grep 'CM062595.1' $GFF_MB | grep 'ID\=gene' > chr29.MB.genes.gff
```
I investigated the functions of various genes of interest in differentiated regions of Chromosome 29 using the scripts in ***1c_identify_genes_of_interest_chr29.md***.


Generate a combined plot of all the information across chromosome 29:
```
module load micromamba
source ~/.bashrc
micromamba activate r_ocelote
Rscript 1d_chr29_heterozygosity.r
/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/snpdensity_chr29/merged.snpden.tsv
```







sbatch 1_call29.sh