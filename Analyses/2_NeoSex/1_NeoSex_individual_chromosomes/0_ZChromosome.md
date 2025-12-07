# Plot statistics across individual chromosomes of interest.

## 0. Z chromosome
### Pi
Set up working directory and environment.
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference

mkdir -p piplot_chrZ
cd piplot_chrZ

module load micromamba vcftools
source ~/.bashrc
micromamba activate r_ocelote

```
Define and prepare input file.
```
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062600.1.qualitysort_filtered_mind2.vcf
# make sure VCF is indexed
bgzip "$VCF"
tabix -p vcf "$VCF.gz"
```
Compute windowed pi stats
```
CHR=CM062600.1
WIN=1000000
STEP=1000000
vcftools --gzvcf "$VCF.gz" --chr "$CHR" \
    --window-pi $WIN --window-pi-step $STEP \
    --out ${CHR}.pi_${WIN}
```
Generate a plot of the output
```
Rscript 0a_chrZ_pi.r
```
### SNP Density
Set up working directory
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference
mkdir -p snpdensity_chrZ
cd snpdensity_chrZ
```
Visualize the density of snps in 1Mb non-overlapping windows
```
VCF="/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062600.1.qualitysort_filtered_mind2.vcf.gz"

while read -r IND; do
vcftools --gzvcf $VCF --SNPdensity 1000000 --indv ${IND} --out snpden.${IND} --non-ref-ac 1 --mac 1 --max-mac 2
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
Rscript 0b_chrZ_snpdensity.r
```

### Heterozygosity
Set up working directory and environment
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference

module load vcftools bcftools samtools htslib micromamba
source ~/.bashrc
micromamba activate r_ocelote

mkdir -p hetplot_chrZ

cd hetplot_chrZ

CHR=CM062600.1
OUT=hets_${CHR}
VCF="/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062600.1.qualitysort_filtered_mind2.vcf.gz"
tabix -p vcf "$VCF"
```

Filter to just keep biallelic SNPs on the chromosome of interest for clean heterozygosity calling.
```
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
Plot it. (Can run interactively with mem=6Gb)
```
Rscript 0c_chrZ_heterozygosity.r
```
