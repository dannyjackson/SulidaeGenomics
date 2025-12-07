# Compute SNP density across chromosomes of interest

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/snpdensity/haploid/chrZ
# visualize the density of snps
VCF="/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_haploid/Z_haploid.qualitysort.vcf.gz"

# 10 kb non-overlapping windows
while read -r IND; do
vcftools --gzvcf $VCF --SNPdensity 10000000 --indv ${IND} --out snpden.${IND} --non-ref-ac 1 --mac 1 --max-mac 2
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

# merge into one big matrix
#!/usr/bin/env bash

# inputs
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
