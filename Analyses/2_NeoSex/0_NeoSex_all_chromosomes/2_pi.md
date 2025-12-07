# Nucleotide diversity across the genome

Set working directory and environment.
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/

mkdir -p piplot_allchr
cd piplot_allchr

module load vcftools samtools htslib

```
Compute pi across each chromosome for each individual.

```
sbatch 2a_pi.sh
```
Compute the average pi for each chromosome for each individual.
```

# Directory containing the *.pi.sites.pi files
PI_DIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_heterozygosity_filtered/hetplot_allchr/pifiles
OUT=mean_pi_matrix.tsv

# Write header
{
  printf "CHROM"
  while read -r IND; do printf "\t%s" "$IND"; done < ${INDFILE}
  printf "\n"
} > "$OUT"

# For each chromosome, compute mean π per individual
while read -r CHROM; do
  printf "%s" "$CHROM" >> "$OUT"
  while read -r IND; do
    file="${PI_DIR}/${CHROM}.${IND}.pi.sites.pi"
    if [[ -f "$file" ]]; then
      mean=$(awk 'NR>1 && $3 ~ /^[0-9.eE+-]+$/ {sum+=$3; n++} END{if(n>0) printf("%.10g", sum/n); else printf("NA")}' "$file")
    else
      mean="NA"
    fi
    printf "\t%s" "$mean" >> "$OUT"
  done < ${INDFILE}
  printf "\n" >> "$OUT"
done < ${CONTIGS}

echo "Wrote mean π matrix to: $OUT"
```
Plot the output.
```
Rscript 2b_pi.r
```