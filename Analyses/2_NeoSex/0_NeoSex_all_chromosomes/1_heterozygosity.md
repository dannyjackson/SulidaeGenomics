# Genome-wide average heterozygosity
This script computes average heterozygosity across each chromosome by individual and then plots these stats based on sex. This allows us to infer which microchromosome in the gannet referenge genome is a neosex chromosome in boobies.

All of these scripts are quick enough that they can be run interactively.

First, set working directory and load in modules.
``` 
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/hetplot_allchr

module load vcftools samtools htslib 
# VCFtools (0.1.16)
# samtools 1.10
# htslib 1.10

module load micromamba vcftools
# 2.0.2
source ~/.bashrc
micromamba activate r_ocelote
```
This version of vcftools breaks if the header defines the vcf as a version > 4.1, so I'm recoding the header to say these versions are 4.3. This seems a bit hacky on the surface, but the vcf format is not different in a way that affects the analyses -- it's just that the vcftools version was written to break if it is passed a vcf of a more recent version than when the vcftools was written.
```
while read -r CHROM; do
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/"${CHROM}".qualitysort_filtered_mind2.vcf        # bgzipped + indexed VCF
echo $VCF
sed -i 's/VCFv4\.3/VCFv4\.1/' ${VCF}
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
```
Compute individual heterozygosity for each chromosome from our vcf using vcftools.
```
while read -r CHROM; do
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/"${CHROM}".qualitysort_filtered_mind2.vcf        # bgzipped + indexed VCF

vcftools --gzvcf "$VCF" --het --out "${CHROM}"

done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
```
Compute average heterozygosity from these output files using this equation:
Observed heterozygosity H_O = (N_Sites - O(HOM)) / N_Sites  (from the .het file)
```
while read -r CHROM; do
  infile="${CHROM}.het"
  outfile="${CHROM}.with_HO.tsv"

  echo "Processing ${infile}..."
  awk 'BEGIN {OFS="\t"} 
       NR==1 {print $0, "H_O"; next} 
       {
         H_O = ($4 - $2) / $4;  # (N_SITES - O(HOM)) / N_SITES
         print $0, H_O
       }' "$infile" > "$outfile"
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
```
Create a general matrix of all chromosomes and individuals for plotting.
```
awk 'FNR==1 && NR!=1 { next }
     {
       split(FILENAME,a,"[./]");
       CHROM=a[1]"."a[2];
       if($1=="INDV") next;
       print CHROM, IND, $0
     }' *.with_HO.tsv \
| awk 'BEGIN{OFS="\t"; print "CHROM","IND","O(HOM)","E(HOM)","N_SITES","F", "H_O"}1' \
> all_het_matrix.tsv

awk '{$1=$1; gsub(/ +/, "\t"); print}' all_het_matrix.tsv > all_het_matrix.tsv.tmp
mv all_het_matrix.tsv.tmp all_het_matrix.tsv

```
Generate plots of heterozygosity by sex by chromosome (run interactively).
```
Rscript 1a_plot_het.r
```