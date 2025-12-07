# Scripts to call and variants
First, call variants:
```
sbatch 0_call.sh
```
Then, remove scaffolds, rename samples, and 

1. filter on quality with bcftools (QUAL>30)
2. filter on depth (--min-meanDP 4 --max-meanDP 20) and remove indels with vcftools
3. filter on call rate (--geno 0.02), minimum individual (--mind 0.2) and minor allele frequency (--maf 0.01), again keeping only SNPs.

```
sbatch 1_filter.sh
```
Finally, split the vcf into individual chromosomes.
```
sbatch --array=1-35 2_split_vcf_chrom.sh 
```