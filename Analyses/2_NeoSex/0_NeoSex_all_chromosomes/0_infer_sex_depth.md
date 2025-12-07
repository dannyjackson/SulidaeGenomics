# Identify sex of individuals based on bamfiles post-filtering

Define variables and set up environment
```
PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
OUTDIR=${PROJDIR}/analyses/sexing_by_depth
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams/

module load micromamba 
micromamba activate r_ocelote
mkdir -p ${OUTDIR}
cd ${OUTDIR}
```
Compute depth for each chromosome of each individual
```
while read -r bird; do
samtools idxstats ${BAMDIR}/"$bird".final.bam | grep 'CM' > ${OUTDIR}/"$bird".idxstats
done <  ${PROJDIR}/referencelists/allsamplecodes.txt 
```
Create a list of all files to use in R.
```
ls "$PWD"/*idxstats > ${PROJDIR}/referencelists/idxindex.txt
```
Generate various plots in R.
```
Rscript 0a_infer_sex_depth.r
````
Notes:

MABO305 is the outlier on chr 33 with high depth. This chromosome is too noisy and should be left out of other analyses.

Chromosome 6 shows some variation across species, with all species having very high depth compared to red-footed boobies. I suppose low depth in red-footed boobies could emerge due to a structural difference in their chromosome 6? Evaluate further.

Follow up: I looked at normalized heterozygosity across the genome and across windows of Chr6 and it appears that there is something a bit strange about Chr6 in red-footed boobies but nothing interpretable. They have higher heterozygosity but not it is sex-associated or anything, and it is even across the chromosome (suggesting there isn't an inversion driving the signal).