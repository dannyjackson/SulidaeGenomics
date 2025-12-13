cd /xdisk/mcnew/dannyjackson/sulidae/analyses/manhattanplots

Plot the genome wide statistics of FST and DFM
```
Rscript CombinedManhattan.oneaxis.Chr2Highlight.r
```
## Identify genes in outlier regions identified by DFM scans. 
First, masked vs Nazca


```
# pull annotations
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

CANDWIN=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo/analyses/f_dM/MABO_AtlCar_MABO_IndoPacific_NABO/5000_200/MABO_AtlCar_MABO_IndoPacific_NABO.f_dM.5000_200.0.999.outliers.tsv
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/${sp}
PREFIX=MABO_AtlCar_MABO_IndoPacific_NABO.f_dM.5000_200.0.999.outliers

awk '{print $2, $5, $6 }' \
${CANDWIN} | tail -n +2 > "$OUTDIR/$PREFIX.outlier.bed"

tr ' ' '\t' < "$OUTDIR/$PREFIX.outlier.bed" > "$OUTDIR/$PREFIX.outlier.bed.tmp" && mv "$OUTDIR/$PREFIX.outlier.bed.tmp" "$OUTDIR/$PREFIX.outlier.bed"
BEDFILE="$OUTDIR/$PREFIX.outlier.bed"
GENEFILE="$OUTDIR/$PREFIX.genelist.txt"
GENENAMES="$OUTDIR/$PREFIX.genenames.txt"
GENEMAPS="$OUTDIR/$PREFIX.genecoords.txt"

# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt 

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > ${GENENAMES}
```
Next, blue-footed vs Peruvian
```

# pull annotations
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

CANDWIN=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/bfbo_pebo/analyses/f_dM/BFBO_GofCA_BFBO_southern_PEBO/5000_200/BFBO_GofCA_BFBO_southern_PEBO.f_dM.5000_200.0.999.outliers.tsv
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/${sp}
PREFIX=BFBO_GofCA_BFBO_southern_PEBO.f_dM.5000_200.0.999

awk '{print $2, $5, $6 }' \
${CANDWIN} | tail -n +2 > "$OUTDIR/$PREFIX.outlier.bed"

tr ' ' '\t' < "$OUTDIR/$PREFIX.outlier.bed" > "$OUTDIR/$PREFIX.outlier.bed.tmp" && mv "$OUTDIR/$PREFIX.outlier.bed.tmp" "$OUTDIR/$PREFIX.outlier.bed"
BEDFILE="$OUTDIR/$PREFIX.outlier.bed"
GENEFILE="$OUTDIR/$PREFIX.genelist.txt"
GENENAMES="$OUTDIR/$PREFIX.genenames.txt"
GENEMAPS="$OUTDIR/$PREFIX.genecoords.txt"

# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt 

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > ${GENENAMES}

```
# Finally, brown boobies and 4 taxon clade
```

# pull annotations
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

CANDWIN=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/brbo_4taxa/analyses/f_dM/BRBO_Pacific_BRBO_AtlCar_4taxa/5000_200/BRBO_Pacific_BRBO_AtlCar_4taxa.f_dM.5000_200.0.999.outliers.tsv
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/${sp}
PREFIX=BRBO_Pacific_BRBO_AtlCar_4taxa.f_dM.5000_200.0.999

awk '{print $2, $5, $6 }' \
${CANDWIN} | tail -n +2 > "$OUTDIR/$PREFIX.outlier.bed"

tr ' ' '\t' < "$OUTDIR/$PREFIX.outlier.bed" > "$OUTDIR/$PREFIX.outlier.bed.tmp" && mv "$OUTDIR/$PREFIX.outlier.bed.tmp" "$OUTDIR/$PREFIX.outlier.bed"
BEDFILE="$OUTDIR/$PREFIX.outlier.bed"
GENEFILE="$OUTDIR/$PREFIX.genelist.txt"
GENENAMES="$OUTDIR/$PREFIX.genenames.txt"
GENEMAPS="$OUTDIR/$PREFIX.genecoords.txt"

# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt 

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > ${GENENAMES}
```











Identify genes in Chr2 outlier region of MABO NABO DFM scan
```
# pull annotations
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

CANDWIN=MABO_NABO_DFM_chrCM062568.1_top0.05_outliers.tsv
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/${sp}
PREFIX=MABO_NABO_DFM_chrCM062568.1_top0.05_outliers

awk '{print $1, $2, $3 }' \
"MABO_NABO_DFM_chrCM062568.1_top0.05_outliers.tsv" | tail -n +2 > "$OUTDIR/$PREFIX.outlier.bed"

tr ' ' '\t' < "$OUTDIR/$PREFIX.outlier.bed" > "$OUTDIR/$PREFIX.outlier.bed.tmp" && mv "$OUTDIR/$PREFIX.outlier.bed.tmp" "$OUTDIR/$PREFIX.outlier.bed"
BEDFILE="$OUTDIR/$PREFIX.outlier.bed"
GENEFILE="$OUTDIR/$PREFIX.genelist.txt"
GENENAMES="$OUTDIR/$PREFIX.genenames.txt"
GENEMAPS="$OUTDIR/$PREFIX.genecoords.txt"

# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt 

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > ${GENENAMES}
```