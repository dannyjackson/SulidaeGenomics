BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams # directory with bamfiles
BAMFILE=${BAMDIR}/${IND}.final.bam

# Generate SAF Files

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs

~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh

mkdir /xdisk/mcnew/dannyjackson/sulidae/paramfiles
cp /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh /xdisk/mcnew/dannyjackson/sulidae/paramfiles/params_preprocessing.sh
cp /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_base.sh /xdisk/mcnew/dannyjackson/sulidae/paramfiles/params_base.sh
cp /home/u15/dannyjackson/programs/DarwinFinches/base_setup.sh /xdisk/mcnew/dannyjackson/sulidae/paramfiles/base_setup.sh

species=( "MABO" "NABO" "BFBO" "PEBO" "RFBO" "BRBO")
for sp in "${species[@]}"; do 
ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams/$sp*.final.bam > /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt
done

ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams/*.final.bam > /xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt

CHR_FILE="/xdisk/mcnew/dannyjackson/sulidae/referencelists/GCA_031468815_chromconversion.txt"
tail -n +2 SequenceReport.txt | awk 'BEGIN{OFS=","} {print $4,$7}' > $CHR_FILE

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# Create list of all SNPs across all individuals to use as input for SAF files for each species

#!/usr/bin/env bash
#SBATCH --job-name=snpID
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=20
#SBATCH --mem=470G
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/snpID.out
#SBATCH --mail-type=ALL


while read -r line; do
  
    OUTBASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/${line}_popgen

    ~/programs/angsd/angsd -GL 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 \
      -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 15 -minQ 30 -minMapQ 30 \
      -minMaf 0.05 \
      -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt \
      -out $OUTBASE \
      -P 20 \
      -r $line
  
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt


zcat /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.mafs.gz | awk '{print $1, $2, $3, $4}' > /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites.mafs

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites.mafs > /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs

~/programs/angsd/angsd sites index /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs

sp=BFBO
~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt -out /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/${sp} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 4 -minQ 30 -minMapQ 30 -anc ${REF} -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs  -nThreads 1 


~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh \
-p /xdisk/mcnew/dannyjackson/sulidae/paramfiles/params_preprocessing.sh \
-n BFBO





#!/usr/bin/env bash
#SBATCH --job-name=safs
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/safs.out
#SBATCH --mail-type=ALL

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

species=/xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt

for sp in `cat $species`; do 
~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt -out /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/${sp} -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 4 -minQ 30 -minMapQ 30 -anc ${REF} -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs  -nThreads 18
done



### Redo with 

#!/usr/bin/env bash
#SBATCH --job-name=snpID
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=20
#SBATCH --mem=470G
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/snpID.out
#SBATCH --mail-type=ALL


while read -r line; do
  
    OUTBASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/${line}_selection

    ~/programs/angsd/angsd -GL 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 \
      -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 5 -minQ 30 -minMapQ 30 \
      -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt \
      -out $OUTBASE \
      -P 20 \
      -r $line
  
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt


################################################
# Run sliding window Fst
################################################
mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/fst
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/fst

sp=( "BFBO_PEBO", "MABO_NABO", "BRBO_RFBO" )

chmod +x ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh 


################################################
# BFBO PEBO
################################################

sp=( "BFBO_PEBO" )

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 500000 -s 500000




source ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh

WIN=500000
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
grep 'CM' "$WIN_OUT" | grep -Ev 'CM062595|CM062599|CM062600|CM062610' > "${WIN_OUT}.chrom"

# replace header (for whatever reason, it lacks a label for the fst column)
echo -e 'region\tchr\tmidPos\tNsites\tfst' > "${WIN_OUT}.chrom.txt"

cat ${WIN_OUT}.chrom >> "${WIN_OUT}.chrom.txt" 
# Replace chromosome names if conversion file is provided
if [ -n "$CHR_FILE" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "${WIN_OUT}.chrom.txt"
    done < "$CHR_FILE"
fi


# z transform windowed data
Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/ztransform_windows.r \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}.chrom.txt" "${WIN}" "${POP1}_${POP2}"

Z_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${POP1}_${POP2}.fst.${WIN}.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/manhattanplot.p.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"




~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 50000 -s 50000



source ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh

WIN=50000
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
grep 'CM' "$WIN_OUT" | grep -Ev 'CM062595|CM062599|CM062600|CM062610' > "${WIN_OUT}.chrom"

# replace header (for whatever reason, it lacks a label for the fst column)
echo -e 'region\tchr\tmidPos\tNsites\tfst' > "${WIN_OUT}.chrom.txt"

cat ${WIN_OUT}.chrom >> "${WIN_OUT}.chrom.txt" 
# Replace chromosome names if conversion file is provided
if [ -n "$CHR_FILE" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "${WIN_OUT}.chrom.txt"
    done < "$CHR_FILE"
fi


# z transform windowed data
Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/ztransform_windows.r \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}.chrom.txt" "${WIN}" "${POP1}_${POP2}"

Z_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${POP1}_${POP2}.fst.${WIN}.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"


# SNP analysis

#!/bin/bash
#SBATCH --job-name=fst_bfbopebo_snps
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=94
#SBATCH --mem=100G 
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/fst_bfbopebo_snps

sp=BFBO_PEBO
FST_INDEX=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO.fst.idx
WIN_OUT=/xdisk/mcnew/dannyjackson/sulidae//analyses/fst/BFBO_PEBO/1/BFBO_PEBO.1.fst

~/programs/angsd/misc/realSFS fst stats2 "$FST_INDEX" -win 1 -step 1 -P 94 > "$WIN_OUT"

awk '
NR == 1 {
    # print header with changed 5th column name
    print $1, $2, $3, $4, "fst";
    next
}
$5 > 0.99
' OFS='\t' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/1/BFBO_PEBO.1.fst \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites

awk 'NR > 1 { 
    start = $3 - 1; 
    end   = $3; 
    print $2, start, end 
}' OFS='\t' \
/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites \
> /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites.bed

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

bedtools intersect \
    -a /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites.bed \
    -b $GFF \
    -wa -wb \
> /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO.fixedsites_in_genes.tsv

wc -l /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO.fixedsites_in_genes.tsv

grep 'CM' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO.fixedsites_in_genes.tsv \
  | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID=gene' \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($12, a, ";");
      id=""; name="";
      for (i=1; i<=n; i++) {
        if (a[i] ~ /^gene=/) id=a[i];
      }
      print id;
    }' |  awk '{FS = "="} {print $2}' | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/BFBO_PEBO/BFBO_PEBO.fixedsites_in_genes.genenames.tsv


#!/bin/bash
#SBATCH --job-name=fst_mabonabo_snps
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=94
#SBATCH --mem=100G 
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/fst_mabonabo_snps

sp=MABO_NABO
FST_INDEX=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO.fst.idx
WIN_OUT=/xdisk/mcnew/dannyjackson/sulidae//analyses/fst/MABO_NABO/1/MABO_NABO.1.fst

~/programs/angsd/misc/realSFS fst stats2 "$FST_INDEX" -win 1 -step 1 -P 94 > "$WIN_OUT"

awk '
NR == 1 {
    # print header with changed 5th column name
    print $1, $2, $3, $4, "fst";
    next
}
$5 > 0.99
' OFS='\t' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/1/MABO_NABO.1.fst \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fixedsites


awk 'NR > 1 { 
    start = $3 - 1; 
    end   = $3; 
    print $2, start, end 
}' OFS='\t' \
/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fixedsites \
> /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fixedsites.bed

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

bedtools intersect \
    -a /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fixedsites.bed \
    -b $GFF \
    -wa -wb \
> /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO.fixedsites_in_genes.tsv

wc -l /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO.fixedsites_in_genes.tsv

grep 'CM' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO.fixedsites_in_genes.tsv \
  | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID=gene' \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($12, a, ";");
      id=""; name="";
      for (i=1; i<=n; i++) {
        if (a[i] ~ /^gene=/) id=a[i];
      }
      print id;
    }' |  awk '{FS = "="} {print $2}' | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MABO_NABO/MABO_NABO.fixedsites_in_genes.genenames.tsv


# pull annotations
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
CANDWIN=BFBO_PEBO.fst_50000.outlier.csv
FSTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/BFBO_PEBO
PREFIX=BFBO_PEBO.fst_50000

awk 'BEGIN {
  FPAT = "([^,]*)|(\"[^\"]+\")"  # treat quoted fields as single units
  OFS = "\t"
}
NR==1 { print "chromo", "position"; next }
{ gsub(/"/, "", $2); gsub(/"/, "", $3); print $2, $3-25000, $3+25000 }' \
"$FSTDIR/$PREFIX.outlier.csv" | tail -n +2 > "$OUTDIR/$PREFIX.outlier.bed"

BEDFILE="$OUTDIR/$PREFIX.outlier.bed"
GENEFILE="$OUTDIR/$PREFIX.genelist.txt"
GENENAMES="$OUTDIR/$PREFIX.genenames.txt"
GENEMAPS="$OUTDIR/$PREFIX.genecoords.txt"

# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt 

python3 ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py $FSTDIR/$PREFIX.outlier.csv 50000 ${BEDFILE} ${CHROM_CONVERSION}

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > ${GENENAMES}

ALX1
CABIN1
CACNA1D
DDT
GSS # https://link.springer.com/article/10.1007/s00439-015-1559-0
DNAJC8
EYA3
LRRIQ1
MYOCD
NCAM1
PCDH7
PCDH9
PTAFR

grep 'ID\=gene' ${GENEFILE} | awk '{print $1,$4,$9}'| awk '{OFS = "\t"} {split($3, arr, ";"); print($1, $2, arr[3])}' | sed 's/Name\=//g' 

CM062568.1      41612805        PCDH9
CM062570.1      62385350        PCDH7
CM062571.1      38046086        ALX1
CM062571.1      38078899        LRRIQ1

CM062577.1      13630623        MYOCD
CM062577.1      14853603        description=glutathione
CM062577.1      14867193        DDT
CM062577.1      14871638        CABIN1

CM062580.1      16596162        CACNA1D
CM062588.1      7406571 NCAM1
CM062590.1      691686  DNAJC8
CM062590.1      715580  PTAFR
CM062590.1      724850  EYA3

# MorBas
awk '{print $1, $2, $4}' gene_window_overlaps.10kb.tsv | sort | uniq



################################################
# MABO NABO
################################################

sp=( "MABO_NABO" )


~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 500000 -s 500000




source ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh

WIN=500000
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
grep 'CM' "$WIN_OUT" | grep -Ev 'CM062595|CM062599|CM062600|CM062610' > "${WIN_OUT}.chrom"

# replace header (for whatever reason, it lacks a label for the fst column)
echo -e 'region\tchr\tmidPos\tNsites\tfst' > "${WIN_OUT}.chrom.txt"

cat ${WIN_OUT}.chrom >> "${WIN_OUT}.chrom.txt" 
# Replace chromosome names if conversion file is provided
if [ -n "$CHR_FILE" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "${WIN_OUT}.chrom.txt"
    done < "$CHR_FILE"
fi


# z transform windowed data
Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/ztransform_windows.r \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}.chrom.txt" "${WIN}" "${POP1}_${POP2}"

Z_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${POP1}_${POP2}.fst.${WIN}.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"




~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 50000 -s 10000



source ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh

WIN=50000
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
grep 'CM' "$WIN_OUT" | grep -Ev 'CM062600|CM062610' > "${WIN_OUT}.chrom"

# replace header (for whatever reason, it lacks a label for the fst column)
echo -e 'region\tchr\tmidPos\tNsites\tfst' > "${WIN_OUT}.chrom.txt"

cat ${WIN_OUT}.chrom >> "${WIN_OUT}.chrom.txt" 
# Replace chromosome names if conversion file is provided
if [ -n "$CHR_FILE" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "${WIN_OUT}.chrom.txt"
    done < "$CHR_FILE"
fi


# z transform windowed data
Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/ztransform_windows.r \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}.chrom.txt" "${WIN}" "${POP1}_${POP2}"

Z_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${POP1}_${POP2}.fst.${WIN}.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"




# pull annotations
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

CANDWIN=${sp}.fst_50000.outlier.csv
FSTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/${sp}
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/${sp}
PREFIX=${sp}.fst_50000

awk 'BEGIN {
  FPAT = "([^,]*)|(\"[^\"]+\")"  # treat quoted fields as single units
  OFS = "\t"
}
NR==1 { print "chromo", "position"; next }
{ gsub(/"/, "", $2); gsub(/"/, "", $3); print $2, $3-25000, $3+25000 }' \
"$FSTDIR/$PREFIX.outlier.csv" | tail -n +2 > "$OUTDIR/$PREFIX.outlier.bed"

BEDFILE="$OUTDIR/$PREFIX.outlier.bed"
GENEFILE="$OUTDIR/$PREFIX.genelist.txt"
GENENAMES="$OUTDIR/$PREFIX.genenames.txt"
GENEMAPS="$OUTDIR/$PREFIX.genecoords.txt"

# CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt
CHROM_CONVERSION=/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt 

python3 ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/outlier_to_bed.py $FSTDIR/$PREFIX.outlier.csv 50000 ${BEDFILE} ${CHROM_CONVERSION}

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > ${GENENAMES}

