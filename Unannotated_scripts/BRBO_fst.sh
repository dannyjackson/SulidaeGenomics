# BRBO FST


#!/usr/bin/env bash
#SBATCH --job-name=SAF_BRBOatlcar
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/SAF_BRBOatlcar.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch SAF_BRBOatlcar.sh
# compiled on puma

BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams # directory with bamfiles
BAMFILE=${BAMDIR}/${IND}.final.bam
CHR_FILE="/xdisk/mcnew/dannyjackson/sulidae/referencelists/GCA_031468815_chromconversion.txt"
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# Generate SAF Files

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs

sp=BRBOatlcar_
ls "${BAMDIR}/BRBO203.final.bam" > /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt
ls "${BAMDIR}/BRBO205.final.bam" >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt \
    -out /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/${sp} \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 \
    -anc ${REF} -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs  \
    -nThreads 1 



#!/usr/bin/env bash
#SBATCH --job-name=SAF_BRBOpac
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/SAF_BRBOpac.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch SAF_BRBOpac.sh
# compiled on puma

BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams # directory with bamfiles
BAMFILE=${BAMDIR}/${IND}.final.bam
CHR_FILE="/xdisk/mcnew/dannyjackson/sulidae/referencelists/GCA_031468815_chromconversion.txt"
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# Generate SAF Files

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs

sp=BRBOpac_
ls "${BAMDIR}/BRBO201.final.bam" > /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt
ls "${BAMDIR}/BRBO202.final.bam" >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp}bams.txt \
    -out /xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/${sp} \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 \
    -anc ${REF} -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs  \
    -nThreads 1 



cd /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO



################################################
# Run sliding window Fst
################################################
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO

chmod +x ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh 



################################################
# BRBO
################################################

#!/usr/bin/env bash
#SBATCH --job-name=BRBO
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBO.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO

sp=( "BRBO" )

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 50000 -s 50000

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 1 -s 1



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



# pull annotations
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_
module load bedtools2
#GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
CANDWIN=BRBOatlcar__BRBOpac_.fst_50000.outlier.csv
FSTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/BRBOatlcar__BRBOpac_
PREFIX=BRBOatlcar__BRBOpac_.fst_50000

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

grep 'ID=gene' ${GENEFILE} \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($9, a, ";");
      id=""; name="";
      for (i=1; i<=n; i++) {
        if (a[i] ~ /^Name=/) id=a[i];
      }
      print id;
    }' |  awk '{FS = "="} {print $2}' | sort -u > ${GENENAMES}



# SNP Analysis

awk '
NR == 1 {
    # print header with changed 5th column name
    print $1, $2, $3, $4, "fst";
    next
}
$5 > 0.99
' OFS='\t' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO/BRBO.fixedsites

awk 'NR > 1 { 
    start = $3 - 1; 
    end   = $3; 
    print $2, start, end 
}' OFS='\t' \
/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO/BRBO.fixedsites \
> /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO/BRBO.fixedsites.bed

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

bedtools intersect \
    -a /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO/BRBO.fixedsites.bed \
    -b $GFF \
    -wa -wb \
> /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO.fixedsites_in_genes.tsv


grep 'CM' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO.fixedsites_in_genes.tsv \
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
    }' |  awk '{FS = "="} {print $2}' | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/BRBO/BRBO.fixedsites_in_genes.genenames.tsv


# Really not a lot of peaks to investigate in the 50kb windows thing but one on chr 13 is interesting
"(79219,79803)(7550188,7599906)(7550000,7600000)",13,7575000,586,0.855958,3.41696522824763,3.49946599712483

echo -e "CM062579.1\t7550000\t7600000" > window.bed

bedtools intersect \
    -a window.bed \
    -b $GFF \
    -wa -wb 
	

/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff