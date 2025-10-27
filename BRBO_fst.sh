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

sp=( "BFBO_PEBO", "MABO_NABO", "BRBO_RFBO" )

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
-w 1 -s 1




# run interactively
sp=( "BRBO" )

source ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh

WIN=1
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
grep 'CM' "$WIN_OUT" | grep -Ev 'CM062595|CM062600|CM062610' > "${WIN_OUT}.chrom"

# replace header (for whatever reason, it lacks a label for the fst column)
echo -e 'region\tchr\tmidPos\tNsites\tfst' > "${WIN_OUT}.chrom.txt"

cat ${WIN_OUT}.chrom >> "${WIN_OUT}.chrom.txt" 

awk 'NR==1 || $NF > 0.95' "${WIN_OUT}.chrom.txt" > "${WIN_OUT}.chrom.fixed.txt" 

# filter to fixed snps
/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/BRBOatlcar__BRBOpac_.fst_1.outlier.csv

/xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs

awk 'NR==FNR {keep[$2":"$3]; next} ($1":"$2) in keep' "${WIN_OUT}.chrom.fixed.txt"  \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs \
  > "${WIN_OUT}.chrom.fixed.mafs"

grep -v "RFBO" /xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBOabba.txt

### MAKE PCA OF SITES THAT ARE FIXED IN BROWN BOOBIES ###

#!/usr/bin/env bash
#SBATCH --job-name=BRBO_sites
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBO_sites.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO_sites.sh
# compiled on elgato

BAMLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBOabba.txt
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/BRBO_abba/analysesfiles/

# ~/programs/angsd_elgato/angsd/angsd sites index /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.mafs

~/programs/angsd_elgato/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 \
    -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.mafs \
    -bam ${BAMLIST} -anc ${REF} -doPlink 2 -out genolike_all -doGeno -1 -dopost 1  -nThreads 12


awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt | grep -v "RFBO" \
    > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO.abba.pca.txt 

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN {OFS="\t"} {print $1, $2, $2}' /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.mafs > BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.range

awk 'BEGIN{OFS="\t"} {print $1, $2, $2, "r"NR}' \
  BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.range \
  > BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.range4

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed \
    --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim \
    --extract range BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.range4 \
    --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO.abba.pca.txt --out genolike_pruned_BRBO.abba --make-bed 


module load python/3.11 plink

pcangsd -p genolike_pruned_BRBO.abba -o genolike_pruned_BRBO.abba -t 12 -e 2 --selection --admix


# plot it in R

# head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO_abba.txt
grep -v 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO_abba.txt

mkdir -p ~/R/library_elgato

R
.libPaths("~/R/library_elgato")

# install.packages("ggplot2", lib="~/R/library_elgato", repos="https://cloud.r-project.org")

library(ggplot2)
C <- as.matrix(read.table("genolike_pruned_BRBO.abba.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO_abba.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)


# --- Scores dataframe for plotting ---
scores <- data.frame(
  PC1 = e$vectors[, 1],
  PC2 = e$vectors[, 2],
  Sample  = labs$Sample,
  Species = labs$Species
)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 6, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO_pac = "#bdbdbd",
    BRBO_atlcar = "#636363",
    MABO = "#E8CA3D",
    NABO = "#B32904"
  )) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)"),
    fill = "Species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("pca_pruned_BRBO.abba.species.pdf", p, width = 6.5, height = 5.5)


p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 6, alpha = 0.9, color = "white", stroke = 0.3) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)"),
    fill = "Species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("pca_pruned_BRBO.abba.sample.pdf", p, width = 6.5, height = 5.5)

# plot admixture
tbl=read.table("genolike_pruned_BRBO.abba.admix.3.Q")
pdf(file = "admix_pruned_BRBO.abba.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()




############################################ 
# Filter MAFS to JUST sites that are fixed between brown and red footed boobies
############################################ 
#!/usr/bin/env bash
#SBATCH --job-name=BRBO_RFBO
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/BRBO_RFBO.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO_RFBO.sh
# compiled on puma

# cd /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO_RFBO

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/BRBO_abba/BRBO_RFBO_fst
sp=BRBO_RFBO
~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh \
-w 1 -s 1


# run interactively
sp=( "BRBO_RFBO" )

source ~/programs/SulidaeGenomics/param_files/${sp}_params_fst.sh

WIN=1
WIN_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${WIN}/${POP1}_${POP2}.${WIN}.fst"
grep 'CM' "$WIN_OUT" | grep -Ev 'CM062595|CM062600|CM062610' > "${WIN_OUT}.chrom"

# replace header (for whatever reason, it lacks a label for the fst column)
echo -e 'region\tchr\tmidPos\tNsites\tfst' > "${WIN_OUT}.chrom.txt"

cat ${WIN_OUT}.chrom >> "${WIN_OUT}.chrom.txt" 

awk 'NR==1 || $NF > 0.95' "${WIN_OUT}.chrom.txt" > "${WIN_OUT}.chrom.fixed.txt" 

# filter to fixed snps
awk 'NR==FNR {keep[$2":"$3]; next} ($1":"$2) in keep' "${WIN_OUT}.chrom.fixed.txt"  \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.mafs \
  > "${WIN_OUT}.chrom.fixed.mafs"





# IDENTIFY SITES FIXED BETWEEN BRBO POPS AND IN RFBO

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/BRBO_abba/analysis
BRBO1=/xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/BRBOatlcar_.mafs.gz
BRBO2=/xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/BRBOpac_.mafs.gz
BRBO12=/xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/BRBO.mafs.gz   # pairwise comparison
RFBO=/xdisk/mcnew/dannyjackson/sulidae/datafiles/safs/RFBO.mafs.gz

# keep chrom, pos, and knownEM
zcat $BRBO1 | cut -f1,2,5 > BRBO1.simple
zcat $BRBO2 | cut -f1,2,5 > BRBO2.simple
zcat $RFBO  | cut -f1,2,5 > RFBO.simple

# Ensure sorted order
sort -k1,1 -k2,2n BRBO1.simple > a
sort -k1,1 -k2,2n BRBO2.simple > b
sort -k1,1 -k2,2n RFBO.simple  > c

# Merge three mafs tables on chrom and pos
awk 'NR==FNR {k=$1"_"$2; a[k]=$3; next}
     FNR==NR {k=$1"_"$2; b[k]=$3; next}
     {k=$1"_"$2; if (k in a && k in b) print $1, $2, a[k], b[k], $3}' a b c > merged.txt

awk '($3<0.05 && $4>0.95 || $3>0.95 && $4<0.05) && $5<0.05' merged.txt > fixed_diff_sites.txt

wc -l fixed_diff_sites.txt

### MAKE PCA OF SITES THAT ARE FIXED IN BROWN BOOBIES AND BETWEEN BRBO AND RFBO ###


#!/usr/bin/env bash
#SBATCH --job-name=BRBORFBO_sites
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/BRBORFBO_sites.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBORFBO_sites.sh
# compiled on elgato

BAMLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBOabba.txt
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/BRBO_abba/analysesfiles/

# ~/programs/angsd_elgato/angsd/angsd sites index /xdisk/mcnew/dannyjackson/sulidae//analyses/fst/BRBO_RFBO/1/BRBO_RFBO.1.fst.chrom.fixed.mafs

~/programs/angsd_elgato/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 \
    -sites /xdisk/mcnew/dannyjackson/sulidae//analyses/fst/BRBO_RFBO/1/BRBO_RFBO.1.fst.chrom.fixed.mafs \
    -bam ${BAMLIST} -anc ${REF} -doPlink 2 -out genolike_BRBORFBO -doGeno -1 -dopost 1  -nThreads 12

gunzip genolike_BRBORFBO.mafs.gz

# intersect these
/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.mafs
genolike_BRBORFBO.mafs

awk 'NR==FNR { key[$1"_"$2]; next } ($1"_"$2) in key' \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBOatlcar__BRBOpac_/1/BRBOatlcar__BRBOpac_.1.fst.chrom.fixed.mafs \
  genolike_BRBORFBO.mafs > intersect.BRBOpac_RFBO.mafs


awk 'BEGIN {OFS="\t"} {print $1, $2, $2}' genolike_BRBORFBO.mafs > BRBOvsBRBOfixed.BRBOvsRFBOfixed.range

awk 'BEGIN{OFS="\t"} {print $1, $2, $2, "r"NR}' \
  BRBOvsBRBOfixed.BRBOvsRFBOfixed.range \
  > BRBOvsBRBOfixed.BRBOvsRFBOfixed.range4

plink --tped /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_all.tped --tfam /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_all.tfam --allow-extra-chr --snps-only 'just-acgt' --out genolike_all2 --make-bed 


plink --bed genolike_all2.bed \
    --bim genolike_all2.bim \
    --extract range BRBOvsBRBOfixed.BRBOvsRFBOfixed.range4 \
    --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO.abba.pca.txt --out genolike_pruned_BRBO.abba --make-bed 


module load python/3.11 plink

pcangsd -p genolike_pruned_BRBO.abba -o genolike_pruned_BRBO.abba -t 12 -e 2 --selection --admix


# plot it in R

# head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO_abba.txt
grep -v 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO_abba.txt

mkdir -p ~/R/library_elgato

R

# install.packages("ggplot2", lib="~/R/library_elgato", repos="https://cloud.r-project.org")

library(ggplot2)
C <- as.matrix(read.table("genolike_pruned_BRBO.abba.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO_abba.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[3]/sum(e$values))*100
PC2.PV.full = (e$values[4]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)


# --- Scores dataframe for plotting ---
scores <- data.frame(
  PC1 = e$vectors[, 3],
  PC2 = e$vectors[, 4],
  Sample  = labs$Sample,
  Species = labs$Species
)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 6, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO_pac = "#bdbdbd",
    BRBO_atlcar = "#636363",
    MABO = "#E8CA3D",
    NABO = "#B32904"
  )) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)"),
    fill = "Species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("pca_pruned_BRBO.abba.species.2.pdf", p, width = 6.5, height = 5.5)


p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 6, alpha = 0.9, color = "white", stroke = 0.3) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)"),
    fill = "Species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("pca_pruned_BRBO.abba.sample.2.pdf", p, width = 6.5, height = 5.5)

# plot admixture
tbl=read.table("genolike_pruned_BRBO.abba.admix.3.Q")
pdf(file = "admix_pruned_BRBO.abba.2.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()

