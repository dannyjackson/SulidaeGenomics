# PCA

# B1_PCA.sh with trans, subset

#!/usr/bin/env bash
#SBATCH --job-name=genolike_all
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/genolike_all%j.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
BAMLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs -bam ${BAMLIST} -anc ${REF} -doPlink 2 -out genolike_all -doGeno -1 -dopost 1  -nThreads 8


#!/usr/bin/env bash
#SBATCH --job-name=genolike_filtering
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --output=slurm_output/genolike_filtering%j.out
#SBATCH --mail-type=ALL

module load plink

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all

plink --tped /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_all.tped --tfam /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_all.tfam --allow-extra-chr --snps-only 'just-acgt' --indep-pairwise 50kb 1 0.5 --out genolike_filtered --bad-ld

plink --tped /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_all.tped --tfam /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_all.tfam --allow-extra-chr --snps-only 'just-acgt' --extract genolike_filtered.prune.in --out genolike_pruned --make-bed 



# Submitted batch job 3637791

#!/usr/bin/env bash
#SBATCH --job-name=pcangsd
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/pcangsd%j.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all
module load python
module load bcftools

# for pop gen
pcangsd -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/pruned -t 12 -e 2 --selection --admix

# save the following as: /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt

sample  species
BFBO501 BFBO
BFBO502 BFBO
BFBO503 BFBO
BFBO504 BFBO
BFBO505 BFBO
BRBO201 BRBO
BRBO202 BRBO
BRBO203 BRBO
BRBO205 BRBO
MABO302 MABO
MABO304 MABO
MABO305 MABO
MABO306 MABO
NABO402 NABO
NABO403 NABO
NABO404 NABO
NABO405 NABO
NABO406 NABO
PEBO601 PEBO
PEBO603 PEBO
PEBO604 PEBO
PEBO605 PEBO
PEBO606 PEBO
RFBO101 RFBO
RFBO102 RFBO
RFBO103 RFBO
RFBO104 RFBO
RFBO105 RFBO
RFBO106 RFBO


# /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt
# plot it in R
cd  /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all

R

C <- as.matrix(read.table("pruned.cov")) # Reads estimated covariance matrix
# D <- as.matrix(read.table("output.selection")) # Reads PC based selection statistics
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt", header = TRUE)
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

# --- Plot with ggplot2 (fill-based legend/colors) ---
library(ggplot2)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 6, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO = "#000000",
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

ggsave("pca_pruned_species.ggplot.pdf", p, width = 6.5, height = 5.5)



# plot admixture
tbl=read.table("pruned.admix.3.Q")
pdf(file = "admix_pruned.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()




# plot PCA by each species
# BFBO
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/BFBO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO.pca.txt 

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO.pca.txt --out genolike_pruned_BFBO --make-bed 

pcangsd -p genolike_pruned_BFBO -o genolike_pruned_BFBO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/BFBO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BFBO.txt
grep 'BFBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BFBO.txt

C <- as.matrix(read.table("genolike_pruned_BFBO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BFBO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
library(ggplot2)
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
    BRBO = "#000000",
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

ggsave("pca_pruned_BFBO.ggplot.pdf", p, width = 6.5, height = 5.5)


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

ggsave("pca_pruned_BFBO.sample.pdf", p, width = 6.5, height = 5.5)


# plot admixture
tbl=read.table("genolike_pruned_BFBO.admix.3.Q")
pdf(file = "admix_pruned_BFBO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()



# PEBO
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/PEBO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO.pca.txt 

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO.pca.txt --out genolike_pruned_PEBO --make-bed 

pcangsd -p genolike_pruned_PEBO -o genolike_pruned_PEBO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/PEBO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO.txt
grep 'PEBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO.txt

C <- as.matrix(read.table("genolike_pruned_PEBO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO.txt", header = TRUE)
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
    BRBO = "#000000",
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

ggsave("pca_pruned_PEBO.species.pdf", p, width = 6.5, height = 5.5)


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

ggsave("pca_pruned_PEBO.sample.pdf", p, width = 6.5, height = 5.5)



# plot admixture
tbl=read.table("genolike_pruned_PEBO.admix.3.Q")
pdf(file = "admix_pruned_PEBO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()


# MABO
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/MABO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO.pca.txt 

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO.pca.txt --out genolike_pruned_MABO --make-bed 

pcangsd -p genolike_pruned_MABO -o genolike_pruned_MABO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/MABO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO.txt
grep 'MABO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO.txt

library(ggplot2)
C <- as.matrix(read.table("genolike_pruned_MABO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO.txt", header = TRUE)
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
    BRBO = "#000000",
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

ggsave("pca_pruned_MABO.species.pdf", p, width = 6.5, height = 5.5)


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

ggsave("pca_pruned_MABO.sample.pdf", p, width = 6.5, height = 5.5)

# plot admixture
tbl=read.table("genolike_pruned_MABO.admix.3.Q")
pdf(file = "admix_pruned_MABO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()




# NABO
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/NABO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO.pca.txt 

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO.pca.txt --out genolike_pruned_NABO --make-bed 

pcangsd -p genolike_pruned_NABO -o genolike_pruned_NABO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/NABO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.NABO.txt
grep 'NABO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.NABO.txt

C <- as.matrix(read.table("genolike_pruned_NABO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.NABO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)

library(ggplot2)

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
    BRBO = "#000000",
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

ggsave("pca_pruned_NABO.species.pdf", p, width = 6.5, height = 5.5)


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

ggsave("pca_pruned_NABO.sample.pdf", p, width = 6.5, height = 5.5)

# plot admixture
tbl=read.table("genolike_pruned_NABO.admix.3.Q")
pdf(file = "admix_pruned_NABO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()





# BRBO
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/BRBO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO.pca.txt 

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO.pca.txt --out genolike_pruned_BRBO --make-bed 

pcangsd -p genolike_pruned_BRBO -o genolike_pruned_BRBO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/BRBO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO.txt
grep 'BRBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO.txt

C <- as.matrix(read.table("genolike_pruned_BRBO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.BRBO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)


library(ggplot2)

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
    BRBO = "#000000",
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

ggsave("pca_pruned_BRBO.species.pdf", p, width = 6.5, height = 5.5)


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

ggsave("pca_pruned_BRBO.sample.pdf", p, width = 6.5, height = 5.5)


# plot admixture
tbl=read.table("genolike_pruned_BRBO.admix.3.Q")
pdf(file = "admix_pruned_BRBO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()


# RFBO

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/RFBO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO.pca.txt 

plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO.pca.txt --out genolike_pruned_RFBO --make-bed 

pcangsd -p genolike_pruned_RFBO -o genolike_pruned_RFBO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/RFBO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO.txt
grep 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO.txt

C <- as.matrix(read.table("genolike_pruned_RFBO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)


library(ggplot2)

# --- Scores dataframe for plotting ---
scores <- data.frame(
  PC1 = e$vectors[, 1],
  PC2 = e$vectors[, 2],
  Sample  = labs$Sample,
  Species = labs$Species
)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO = "#000000",
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

ggsave("pca_pruned_RFBO.species.pdf", p, width = 6.5, height = 5.5)


p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
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

ggsave("pca_pruned_RFBO.sample.pdf", p, width = 6.5, height = 5.5)



# plot admixture
tbl=read.table("genolike_pruned_RFBO.admix.3.Q")
pdf(file = "admix_pruned_RFBO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()



### Plot species pairs

## RFBO vs BRBO

mkdir -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/RFBO_BRBO

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/RFBO_BRBO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_BRBO.pca.txt 
awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO_samplecodes.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_BRBO.pca.txt 


plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_BRBO.pca.txt --out genolike_pruned_RFBO_BRBO --make-bed 

pcangsd -p genolike_pruned_RFBO_BRBO -o genolike_pruned_RFBO_BRBO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/RFBO_BRBO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO_BRBO.txt
grep 'BRBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO_BRBO.txt
grep 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO_BRBO.txt

C <- as.matrix(read.table("genolike_pruned_RFBO_BRBO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.RFBO_BRBO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)



library(ggplot2)

# --- Scores dataframe for plotting ---
scores <- data.frame(
  PC1 = e$vectors[, 1],
  PC2 = e$vectors[, 2],
  Sample  = labs$Sample,
  Species = labs$Species
)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO = "#000000",
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
    legend.position = "none"
  )

ggsave("pca_pruned_RFBO_BRBO.species.pdf", p, width = 6.5, height = 5.5)


p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
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

ggsave("pca_pruned_RFBO_BRBO.sample.pdf", p, width = 6.5, height = 5.5)



# plot admixture
tbl=read.table("genolike_pruned_RFBO_BRBO.admix.3.Q")
pdf(file = "admix_pruned_RFBO_BRBO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()


# MABO vs NABO

mkdir -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/MABO_NABO

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/MABO_NABO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_NABO.pca.txt 
awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO_samplecodes.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_NABO.pca.txt 


plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_NABO.pca.txt --out genolike_pruned_MABO_NABO --make-bed 

pcangsd -p genolike_pruned_MABO_NABO -o genolike_pruned_MABO_NABO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/MABO_NABO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO_NABO.txt
grep 'MABO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO_NABO.txt
grep 'NABO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO_NABO.txt

C <- as.matrix(read.table("genolike_pruned_MABO_NABO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.MABO_NABO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)



library(ggplot2)

# --- Scores dataframe for plotting ---
scores <- data.frame(
  PC1 = e$vectors[, 1],
  PC2 = e$vectors[, 2],
  Sample  = labs$Sample,
  Species = labs$Species
)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO = "#000000",
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
    legend.position = "none"
  )

ggsave("pca_pruned_MABO_NABO.species.pdf", p, width = 6.5, height = 5.5)


p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
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

ggsave("pca_pruned_MABO_NABO.sample.pdf", p, width = 6.5, height = 5.5)



# plot admixture
tbl=read.table("genolike_pruned_MABO_NABO.admix.3.Q")
pdf(file = "admix_pruned_MABO_NABO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()


### PEBO vs BFBO

mkdir -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/PEBO_BFBO

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/PEBO_BFBO

tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_BFBO.pca.txt 
awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO_samplecodes.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_BFBO.pca.txt 


plink --bed /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_BFBO.pca.txt --out genolike_pruned_PEBO_BFBO --make-bed 

pcangsd -p genolike_pruned_PEBO_BFBO -o genolike_pruned_PEBO_BFBO -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods/PEBO_BFBO

head -n 1 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO_BFBO.txt
grep 'BFBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO_BFBO.txt
grep 'PEBO' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO_BFBO.txt

C <- as.matrix(read.table("genolike_pruned_PEBO_BFBO.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.PEBO_BFBO.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Sample <- factor(labs$sample)



library(ggplot2)

# --- Scores dataframe for plotting ---
scores <- data.frame(
  PC1 = e$vectors[, 1],
  PC2 = e$vectors[, 2],
  Sample  = labs$Sample,
  Species = labs$Species
)

p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = c(
    BFBO = "#4FE98C",
    PEBO = "#2894FF",
    RFBO = "#F3447D",
    BRBO = "#000000",
    MABO = "#E8CA3D",
    NABO = "#B32904"
  )) +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)"),
    fill = "Species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave("pca_pruned_PEBO_BFBO.species.pdf", p, width = 5.5, height = 5)


p <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 12, alpha = 0.9, color = "white", stroke = 0.3) +
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

ggsave("pca_pruned_PEBO_BFBO.sample.pdf", p, width = 5.5, height = 5.5)



# plot admixture
tbl=read.table("genolike_pruned_PEBO_BFBO.admix.3.Q")
pdf(file = "admix_pruned_PEBO_BFBO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()




# Plot all PCAs
# run_all_pca_admix.sh

BASE="/xdisk/mcnew/dannyjackson/sulidae/datafiles/geno_likelihoods"
R_SCRIPT="pca_admix_plot.multisp.R"

groups=( MABO_NABO PEBO_BFBO RFBO_BRBO)

for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  Rscript "$R_SCRIPT" "$g" "$BASE"
done


R_SCRIPT="pca_admix_plot.nonumbers.R"
groups=(all)

for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  Rscript "$R_SCRIPT" "$g" "$BASE"
done

groups=(BFBO BRBO MABO NABO PEBO RFBO)
R_SCRIPT="pca_admix_plot.onesp.R"


for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  Rscript "$R_SCRIPT" "$g" "$BASE"
done