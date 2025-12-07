cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes


# chicken, Gallus gallus, GalGal
~/programs/datasets download genome accession GCF_000002315.5 --include gff3,rna,cds,protein,genome,seq-report
# /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/ncbi_dataset/data/GCF_000002315.5/GCF_000002315.5_GRCg6a_genomic.fna

# crested ibis, Nipponia nippon, NipNip
~/programs/datasets download genome accession GCA_035839065.1 --include gff3,rna,cds,protein,genome,seq-report
# /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/NipNip/ncbi_dataset/data/GCA_035839065.1/GCA_035839065.1_ASM3583906v1_genomic.fna

# zebra finch, Taeniopygia guttata, TaeGut
~/programs/datasets download genome accession GCA_048771995.1 --include gff3,rna,cds,protein,genome,seq-report
# /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/ncbi_dataset/data/GCA_048771995.1/GCA_048771995.1_bTaeGut7.mat_genomic.fna
# Following indications of sex linkage and unexpected sex-biased gene expression in warblers (Sylvioidea; Passeriformes), we have conducted an extensive marker analysis targeting 31 orthologues of loci on zebra finch chromosome 4a in five species, representative of independent branches of Passerida.


# Sooty myzomela, Myzomela tristrami, SooMyz
~/programs/datasets download genome accession GCA_051294455.1 --include gff3,rna,cds,protein,genome,seq-report

# Myzomela cardinalis, cardinal Myzomela, CarMyz
~/programs/datasets download genome accession GCA_051294465.1 --include gff3,rna,cds,protein,genome,seq-report


# Adelie penguins also have a neo sex Z1Z2 system but don't have a chromosome level genome.
# They do however have seq data... maybe align to king penguin and ID to that? cool comparative paper later on
# king penguin, Aptenodytes patagonicus, AptPat
# ~/programs/datasets download genome accession GCA_048771995.1 --include gff3,rna,cds,protein,genome,seq-report

GCA_965638725.1

module load micromamba
module load minimap2
# module load R


micromamba activate r_elgato


#######################################################################
# installation
#######################################################################

R 

# .libPaths("/xdisk/mcnew/dannyjackson/R/library_elgato")
.libPaths("/xdisk/mcnew/dannyjackson/R/library_ocelote")
.libPaths("/xdisk/mcnew/dannyjackson/R/library_puma")


if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/DEEPSPACE")

devtools::load_all("/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/DEEPSPACE")

q()


### Install GENESPACE
devtools::install_github("jtlovell/GENESPACE")


#######################################################################
# Analysis -- Five genomes
#######################################################################

# 0) env (conda/mamba)
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/MorBas.filtered.fna fastafiles/MorBas.fa
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/PhaCar.filtered.fna fastafiles/PhaCar.fa
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/GalGal/ncbi_dataset/data/GCF_000002315.5/GCF_000002315.5_GRCg6a_genomic.fna fastafiles/GalGal
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/NipNip/ncbi_dataset/data/GCA_035839065.1/GCA_035839065.1_ASM3583906v1_genomic.fna fastafiles/NipNip.fa
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/TaeGut/ncbi_dataset/data/GCA_048771995.1/GCA_048771995.1_bTaeGut7.mat_genomic.fna fastafiles/TaeGut.fa

cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/SooMyz/ncbi_dataset/data/GCA_051294455.1/GCA_051294455.1_MyzoTris1.0.pri_genomic.fna fastafiles/SooMyz.fa

cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/CarMyz/ncbi_dataset/data/GCA_051294465.1/GCA_051294465.1_MyzoCard1.0.pri_genomic.fna fastafiles/CarMyz.fa

# clean fastas by removing scaffolds
micromamba activate genomics_elgato

cd fastafiles

seqkit grep -v -r -p "^NW_" GalGal.fa > filtered/GalGal_noNW.fa
seqkit grep -v -r -p "^JAS" NipNip.fa > filtered/NipNip_noNW.fa
seqkit grep -v -r -p "^JBM" SooMyz.fa > filtered/SooMyz_noNW.fa
seqkit grep -v -r -p "^JBM" CarMyz.fa > filtered/CarMyz_noNW.fa



# cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/sevengenomes

module load micromamba
module load minimap2

# micromamba activate r_elgato
micromamba activate r_puma
R

# .libPaths("/xdisk/mcnew/dannyjackson/R/library_elgato")
# .libPaths("/xdisk/mcnew/dannyjackson/R/library_ocelote")
# .libPaths("/xdisk/mcnew/dannyjackson/R/library_puma")

devtools::load_all("/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/DEEPSPACE")

library(DEEPSPACE)

wd <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/sevengenomes"
tmpDir <- file.path(wd, "DEEPSPACEtmpwd")

# 1) Pre-create tmpDir so deepspace will NOT auto-delete it
dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)

fastaFiles <- c(
  MorBas = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/MorBas.fa",
  PhaCar = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/PhaCar.fa",
  NipNip = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/NipNip_noNW.fa",
  GalGal = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/GalGal_noNW.fa",
  TaeGut = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/TaeGut.fa",
  SooMyz = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/SooMyz_noNW.fa",
  CarMyz = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/fastafiles/filtered/CarMyz_noNW.fa"
)

fiveG <- clean_windows(
  faFiles = fastaFiles,
  genomeIDs = c("MorBas", "PhaCar", "NipNip", "GalGal", "TaeGut", "SooMyz", "CarMyz"),
  wd = wd,
  
  preset = "dist fast",
  
  stripChrname = ".*chromosome\\s+|,.*",
  minChrLen = 2e5, 
  
  nCores = 1,
  MCScanX_hCall = Sys.which("MCScanX_h"),
  minimap2call = Sys.which("minimap2"))

"GalGal_windows__vs__MorBas.syn.paf.gz",
"NipNip_windows__vs__MorBas.syn.paf.gz",
"PhaCar_windows__vs__MorBas.syn.paf.gz",
"TaeGut_windows__vs__MorBas.syn.paf.gz",
"NipNip_windows__vs__GalGal.syn.paf.gz",
"NipNip_windows__vs__PhaCar.syn.paf.gz",
"TaeGut_windows__vs__GalGal.syn.paf.gz"

NipNip_windows__vs__MorBas.syn.paf.gz

pdf("riparian_paf.wide.pdf", width = 16, height = 6) 
riparian_paf(
  pafFiles = c(
    "GalGal_windows__vs__MorBas.syn.paf.gz",
    "NipNip_windows__vs__MorBas.syn.paf.gz",
    "PhaCar_windows__vs__MorBas.syn.paf.gz",
    "TaeGut_windows__vs__MorBas.syn.paf.gz",
    "NipNip_windows__vs__GalGal.syn.paf.gz",
    "NipNip_windows__vs__PhaCar.syn.paf.gz",
    "TaeGut_windows__vs__GalGal.syn.paf.gz"
  ),
  refGenome = "MorBas",
  genomeIDs = c("MorBas", "PhaCar", "NipNip", "GalGal", "TaeGut"),
  orderyBySynteny = FALSE,
  chrOutlineColor = "#f0f0f0ff",
  chrOutlineWidth = 1
)
dev.off()





#!/usr/bin/env bash
#SBATCH --job-name=fivegenomes_sbatch
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/fivegenomes_sbatch.out
#SBATCH --mail-type=ALL
# sbatch fivegenomes.sh

module load micromamba
module load minimap2
source ~/.bashrc
micromamba activate r_elgato

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes_sbatch
Rscript fivegenomes.R

#######################################################################
# Analysis -- three genomes
#######################################################################

minimap2 -x asm20 -t 32 /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/MorBas.filtered.fna /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/PhaCar.filtered.fna > MorBas_vs_PhaCar.paf
# use -x asm5 if nearly identical; asm10 is safer for modest divergence
# minimap2 -x asm20 MorBas.fna PhaCar.fna --no-long-join -r 200 | cut -f 1-12 > MorBas_vs_PhaCar.sam

# 0) env (conda/mamba)
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/MorBas.filtered.fna fastafiles/MorBas.fa
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/PhaCar.filtered.fna fastafiles/PhaCar.fa

R
library(DEEPSPACE)

wd <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/MorBas_vs_PhaCar_all"
tmpDir <- file.path(wd, "DEEPSPACEtmpwd")

# 1) Pre-create tmpDir so deepspace will NOT auto-delete it
dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)

fastaFiles <- c(
  MorBas = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fastafiles/MorBas.fa",
  PhaCar = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fastafiles/PhaCar.fa"
)

MbVPc <- clean_windows(
  faFiles = fastaFiles,
  genomeIDs = c("MorBas", "PhaCar"),
  wd = wd,
  
  preset = "dist fast",
  
  stripChrname = ".*chromosome\\s+|,.*",
  minChrLen = 2e5, 
  
  nCores = 1,
  MCScanX_hCall = Sys.which("MCScanX_h"),
  minimap2call = Sys.which("minimap2"))


#######################################################################
# Analysis -- NipNip vs TaeGut
#######################################################################
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/NipNip_vs_TaeGut

R
library(DEEPSPACE)

wd <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/NipNip_vs_TaeGut"
tmpDir <- file.path(wd, "DEEPSPACEtmpwd")

# 1) Pre-create tmpDir so deepspace will NOT auto-delete it
dir.create(tmpDir, recursive = TRUE, showWarnings = FALSE)

fastaFiles <- c(
  NipNip = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fastafiles/PhaCar.fa",
  TaeGut = "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fastafiles/MorBas.fa"
)

NipNip_vs_TaeGut <- clean_windows(
  faFiles = fastaFiles,
  genomeIDs = c("NipNip", "TaeGut"),
  wd = wd,
  
  preset = "dist fast",
  
  stripChrname = ".*chromosome\\s+|,.*",
  minChrLen = 2e5, 
  
  nCores = 1,
  MCScanX_hCall = Sys.which("MCScanX_h"),
  minimap2call = Sys.which("minimap2"))

