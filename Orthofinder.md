# Orthofinder

##################################################
# Installation
##################################################
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder

micromamba create --name orthofinder_ocelote orthofinder


##################################################
# Prepare data
##################################################
micromamba activate orthofinder_ocelote
micromamba install gffread -c bioconda

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
FASTA=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
gffread $GFF -g $FASTA -y files/bMorBas2.proteins.faa   # proteins
awk '
  # When we hit a header:
  /^>/ {
    # If this isn’t the first record, decide whether to print the previous one
    if (NR > 1 && seq !~ /\./) {
      print header
      # print sequence in wrapped lines of length 60 (optional)
      for (i=1; i<=length(seq); i+=60)
        print substr(seq, i, 60)
    }
    header = $0
    seq = ""
    next
  }

  # Sequence lines: concatenate
  {
    seq = seq $0
  }

  END {
    # Handle last record
    if (seq !~ /\./) {
      print header
      for (i=1; i<=length(seq); i+=60)
        print substr(seq, i, 60)
    }
  }
' bMorBas2.proteins.faa > bMorBas2.proteins.nodots.faa

# gffread $GFF -g $FASTA  -x files/bMorBas2.transcripts.cds # CDS nucleotide (optional)
cp /xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/protein.faa files/bPhaCar2.protein.faa
##################################################
# Run OrthoFinder
##################################################
#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=94
#SBATCH --mem=470G 
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/orthofinder

module load micromamba 

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder

micromamba run -n orthofinder_ocelote \
    orthofinder -f files -t $SLURM_CPUS_PER_TASK









##################################################
# See if GENESPACE runs into the same memory issue
##################################################

module load bedtools2 bedops minimap2


GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
FASTA=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
GFF2=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff
bedtools2


awk -v FS='\t' -v OFS='\t' '$3=="mRNA" {
  id=""; name="";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^ID=/)   id=substr(a[i],4);
    if(a[i] ~ /^Name=/) name=substr(a[i],6);
  }
  if(name=="") name=id;
  print $1, $4-1, $5, name
}' "$GFF" > bMorBas2.bed

awk -v FS='\t' -v OFS='\t' '$3=="CDS" {
  id=""; name="";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    if(a[i] ~ /^ID=/)   id=substr(a[i],8);
    if(a[i] ~ /^Name=/) name=substr(a[i],6);
  }
  if(name=="") name=id;
  print $1, $4-1, $5, id
}' "$GFF2" > bPhaCar2.bed

awk '{print $1}' bPhaCar2.protein.faa > bPhaCar2.protein.clean.faa



.libPaths("/xdisk/mcnew/dannyjackson/R/library_ocelote")

library(GENESPACE)
gpar <- init_genespace(
  wd = "/xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder/files/for_genespace", 
  path2mcscanx = "/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/r_ocelote/bin/")


gpar <- run_genespace(gsParam = gpar)



#!/usr/bin/env bash 
#SBATCH --job-name=genespace
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --ntasks=12
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/genespace
#SBATCH --mail-type=ALL

module load bedtools2 bedops minimap2 micromamba
micromamba activate r_ocelote

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder

Rscript genespace.R