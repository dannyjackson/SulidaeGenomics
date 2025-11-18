# Identify TE and repeats

# Install once (conda works well)
module load micromamba
# micromamba create -n repeats -c conda-forge -c bioconda repeatmasker rmblast trf
micromamba activate repeats
# micromamba install repeatmodeler

# Run on your assembly

BASE=/xdisk/mcnew/dannyjackson/sulidae
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna


samtools faidx -r ChrsOfInterest.txt $REF > ChrsOfInterest.fa
##############################################################################
# Build your own database (took forever and didnt work)
##############################################################################


# 1) Build the RM2 database
BuildDatabase -name MorusDB -engine ncbi $REF

# 2) De novo discovery (+ LTR structure inference)


#!/usr/bin/env bash
#SBATCH --job-name=repeats_DB
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_output/repeats_DB.%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/repeats
module load micromamba

# Install once (conda works well)
# micromamba create -n repeats -c conda-forge -c bioconda repeatmasker rmblast trf
micromamba activate repeats

RepeatModeler -database MorusDB -pa 16 -engine ncbi 

# 3) Coordinates: mask using ONLY your new library (no taxon lookups)
LIB=$(ls -1 RM_*/consensi.fa.classified | tail -n1)
RepeatMasker -pa 16 -e rmblast -a -gff -lib "$LIB" -dir RM_from_RM2 ChrsOfInterest.fa

# 4) Handy BED from GFF
awk 'BEGIN{OFS="\t"} !/^#/{
  print $1, $4-1, $5, $9, ".", $7
}' RM_from_RM2/genome.fasta.out.gff > RM2_TEcoords.bed




##############################################################################
# Use aves database
##############################################################################

export MAMBA_ROOT_PREFIX=/xdisk/mcnew/dannyjackson/.local/share/micromamba

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/repeats
module load micromamba

# Install once (conda works well)
# micromamba create -n repeats -c conda-forge -c bioconda repeatmasker rmblast trf
micromamba activate repeats
micromamba update repeatmasker
cd archelosauria
RepeatMasker -pa 16 -e ncbi -species "Archelosauria" -a -gff ../ChrsOfInterest.fa

cd vertebrata
RepeatMasker -pa 16 -e ncbi -species "Vertebrata" -a -gff ../ChrsOfInterest.fa

RepeatMasker \
  -pa 16 \                # threads
  -e rmblast \            # search engine (rmblast or nhmmer)
  -species "aves" \       # or your clade/species; else omit to use Dfam broadly
  -a -gff -dir RM_out \   # save alignments and GFF in this folder
  genome.fasta

cd  /home/u15/dannyjackson/.local/share/mamba/envs/repeats/share/RepeatMasker/Libraries/famdb
wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.2.h5.gz
wget https://www.dfam.org/releases/current/families/FamDB/dfam39_full.2.h5.gz.md5