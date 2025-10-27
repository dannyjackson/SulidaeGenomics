# Annotating

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff
MaskedBoobyFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/masked_booby/ncbi_dataset/data/GCA_013389905.1/GCA_013389905.1_ASM1338990v1_genomic.fna
MaskedBoobyGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/masked_booby/ncbi_dataset/data/GCA_013389905.1/genomic.gff

# masked booby b10k
~/programs/datasets download genome accession GCA_013389905.1 --include gff3,rna,cds,protein,genome,seq-report

# great cormorant
~/programs/datasets download genome accession GCF_963921805.1 --include gff3,rna,cds,protein,genome,seq-report


conda install -c bioconda liftoff


#!/usr/bin/env bash
#SBATCH --job-name=LiftOff_PhaCar
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/LiftOff_PhaCar.%A_%a.out
#SBATCH --mail-type=ALL

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"
# micromamba create -n LiftGFF -c conda-forge -c bioconda liftoff 
micromamba activate LiftGFF 

liftoff -g $PhaCarGFF -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff -p 24 $REF $PhaCarFasta





#!/usr/bin/env bash
#SBATCH --job-name=LiftOff_MABO
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/LiftOff_MABO.%A_%a.out
#SBATCH --mail-type=ALL

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
MaskedBoobyFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/masked_booby/ncbi_dataset/data/GCA_013389905.1/GCA_013389905.1_ASM1338990v1_genomic.fna
MaskedBoobyGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/masked_booby/ncbi_dataset/data/GCA_013389905.1/genomic.gff


source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"
# micromamba create -n LiftGFF -c conda-forge -c bioconda liftoff 
micromamba activate LiftGFF 

liftoff -g $MaskedBoobyGFF -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.MABO.hap2_genomic_lifted.gff -p 24 $REF $MaskedBoobyFasta

# make background gene list 
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

awk -F'\t' 'BEGIN{OFS="\t"}
  $1 ~ /^CM/ && $1 !~ /^CM062600\.1$/ && $1 !~ /^CM062610\.1$/ && $9 ~ /ID=gene/ {
    split($9, a, ";");
    sub(/^ID=gene-/, "", a[1]);
    print a[1]
  }' GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff | sort -u > GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.backgroundgenes.autosomes.txt



#!/usr/bin/env bash
#SBATCH --job-name=LiftOn_PhaCar
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/LiftOn_PhaCar.%A_%a.out
#SBATCH --mail-type=ALL

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"
# micromamba create -n LiftGFF -c conda-forge -c bioconda liftoff 
micromamba activate LiftGFF 

source ~/.bashrc

# pip3 install lifton
# micromamba install miniprot 

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

lifton -g $PhaCarGFF \
  -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff \
  -copies -t 24 \
  $REF $PhaCarFasta



REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"
# micromamba create -n LiftGFF -c conda-forge -c bioconda liftoff 
micromamba activate LiftGFF 

liftoff -g $PhaCarGFF -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff -p 24 $REF $PhaCarFasta




# Are we confident about these?
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
FASTA=/xdisk/mcnew/dannyjackson/sulidae/analyses/working/dct.fasta
module load blast/2.15.0

# make blast db
makeblastdb -in "$REF" -dbtype nucl -parse_seqids -out "${REF%.fna}"

# high-identity search
blastn -task megablast \
  -query "$FASTA" \
  -db "${REF%.fna}" \
  -evalue 1e-20 \
  -max_target_seqs 50 \
  -max_hsps 100 \
  -outfmt '6 qseqid sseqid sstart send sstrand pident length qlen evalue bitscore' \
  > dct_vs_ref.tsv


# perform synteny analysis
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny

module load micromamba

# 1) Create two separate envs (orthofinder and genespace can't share one, it seems)
# confirm you’re in bash; if blank or not bash, start one
echo $BASH_VERSION || exec bash

# set micromamba root
export MAMBA_ROOT_PREFIX="$HOME/.micromamba"

# ✅ correct hook: no nested $(which …), no extra braces
eval "$(micromamba shell hook -s bash)"


# 1.1) Tools env (OrthoFinder/DIAMOND/MCScanX)
micromamba create -n gs-tools -c bioconda -c conda-forge \
  orthofinder=2.5.5 diamond=2.0.15 mcscanx -y

# 1.2) Genespace
micromamba create -n gs-r -c conda-forge \
  r-base=4.4 r-ragg r-remotes r-biocmanager \
  pkg-config freetype libpng harfbuzz fribidi cairo -y

# 2) Run GENESPACE with tools on PATH
# activate R env
# point to where micromamba is installed (adjust path to your binary if needed)
export MAMBA_ROOT_PREFIX="$HOME/.micromamba"
eval "$(micromamba shell hook -s bash)"

micromamba activate gs-r
module load R
# add tools bin to PATH for this session
export PATH="$HOME/.micromamba/envs/gs-tools/bin:$PATH"

# sanity check
which orthofinder
which MCScanX
which diamond
R -q -e 'library(ragg); cat("ragg OK\n")'



# MCScanX (binary in your PATH)
micromamba install -c bioconda mcscanx -y
micromamba install conda-forge::r-devtools

# After having it all installed, run via:
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny

module load micromamba
module load R
# 1) Create two separate envs (orthofinder and genespace can't share one, it seems)
# confirm you’re in bash; if blank or not bash, start one
echo $BASH_VERSION || exec bash

# set micromamba root
export MAMBA_ROOT_PREFIX="$HOME/.micromamba"

# ✅ correct hook: no nested $(which …), no extra braces
eval "$(micromamba shell hook -s bash)"


micromamba activate gs-r

# add tools bin to PATH for this session
export PATH="$HOME/.micromamba/envs/gs-tools/bin:$PATH"


REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
REFGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

mkdir genomeRepo
cd genomeRepo
cp $REF MorBas/MorBas.fna
cp $REFGFF MorBas/MorBas.gff
cp $PhaCarFasta PhaCar/PhaCar.fna
cp $PhaCarGFF PhaCar/PhaCar.gff


# Analyze in R

R 

# Install genespace

# install biostrings and rtracklayer
install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c("Biostrings","rtracklayer"), ask=FALSE, update=FALSE)

# install genespace
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")
remotes::install_github("jtlovell/GENESPACE", upgrade="never")

# install dependencies
BiocManager::install(c(
    "Matrix", "Biostrings", "rtracklayer",     # core Bioconductor libs
    "ape", "ggplot2", "dplyr", "data.table",   # common CRAN deps
    "jsonlite", "Rcpp", "cowplot", "tidyr",
    "reshape2", "scales", "ggrepel", "viridis"
), ask = FALSE, update = FALSE)

# peptide and cds file generation
gffread -g MorBas/MorBas.fna -y MorBas/MorBas.pep.fa -x MorBas/MorBas.cds.fa \
        -S -F MorBas/MorBas.gff

gffread -g PhaCar/PhaCar.fna -y PhaCar/PhaCar.pep.fa -x PhaCar/PhaCar.cds.fa \
        -S -F PhaCar/PhaCar.gff

library(GENESPACE)

genomeRepo <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/genomeRepo"
wd         <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/gs_workspace"   # will hold derived files
dir.create(wd, showWarnings = FALSE)


parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo, 
  genomeDirs = c("MorBas", "PhaCar"),
  genomeIDs = c("MorBas", "PhaCar"),
  gffString = "gff",
  faString = "pep.fa",
  presets = "ncbi", 
  genespaceWd = wd)



parsed <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs    = c("MorBas","PhaCar"),   # subfolders you created
  genomeIDs     = c("MorBas","PhaCar"),       # short IDs to appear on plots
  presets       = "ncbi",                 # works well for standard GFF3
  genespaceWd   = wd,
  pepSubstr     = "pep.fa",               # ensure it picks your peptide FASTAs
  gffSubstr     = "gff",                 # use your provided GFFs
  minPepLen     = 50
)




















# Parse annotations manually
# Extract CDS and translate to proteins from GFF + genome
# (gffread chooses transcript models; we’ll collapse to one per gene below)
micromamba install -n gs-tools -c bioconda gffread seqkit -y

gffread MorBas/MorBas.gff -g MorBas/MorBas.fa  -x MorBas/cds.fa -y MorBas/proteins.fa
gffread PhaCar/PhaCar.gff -g PhaCar/PhaCar.fa -x PhaCar/cds.fa -y PhaCar/proteins.fa

# Make a BED of gene coordinates (one line per gene)
# Pull gene features; edit 'gene'→'mRNA' if your IDs live on mRNA rows instead.
awk -F'\t' '$3=="gene"{ 
  split($9,a,";"); id="";
  for(i in a){ if(a[i]~/^ID=/){sub(/^ID=/,"",a[i]); id=a[i]} }
  if(id!=""){ printf "%s\t%d\t%d\t%s\t0\t%s\n",$1,$4-1,$5,id,$7 }
}' MorBas/MorBas.gff > NOGA.genes.bed

# OPTIONAL: collapse proteins to one per gene (longest) and rename headers to gene IDs
# Create transcript→gene map from GFF:
# mRNA lines usually carry ID=<transcript> and Parent=<gene>
awk -F'\t' '$3=="mRNA"{
  split($9,a,";"); id=""; parent="";
  for(i in a){
    if(a[i]~/^ID=/){id=substr(a[i],4)}
    else if(a[i]~/^Parent=/){parent=substr(a[i],8)}
  }
  if(id!="" && parent!=""){
    sub(/\.[0-9]+$/,"",id);      # strip .1/.2 version suffix from tx
    sub(/\.[0-9]+$/,"",parent);  # strip from gene too (if present)
    print id "\t" parent
  }
}' MorBas/MorBas.gff > MorBas.txnv2genenv.tsv

# (seqkit handy to rank by length)
# Full header lengths: (full_tx\tlen)
seqkit fx2tab -n -l MorBas/proteins.fa > MorBas.proteome.len.tsv

# Map full header -> no-version transcript ID
awk '{
  tx_full=$1; tx_nv=$1; sub(/\.[0-9]+$/,"",tx_nv);
  print tx_full "\t" tx_nv
}' MorBas.proteome.len.tsv > MorBas.full2nov.tsv

# Attach tx_nv to the length table: (full_tx\tlen\ttx_nv)
awk 'BEGIN{FS=OFS="\t"} NR==FNR{nov[$1]=$2; next} {print $1,$2,nov[$1]}' \
  MorBas.full2nov.tsv MorBas.proteome.len.tsv > MorBas.len_with_nov.tsv

# Pick longest isoform per gene and rewrite protein headers to the gene ID

# Use the GFF-derived map (tx_nv -> gene_nv) to pick the longest per gene
awk -F'\t' 'NR==FNR{gene[$1]=$2; next}{
  full=$1; len=$2; txnv=$3; g=gene[txnv];
  if(g!=""){
    if(len>best[g]){best[g]=len; pick[g]=full}
  }
} END{ for(g in pick) print pick[g] }' \
  MorBas.txnv2genenv.tsv MorBas.len_with_nov.tsv > MorBas.keep.fulltx

# Produce (full_tx \t gene_nv) for all transcripts
awk -F'\t' 'NR==FNR{gene[$1]=$2; next} {print $1"\t"gene[$3]}' \
  MorBas.txnv2genenv.tsv MorBas.len_with_nov.tsv > MorBas.full2gene.tsv

# Keep only the chosen (longest) transcripts
grep -F -f MorBas.keep.fulltx MorBas.full2gene.tsv > MorBas.keep.full2gene.tsv

# Extract the sequences of the chosen full transcript IDs
seqkit grep -n -f MorBas.keep.fulltx MorBas/proteins.fa \
| awk 'NR==FNR{m[$1]=$2; next}
       /^>/{h=substr($0,2); print ">" m[h]; next}
       {print}' MorBas.keep.full2gene.tsv - > NOGA.fa

# MorBas genes BED
awk -F'\t' '$3=="gene"{
  split($9,a,";"); id="";
  for(i in a){ if(a[i]~/^ID=/){id=substr(a[i],4)} }
  if(id!=""){ printf "%s\t%d\t%d\t%s\t0\t%s\n",$1,$4-1,$5,id,$7 }
}' MorBas/MorBas.gff > NOGA.bed

# Repeat for GRCO
awk -F'\t' '$3=="gene"{ 
  split($9,a,";"); id="";
  for(i in a){ if(a[i]~(/^ID=/)){sub(/^ID=/,"",a[i]); id=a[i]} }
  if(id!=""){ printf "%s\t%d\t%d\t%s\t0\t%s\n",$1,$4-1,$5,id,$7 }
}' PhaCar/PhaCar.gff > GRCO.genes.bed


# OPTIONAL: collapse proteins to one per gene (longest) and rename headers to gene IDs
# Create transcript→gene map from GFF:
# mRNA lines usually carry ID=<transcript> and Parent=<gene>
awk -F'\t' '$3=="mRNA"{
  split($9,a,";"); id=""; parent="";
  for(i in a){
    if(a[i]~/^ID=/){id=substr(a[i],4)}
    else if(a[i]~/^Parent=/){parent=substr(a[i],8)}
  }
  if(id!="" && parent!=""){
    sub(/\.[0-9]+$/,"",id);      # strip .1/.2 version suffix from tx
    sub(/\.[0-9]+$/,"",parent);  # strip from gene too (if present)
    print id "\t" parent
  }
}' PhaCar/PhaCar.gff > PhaCar.txnv2genenv.tsv

# (seqkit handy to rank by length)
# Full header lengths: (full_tx\tlen)
seqkit fx2tab -n -l PhaCar/proteins.fa > PhaCar.proteome.len.tsv

# Map full header -> no-version transcript ID
awk '{
  tx_full=$1; tx_nv=$1; sub(/\.[0-9]+$/,"",tx_nv);
  print tx_full "\t" tx_nv
}' PhaCar.proteome.len.tsv > PhaCar.full2nov.tsv

# Attach tx_nv to the length table: (full_tx\tlen\ttx_nv)
awk 'BEGIN{FS=OFS="\t"} NR==FNR{nov[$1]=$2; next} {print $1,$2,nov[$1]}' \
  PhaCar.full2nov.tsv PhaCar.proteome.len.tsv > PhaCar.len_with_nov.tsv

# Pick longest isoform per gene and rewrite protein headers to the gene ID

# Use the GFF-derived map (tx_nv -> gene_nv) to pick the longest per gene
awk -F'\t' 'NR==FNR{gene[$1]=$2; next}{
  full=$1; len=$2; txnv=$3; g=gene[txnv];
  if(g!=""){
    if(len>best[g]){best[g]=len; pick[g]=full}
  }
} END{ for(g in pick) print pick[g] }' \
  PhaCar.txnv2genenv.tsv PhaCar.len_with_nov.tsv > PhaCar.keep.fulltx

# Produce (full_tx \t gene_nv) for all transcripts
awk -F'\t' 'NR==FNR{gene[$1]=$2; next} {print $1"\t"gene[$3]}' \
  PhaCar.txnv2genenv.tsv PhaCar.len_with_nov.tsv > PhaCar.full2gene.tsv

# Keep only the chosen (longest) transcripts
grep -F -f PhaCar.keep.fulltx PhaCar.full2gene.tsv > PhaCar.keep.full2gene.tsv

# Extract the sequences of the chosen full transcript IDs
seqkit grep -n -f PhaCar.keep.fulltx PhaCar/proteins.fa \
| awk 'NR==FNR{m[$1]=$2; next}
       /^>/{h=substr($0,2); print ">" m[h]; next}
       {print}' PhaCar.keep.full2gene.tsv - > GRCO.fa


# PhaCar genes BED
awk -F'\t' '$3=="gene"{
  split($9,a,";"); id="";
  for(i in a){ if(a[i]~/^ID=/){id=substr(a[i],4)} }
  if(id!=""){ printf "%s\t%d\t%d\t%s\t0\t%s\n",$1,$4-1,$5,id,$7 }
}' PhaCar/PhaCar.gff > GRCO.bed


library(GENESPACE)

genomeRepo <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/genomeRepo"
wd         <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/gs_workspace"
dir.create(wd, showWarnings = FALSE)

gpar <- init_genespace(wd = wd, path2mcscanx = Sys.which("MCScanX"))
gpar <- run_genespace(gsParam = gpar)

plot_riparian(gsParam = gpar,
              genomeIDs = c("NOGA","GRCO"),
              refGenome = "NOGA",
              useOrder  = TRUE,
              pdfFile   = file.path(wd, "riparian.NOGA_vs_GRCO.pdf"))















# Nucmer attempt
eval "$(micromamba shell hook -s bash)"

micromamba create -n synteny -c bioconda -c conda-forge \
  mummer syri plotsr seqkit python=3.8 -y

micromamba activate synteny

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/genomeRepo/MorBas/MorBas.fna .
cp /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/genomeRepo/PhaCar/PhaCar.fna .
# Filter out non-contigs
seqkit grep -f /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt MorBas.fna > MorBas.filtered.fna
while read -r line; do 
echo ${line}
seqkit grep -r -p ${line} ../MorBas.fna > MorBas.${line}.fna
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt

grep 'NC_' /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/genomeRepo/PhaCar/PhaCar.fna | awk '{print $1}' | sed 's/>//g' > /xdisk/mcnew/dannyjackson/sulidae/referencelists/PhaCar.CONTIGS.txt
seqkit grep -f /xdisk/mcnew/dannyjackson/sulidae/referencelists/PhaCar.CONTIGS.txt PhaCar.fna > PhaCar.filtered.fna

# Align genomes with MUMmer (nucmer)
nucmer --maxmatch -l 100 -c 500 inputfiles/MorBas.CM062598.1.fna PhaCar.filtered.fna -p outputfiles/MorBas_CM062598.1_vs_PhaCar



#!/usr/bin/env bash
#SBATCH --job-name=nucmer_chroms
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=400G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/nucmer_chroms.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-35 nucmer_chroms.sh 

set -euo pipefail

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/

INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/inputfiles
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/outputfiles_slurm
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

module load micromamba
eval "$(micromamba shell hook -s bash)"
micromamba activate synteny2

nucmer -l 100 -c 500  ${INDIR}/MorBas.${CHROM}.fna ${INDIR}/PhaCar.filtered.fna -p ${OUTDIR}/MorBas_${CHROM}_vs_PhaCar

#!/usr/bin/env bash
#SBATCH --job-name=nucmer
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=470G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_output/nucmer.out
#SBATCH --mail-type=ALL

INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/nucmer/
OUTDIR=${INDIR}/wgs_output
nucmer --maxmatch -c 100 -b 500 -l 50 ${INDIR}/MorBas.filtered.fna ${INDIR}/PhaCar.filtered.fna -p ${OUTDIR}/MorBas_vs_PhaCar_wgs
mummer merge -p merged outputfiles_slurm/*_vs_PhaCar.delta


for f in outputfiles_slurm/*_vs_PhaCar.delta; do
  show-coords -THrd "$f"
done > merged.filtered.coords

# Ref & query IDs present in the coords
awk '{r[$10]=1} END{for (k in r) print k}' merged.filtered.coords | sort > ref.ids
awk '{q[$11]=1} END{for (k in q) print k}' merged.filtered.coords | sort > qry.ids

# Quick sanity check
echo "ref IDs:  $(wc -l < ref.ids)"
echo "qry IDs:  $(wc -l < qry.ids)"

# Ensure both FASTAs are indexed
samtools faidx MorBas.fna
samtools faidx PhaCar.fna

# Build the subset FASTAs in the **same order as the coords’ ID lists**
xargs -a ref.ids  samtools faidx MorBas.fna > MorBas.syri.fna
xargs -a qry.ids  samtools faidx PhaCar.fna > PhaCar.syri.fna

syri -c merged.filtered.coords \
     -d <(cat outputfiles_slurm/*_vs_PhaCar.delta) \
     -r MorBas.syri.fna \
     -q PhaCar.syri.fna \
     -F T -k \
     --prefix syri.merged

syri -c merged.filtered.coords \
     -d <(cat outputfiles_slurm/*_vs_PhaCar.delta) \
     -r MorBas.fna \
     -q PhaCar.fna \
     -F T -k \
     --prefix syri.merged


head -n 5 outputfiles_slurm/MorBas_CM062567.1_vs_PhaCar.delta | tail -n +2 > merged.delta
for f in outputfiles_slurm/*.delta; do
    tail -n +6 "$f" >> merged.delta
done
delta-filter -m -i 90 -l 100 merged.delta > merged.filtered.delta

# Filter the alignments
while read -r CHROM; do
delta-filter -1 outputfiles_slurm/MorBas_${CHROM}_vs_PhaCar.delta > filteredfiles/MorBas_${CHROM}_vs_PhaCar.delta
show-coords -THrd filteredfiles/MorBas_${CHROM}_vs_PhaCar.delta > filteredfiles/MorBas_${CHROM}_vs_PhaCar.coords

done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt

# make a table of alignments
echo -e "Ref_Chrom\tQuery_Contigs\tRef_Contigs" > alignments_table.tsv

while read -r CHROM; do
  ref_list=$(awk '{print $10}' MorBas_${CHROM}_vs_PhaCar.coords | sort -u | paste -sd "," -)
  qry_list=$(awk '{print $11}' MorBas_${CHROM}_vs_PhaCar.coords | sort -u | paste -sd "," -)
  echo -e "${CHROM}\t${qry_list}\t${ref_list}" >> alignments_table.tsv
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt


# cat filteredfiles/MorBas_*_vs_PhaCar.delta > MorBas_vs_PhaCar.delta
# cat filteredfiles/MorBas_*_vs_PhaCar.coords > MorBas_vs_PhaCar.coords

syri -c filteredfiles/MorBas_CM062567.1_vs_PhaCar.coords \
     -d filteredfiles/MorBas_CM062567.1_vs_PhaCar.delta \
     -r MorBas.fna \
     -q PhaCar.fna \
     -F T -k \
     --prefix syri.CM062567.1 --no-chrmatch

syri -c filteredfiles/MorBas_${CHROM}_vs_PhaCar.coords \
     -d filteredfiles/MorBas_${CHROM}_vs_PhaCar.delta \
     -r MorBas.fna \
     -q PhaCar.fna \
     -F T -k \
     --prefix syri.${CHROM}


# run SyRI to call synteny and rearrangements

#  Generate sequence length tables (for plotsr)
seqkit fx2tab -n -l MorBas.fna > MorBas.sizes
seqkit fx2tab -n -l PhaCar.fna > PhaCar.sizes

# Plot with plotsr
plotsr \
  --assembly1 MorBas.fna --asm1name MorBas --asm1size MorBas.sizes \
  --assembly2 PhaCar.fna --asm2name PhaCar --asm2size PhaCar.sizes \
  --syri syri.out \
  --out MorBas_vs_PhaCar_synteny.png

# Inspect outputs
MorBas_vs_PhaCar.1delta        # filtered alignment
syri.out, syri.vcf, *.bed      # synteny/rearrangement tables
MorBas_vs_PhaCar_synteny.png   # final ribbon plot












#######################################################################
#######################################################################
# module load minimap2

# minimap2 -x asm10 -t 32 MorBas.fna PhaCar.fna > MorBas_vs_PhaCar.paf
# use -x asm5 if nearly identical; asm10 is safer for modest divergence
# minimap2 -x asm20 MorBas.fna PhaCar.fna --no-long-join -r 200 | cut -f 1-12 > MorBas_vs_PhaCar.sam

# 0) env (conda/mamba)
