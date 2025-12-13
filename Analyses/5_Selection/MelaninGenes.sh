Melanin genes in blue-footed, peruvian, and brown boobies


# RAxML of Candidate Introgressed Genes
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees


# RAB11A
module load micromamba bcftools raxml-ng
micromamba activate r_ocelote

mkdir RAB11A
cd RAB11A
mkdir trees

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

REGION=$(grep 'RAB11A' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')


INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/

bcftools view -r "${REGION}" "${INDIR}"/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz > RAB11A.vcf


python3 ~/programs/vcf2phylip.py3 \
    -i RAB11A.vcf \
    -o RAB11A.phy \
    -r 

python3 ~/programs/ascbias.py -p RAB11A.phy -o RAB11A.noINV.phy

module load raxml-ng

raxml-ng --all \
    --msa RAB11A.noINV.phy\
    --model GTR+G+ASC_LEWIS \
    --prefix trees/RAB11A \
    --bs-trees 200 \
    --seed 13 \
    --redo

# EDNRA

mkdir EDNRA
cd EDNRA
mkdir trees

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

REGION=$(grep 'EDNRA' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')


INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/

bcftools view -r "${REGION}" "${INDIR}"/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz > EDNRA.vcf


python3 ~/programs/vcf2phylip.py3 \
    -i EDNRA.vcf \
    -o EDNRA.phy \
    -r 

python3 ~/programs/ascbias.py -p EDNRA.phy -o EDNRA.noINV.phy


raxml-ng --all \
    --msa EDNRA.noINV.phy\
    --model GTR+G+ASC_LEWIS \
    --prefix trees/EDNRA \
    --bs-trees 200 \
    --seed 13 \
    --redo


# RAB32

mkdir RAB32
cd RAB32
mkdir trees

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

REGION=$(grep 'RAB32' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')


INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/

bcftools view -r "${REGION}" "${INDIR}"/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz > RAB32.vcf


python3 ~/programs/vcf2phylip.py3 \
    -i RAB32.vcf \
    -o RAB32.phy \
    -r 

python3 ~/programs/ascbias.py -p RAB32.phy -o RAB32.noINV.phy


raxml-ng --all \
    --msa RAB32.noINV.phy\
    --model GTR+G+ASC_LEWIS \
    --prefix trees/RAB32 \
    --bs-trees 200 \
    --seed 13 \
    --redo


# functional prediction

# subset VCFs to just BFBO and PEBO
bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 EDNRA/EDNRA.vcf > EDNRA/BFBO_PEBO.EDNRA.vcf
bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 RAB32/RAB32.vcf > RAB32/BFBO_PEBO.RAB32.vcf
bcftools view -s BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 RAB11A/RAB11A.vcf > RAB11A/BFBO_PEBO.RAB11A.vcf

# Create directory for this new genome
cd ~/programs/snpEff/
mkdir bMorBas
cd bMorBas

# Get annotation files
cp /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff genes.gff
bgzip genes.gff

# Get the genome reference sequence file
cd  ~/programs/snpEff/data/genomes
cp /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna bMorBas.fa
bgzip bMorBas.fa

echo 'bMorBas.genome : bMorBas' >> ~/programs/snpEff/snpEff.config
echo 'data_dir  : /home/u15/dannyjackson/programs/snpEff/data/' >> snpEff.config

/home/u15/dannyjackson/programs/snpEff/./data/bMorBas/genes.gff
cd ~/programs/snpEff

java -jar ~/programs/snpEff/snpEff.jar build -gff2 -v bMorBas


~/programs/jdk-23.0.2/bin/java -Xmx8g -jar snpEff.jar \
    build \
    -gff3 \
    -c snpEff.config \
    -v bMorBas \
    -noCheckCds -noCheckProtein 
    
# Annotating VCFs


~/programs/jdk-23.0.2/bin/java -Xmx8g -jar ~/programs/snpEff/snpEff.jar \
    -c ~/programs/snpEff/snpEff.config bMorBas \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/GPR148/GPR148.vcf \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/GPR148/GPR148.ann.vcf


~/programs/jdk-23.0.2/bin/java -Xmx8g -jar ~/programs/snpEff/snpEff.jar \
    -c ~/programs/snpEff/snpEff.config bMorBas \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/RAB11A/BFBO_PEBO.RAB11A.vcf \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/RAB11A/BFBO_PEBO.RAB11A.ann.vcf


REGION=$(grep 'RAB11A' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')

FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST" > "${FST%.txt}.RAB11A.txt"
# (14135,14135)(3254642,3254642)(3254642,3254643) CM062581.1      3254642 2       0.990177

VCF_ANN=/xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/RAB11A/BFBO_PEBO.RAB11A.ann.vcf

grep '3254642' $VCF_ANN # G-A
samtools faidx "$REF" "$REGION" > "${REGION}.fa"
samtools faidx "$REF" CM062581.1:3254641-3254643 #CGC ... in BFBO is now CAC
samtools faidx "$REF" CM062581.1:3254635-3254650 #TCCAAACGCCTCTGTA ... in BFBO is now TCCAAACACCTCTGTA... a CPG site

grep 'RAB11A' $GFF # it's in the negative direciton but it doesn't matter. CGC vs GCG is still a CpG site



~/programs/jdk-23.0.2/bin/java -Xmx8g -jar ~/programs/snpEff/snpEff.jar \
    -c ~/programs/snpEff/snpEff.config bMorBas \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/EDNRA/BFBO_PEBO.EDNRA.vcf \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/EDNRA/BFBO_PEBO.EDNRA.ann.vcf

REGION=$(grep 'EDNRA' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')

FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST" > "${FST%.txt}.EDNRA.txt"

# region  chr     midPos  Nsites  fst
# (58094,58094)(17154819,17154819)(17154819,17154820)     CM062570.1      17154819        2       0.990700
# (58169,58169)(17175158,17175158)(17175158,17175159)     CM062570.1      17175158        2       0.990632

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna 
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
VCF_ANN=/xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/EDNRA/BFBO_PEBO.EDNRA.ann.vcf

grep '17154819' $VCF_ANN # G-A
samtools faidx "$REF" "$REGION" > "${REGION}.fa"
samtools faidx "$REF" CM062570.1:17154818-17154820 #TGC ... in PEBO is now TAC
samtools faidx "$REF" CM062570.1:17154810-17154825 #GGCTTGGCTGCGCTAA ... in PEBO is now GGCTTGGCTACGCTAA... evaluate for potential 

grep '17175158' $VCF_ANN # C-T
samtools faidx "$REF" CM062570.1:17175157-17175159 #TCG ... in PEBO is now TTG... now a CpG

grep 'EDNRA' $GFF


~/programs/jdk-23.0.2/bin/java -Xmx8g -jar ~/programs/snpEff/snpEff.jar \
    -c ~/programs/snpEff/snpEff.config bMorBas \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/RAB32/BFBO_PEBO.RAB32.vcf \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/RAB32/BFBO_PEBO.RAB32.ann.vcf


REGION=$(grep 'RAB32' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')

FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST" > "${FST%.txt}.RAB32.txt"
# (176099,176099)(43145479,43145479)(43145479,43145480)   CM062569.1      43145479        2       0.990222

VCF_ANN=/xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/RAB32/BFBO_PEBO.RAB32.ann.vcf


bcftools view -r CM062569.1:43145478-43145480 /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus_snps_multiallelic.vcf.gz
bcftools index /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus_snps_multiallelic.vcf.gz 

grep '43145479' $VCF_ANN # C-G
samtools faidx "$REF" "$REGION" > "${REGION}.fa"
samtools faidx "$REF" CM062569.1:43145478-43145480 # ACA is mutated to AGA
samtools faidx "$REF" CM062569.1:43145470-43145490 # TCCAAGCTACATGTCATTCCC TCCAAGCTAGATGTCATTCCC
samtools faidx "$REF" CM062569.1:43145460-43145500 
# EDNRA
echo -e 'bMorBas2_egapxtmp_002163\t17154819,17175158' > genes_snps.tsv

# RAB11A
echo -e 'bMorBas2_egapxtmp_009595\t3254642' >> genes_snps.tsv

# RAB32
echo -e 'bMorBas2_egapxtmp_000187\t43145479' >> genes_snps.tsv

### Visualize genes with SNPs (one panel per gene, stacked via par)
library(Gviz)
library(rtracklayer)

## =========================
## Inputs
## =========================
gff_file      <- "/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff"
gene_snp_file <- "genes_snps.tsv"  # tab-separated: geneID \t pos1,pos2,pos3,...
genome_name   <- "bMorBas2"
out_pdf       <- "genes_with_SNPs.stack.pdf"

## =========================
## Read GFF and gene/SNP table
## =========================
gff  <- import(gff_file)
meta <- mcols(gff)

id_vec <- as.character(meta$ID)

if ("Parent" %in% colnames(meta)) {
  parent_vec <- vapply(meta$Parent, function(x) paste(x, collapse = ","), character(1))
} else {
  parent_vec <- rep("", length(id_vec))
}

features_to_keep <- c("gene", "mRNA", "exon", "CDS")

# gene/SNP table: col1 = gene ID, col2 = comma-separated SNP positions
gene_snp_df <- read.delim(
  gene_snp_file,
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(gene_snp_df) <- c("gene_id", "snp_string")

# Parse SNP positions into list-column
gene_snp_df$snp_pos <- lapply(
  gene_snp_df$snp_string,
  function(x) {
    x <- trimws(unlist(strsplit(x, ",")))
    as.numeric(x[nzchar(x)])
  }
)

## =========================
## Build per-gene objects to plot
## =========================
options(ucscChromosomeNames = FALSE)

gene_feats_list <- list()
gene_chr        <- character(0)
gene_from       <- numeric(0)
gene_to         <- numeric(0)
gene_snps       <- list()
gene_ids_kept   <- character(0)

for (i in seq_len(nrow(gene_snp_df))) {
  gene_id <- gene_snp_df$gene_id[i]
  snp_pos <- gene_snp_df$snp_pos[[i]]

  # mask for this gene: keep basic features and anything whose ID or Parent contains the gene_id
  mask <- meta$type %in% features_to_keep &
          ( grepl(gene_id, id_vec, fixed = TRUE) |
            grepl(gene_id, parent_vec, fixed = TRUE) )

  gene_feats <- gff[mask]

  if (length(gene_feats) == 0) {
    warning("No features found for gene_id: ", gene_id)
    next
  }

  chr_vec <- as.character(unique(seqnames(gene_feats)))
  if (length(chr_vec) != 1) {
    warning("Multiple chromosomes for gene_id ", gene_id, ": ",
            paste(chr_vec, collapse = ","))
    next
  }
  chr <- chr_vec[1]

  rng  <- range(gene_feats)
  from <- start(rng)
  to   <- end(rng)

  # optional: add a little padding around the gene
  pad  <- max(0, round((to - from) * 0.05))
  from <- from - pad
  to   <- to + pad

  gene_feats_list[[length(gene_feats_list) + 1]] <- gene_feats
  gene_chr        <- c(gene_chr, chr)
  gene_from       <- c(gene_from, from)
  gene_to         <- c(gene_to, to)
  gene_snps[[length(gene_snps) + 1]] <- snp_pos
  gene_ids_kept   <- c(gene_ids_kept, gene_id)
}

n_plots <- length(gene_feats_list)

if (n_plots == 0) {
  stop("No valid gene tracks were created; nothing to plot.")
}

## =========================
## Plot: one panel per gene, stacked with par(mfrow)
## =========================
pdf(out_pdf, width = 10, height = max(3, 1.5 * n_plots))
par(mfrow = c(n_plots, 1), mar = c(4, 4, 3, 2))

for (j in seq_len(n_plots)) {
  gene_feats <- gene_feats_list[[j]]
  chr        <- gene_chr[j]
  from       <- gene_from[j]
  to         <- gene_to[j]
  snp_pos    <- gene_snps[[j]]
  gene_id    <- gene_ids_kept[j]

  axisTrack <- GenomeAxisTrack()

  # Gene track: same style as before, but no "unknown" labels
  geneTrack <- GeneRegionTrack(
    gene_feats,
    genome     = genome_name,
    chromosome = chr,
    name       = gene_id,   # label on the left
    showId     = FALSE,     # do not print IDs on exons
    featureAnnotation = NULL
  )

  # SNP track under this gene
  snpTrack <- AnnotationTrack(
    start      = snp_pos,
    end        = snp_pos,
    chromosome = chr,
    genome     = genome_name,
    name       = "SNPs",
    shape      = "box"
  )

  plotTracks(
    list(axisTrack, geneTrack, snpTrack),
    from = from,
    to   = to,
    chromosome = chr,
    main = gene_id
  )
}

dev.off()




EDNRA
RAB11A
EDNRA
RAB32
PIKFYVE
HPS5
EDNRA
RAB32
PIKFYVE
HPS5
ASIP
RAB32
EDNRA
HPS5
EDNRA
HPS5
CPOX
TYR
COX10
ALAS1
CPOX
TYR
COX10
ALAS1


EDNRA was identified in comparisons of masked vs. Nazca, blue-footed vs. Peruvian, and across populations of masked boobies, 

~/programs/jdk-23.0.2/bin/java -Xmx8g -jar ~/programs/snpEff/snpEff.jar \
    -c ~/programs/snpEff/snpEff.config bMorBas \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/EDNRA/EDNRA.vcf \
    > /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/EDNRA/EDNRA.ann.vcf

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

REGION=$(grep 'EDNRA' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')

FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BFBO_PEBO/BFBO_PEBO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST" > "${FST%.txt}.RAB11A.txt"
# region  chr     midPos  Nsites  fst
# (58094,58094)(17154819,17154819)(17154819,17154820)     CM062570.1      17154819        2       0.990700
# (58169,58169)(17175158,17175158)(17175158,17175159)     CM062570.1      17175158        2       0.990632

FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO_NABO/MABO_NABO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST"
# (27481,27481)(17179448,17179448)(17179448,17179449)     CM062570.1      17179448        2       0.998548
VCF_ANN=/xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees/EDNRA/EDNRA.ann.vcf
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna 

grep '17179448' $VCF_ANN # T-C
samtools faidx "$REF" CM062570.1:17179447-17179449 #GTA ... in PEBO is now GCA... unknown

# HPS5 PIKFYVE 
# RAB32

REGION=$(grep 'RAB32' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')


FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/MABO/MABO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST" 
# (455065,455065)(43148400,43148400)(43148400,43148401)   CM062569.1      43148400        2       0.994088
VCF_ANN=RAB32.ann.vcf
grep '43148400' $VCF_ANN 

FST=/xdisk/mcnew/dannyjackson/sulidae/analyses/fst/BRBO/BRBO.fixedsites
awk -v R="$REGION" -v OFS="\t" '
BEGIN {
    split(R, a, "[:\\-]")
    chr = a[1]
    start = a[2]
    end = a[3]
}
NR==1 { print; next }   # keep header
($2 == chr && $3 >= start && $3 <= end)
' "$FST" 
# 804,438804)(43158309,43158309)(43158309,43158310)   CM062569.1      43158309        2       0.995212 | C - T
# (438835,438835)(43161505,43161505)(43161505,43161506)   CM062569.1      43161505        2       0.996971 | G - C
grep '43158309' $VCF_ANN 
grep '43161505' $VCF_ANN 

samtools faidx "$REF" CM062570.1:43158308-43158310 #AGC ... in PEBO is now ATC... unknown
samtools faidx "$REF" CM062570.1:43161504-43161506 #TGT ... in PEBO is now TCT... unknown
