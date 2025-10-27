# run interactively
BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
MASK_DIR=$BASE/datafiles/snpable_masks/mask
OUT=$BASE/analyses/raxml_windows_5kb_snp

# Make 5kb windows per scaffold
# samtools faidx $REF
cut -f1,2 ${REF}.fai > $OUT/genome.sizes

module load bedtools/2.31.1
module load samtools

## non-overlapping 5kb windows
bedtools makewindows -g $OUT/genome.sizes -w 5000 -s 20000 > $OUT/windows/all.5kb.bed

# Intersect windows with callable mask
## Z chromosome mask
zcat $BASE/datafiles/snpable_masks/mask/GCA_031468815*CM062600*mask.150.50.bed.gz > $MASK_DIR/GCA_031468815.1_Zchr_callable.bed

### Example: a genome-wide mask $MASK_DIR/genome_callable.bed
bedtools coverage -a $OUT/windows/all.5kb.bed -b $MASK_DIR/GCA_031468815.1_Zchr_callable.bed \
  | awk '$7>=0.8' OFS="\t" > $OUT/Zphylogeny/windows/Zchr.5kb.callable80.bed

# In each of the split 5-kbp-window alignments, taxa and sites that only contain ambiguities (N or n) or gaps (−) were removed

WINBED=$OUT/Zphylogeny/windows/Zchr.5kb.callable80.bed
### If no mask, just use:
### WINBED=$OUT/windows/all.5kb.bed




## Chromosome 4 mask
zcat $BASE/datafiles/snpable_masks/mask/GCA_031468815*CM062570*mask.150.50.bed.gz > $MASK_DIR/GCA_031468815.1_Chr4_callable.bed


### Example: a genome-wide mask $MASK_DIR/genome_callable.bed
bedtools coverage -a $OUT/windows/all.5kb.bed -b $MASK_DIR/GCA_031468815.1_Chr4_callable.bed \
  | awk '$7>=0.8' OFS="\t" > $OUT/Chr4phylogeny/windows/Chr4.5kb.callable80.bed

# In each of the split 5-kbp-window alignments, taxa and sites that only contain ambiguities (N or n) or gaps (−) were removed

WINBED=$OUT/Chr4phylogeny/windows/Chr4.5kb.callable80.bed
### If no mask, just use:
### WINBED=$OUT/windows/all.5kb.bed