raxml_windowed

# run interactively
BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
CONTIG_FILE=${BASE}/referencelists/AUTOSOMES.txt
MASK_DIR=$BASE/datafiles/snpable_masks/mask
OUT=$BASE/analyses/raxml_windows_5kb_snp

# Make 5kb windows per scaffold
# samtools faidx $REF
cut -f1,2 ${REF}.fai > $OUT/genome.sizes

module load bedtools2
module load samtools

## non-overlapping 5kb windows
bedtools makewindows -g $OUT/genome.sizes -w 5000 -s 20000 > $OUT/windows/all.5kb.bed

# Intersect windows with callable mask
zcat $BASE/datafiles/snpable_masks/mask/GCA_031468815*CM*mask.150.50.bed.gz > $MASK_DIR/GCA_031468815.1_callable.bed

### Example: a genome-wide mask $MASK_DIR/genome_callable.bed
bedtools coverage -a $OUT/windows/all.5kb.bed -b $MASK_DIR/GCA_031468815.1_callable.bed \
  | awk '$7>=0.8' OFS="\t" > $OUT/windows/all.5kb.callable80.bed

# In each of the split 5-kbp-window alignments, taxa and sites that only contain ambiguities (N or n) or gaps (−) were removed

WINBED=$OUT/windows/all.5kb.callable80.bed
### If no mask, just use:
### WINBED=$OUT/windows/all.5kb.bed

# Convert windows to phylip file
module load python  # if needed

## Prepare per-scaffold window lists
awk '{print > "'$OUT'/windows/by_scaffold/"$1".w5kb.bed"}' $WINBED



# Run this:
#!/usr/bin/env bash
#SBATCH --job-name=PreparePhylips
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/PreparePhylips.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-33 PreparePhylips.sh

BASE=/xdisk/mcnew/dannyjackson/sulidae
OUT=$BASE/analyses/raxml_windows_5kb_snp
CONTIG_FILE=${BASE}/referencelists/AUTOSOMES.txt
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_FILE" | tr -d '\r')
mkdir -p $OUT/windows/by_scaffold $OUT/phy/by_scaffold

export PATH_TO_VCF2PHYLIP=~/programs/vcf2phylip/vcf2phylip.py   # set this!

module load bcftools
module load samtools
module load parallel

WINBED="$OUT/windows/by_scaffold/"$SCAF".w5kb.bed";
VCF="$OUT/merged_scaffold_vcfs/"$SCAF".merged.vcf.gz";
[ -s "$WINBED" ] || exit 0
mkdir -p "$OUT/phy/by_scaffold/"$SCAF

# For each window, extract SNPs and convert
while IFS=$'\t' read -r sc s e _; do
REGION="${sc}:${s}-${e}"
WINVCF="$OUT/phy/by_scaffold/"$SCAF"/${sc}_${s}_${e}.vcf.gz"
bcftools view -r "$REGION" -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' -Oz -o "$WINVCF" "$VCF"
tabix -p vcf "$WINVCF"

# Convert to phylip (relaxed or strict) — SNPs only
python3 $PATH_TO_VCF2PHYLIP --input "$WINVCF" --output-folder "$OUT/phy/by_scaffold/"$SCAF \
    --output-prefix ${sc}_${s}_${e}.phy

# (optional) drop empty windows (no SNPs) to avoid RAxML failures
if [ ! -s "$OUT/phy/by_scaffold/"$SCAF"/${sc}_${s}_${e}.phy" ]; then
    rm -f "$OUT/phy/by_scaffold/"$SCAF"/${sc}_${s}_${e}.phy"
fi

rm "$WINVCF"
rm "$WINVCF.tbi"
done < "$WINBED"

echo done with $SCAF




#######
# Repeat the above with a filtering step
#######

#!/usr/bin/env bash
#SBATCH --job-name=PreparePhylips
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/PreparePhylips.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-33 PreparePhylips.sh
# ---- knobs ----
MISSING_FRAC=0.8    # keep sites with <= 50% missing genotypes
MIN_SNPS=3         # require at least this many SNPs in the window after filtering
# ----------------

BASE=/xdisk/mcnew/dannyjackson/sulidae
OUT=$BASE/analyses/raxml_windows_5kb_snp
CONTIG_FILE=${BASE}/referencelists/AUTOSOMES.txt
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_FILE" | tr -d '\r')
# SCAF=CM062567.1 for testing
export PATH_TO_VCF2PHYLIP=~/programs/vcf2phylip/vcf2phylip.py

module load bcftools
module load samtools
module load parallel

mkdir -p "$OUT/windows/by_scaffold" "$OUT/phy/by_scaffold/$SCAF"

WINBED="$OUT/windows/by_scaffold/${SCAF}.w5kb.bed"
VCF="$OUT/merged_scaffold_vcfs/${SCAF}.merged.vcf.gz"
[ -s "$WINBED" ] || exit 0

while IFS=$'\t' read -r sc s e _; do
  REGION="${sc}:${s}-${e}"

  mkdir -p $OUT/phy/by_scaffold/$SCAF/

  # temp files for this window
  RAWVCF="$OUT/phy/by_scaffold/$SCAF/${sc}_${s}_${e}.raw.vcf.gz"
  FILTVCF="$OUT/phy/by_scaffold/$SCAF/${sc}_${s}_${e}.flt.vcf.gz"
  PHYOUT="$OUT/phy/by_scaffold/$SCAF/${sc}_${s}_${e}.phy"

  # 1) Extract window
  bcftools view -r "$REGION" -Oz -o "$RAWVCF" "$VCF"
  tabix -p vcf "$RAWVCF"

  # 2) ➤ Fill tags and filter sites by MAF>0 (polymorphic among called) AND missingness
  bcftools view \
    -i "(COUNT(GT='./.')/N_SAMPLES) <= ${MISSING_FRAC} && (COUNT(GT='RA')>0 || COUNT(GT='AA')>0)" \
    -Oz -o "$FILTVCF" "$RAWVCF" -m2 -M2 -v snps

  tabix -p vcf "$FILTVCF"

  # 3) ➤ Require a minimum number of SNPs left in the window; otherwise skip the window
  NSNPS=$(bcftools view -H "$FILTVCF" | wc -l)
  if [ "$NSNPS" -lt "$MIN_SNPS" ]; then
    rm -f "$RAWVCF" "$RAWVCF.tbi" "$FILTVCF" "$FILTVCF.tbi"
    echo "Skipping $REGION: only $NSNPS SNPs (need >= $MIN_SNPS)" >> "$OUT/windows/by_scaffold/removedwindows.txt"
    continue
  fi

  # 4) Convert to PHYLIP
  python3 "$PATH_TO_VCF2PHYLIP" --input "$FILTVCF" --output-folder "$OUT/phy/by_scaffold/$SCAF" \
      --output-prefix "${sc}_${s}_${e}.phy"

  # 5) If conversion produced nothing (e.g., all Ns), drop it
  if [ ! -s "$PHYOUT" ]; then
    rm -f "$PHYOUT"
  fi

  # cleanup temps
  rm -f "$RAWVCF" "$RAWVCF.tbi" "$FILTVCF" "$FILTVCF.tbi"
done < "$WINBED"


ls /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/by_scaffold/*/*.phy | wc -l

wc -l /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/windows/by_scaffold/removedwindows.txt



#!/usr/bin/env bash
#SBATCH --job-name=PreparePhylips
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/PreparePhylips.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-33 PreparePhylips.sh
# ---- knobs ----
MISSING_FRAC=0.8    # keep sites with <= 50% missing genotypes
MIN_SNPS=3         # require at least this many SNPs in the window after filtering
# ----------------

BASE=/xdisk/mcnew/dannyjackson/sulidae
OUT=$BASE/analyses/raxml_windows_5kb_snp
CONTIG_FILE=${BASE}/referencelists/AUTOSOMES.txt
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_FILE" | tr -d '\r')
# SCAF=CM062567.1 # for testing

module load bcftools
module load samtools
module load parallel

mkdir -p "$OUT/windows/by_scaffold" "$OUT/phy/by_scaffold/$SCAF"

WINBED="$OUT/windows/by_scaffold/${SCAF}.w5kb.bed"
[ -s "$WINBED" ] || exit 0

# Filter by invariant sites



while IFS=$'\t' read -r sc s e _; do

    echo "Processing window ${sc}:${s}-${e}"

    python ~/programs/filter_invariants_all.py $OUT/phy/by_scaffold/$SCAF/${sc}_${s}_${e}.phy.min4.phy

    mv variantsites.phy $OUT/phy/invariantsites/phylips/${sc}_${s}_${e}.invariant.phy
    mv variantsites_kept.txt  $OUT/phy/invariantsites/infofiles/${sc}_${s}_${e}.variantsites_mind2_kept.txt

    raxml-ng --parse --msa input.phy --model GTR+G --msa-format PHYLIP --prefix filtered --filter

done < "$WINBED"


# Run RAxML 
BASE=/xdisk/mcnew/dannyjackson/sulidae
OUT=$BASE/analyses/raxml_windows_5kb_snp
find /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/by_scaffold/ -name "*.phy" > $OUT/phy/windows.list

pip3 install --user pandas biopython


while read -r PHY; do
    BAS=$(basename "$PHY" .phy)
    python3 ~/programs/raxml_ascbias/ascbias.py -p "$PHY" -o "$OUT/phy/invariantsites/phylips/$BAS.noinv.phy"
done < $OUT/phy/windows.list

find /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/invariantsites/ -name "*.phy" > $OUT/phy/windows.invariant.list

raxml-ng --check \
  --msa /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/invariantsites/phylips/CM062567.1_42100000_42105000.phy.min4.noinv.phy \
  --msa-format PHYLIP \
  --model GTR+G+ASC_FELS{1000} \
  --prefix chk

while read -r PHY; do
    BAS=$(basename "$PHY" .phy)

    raxml-ng --all \
        --msa "$PHY" \
        --model GTR+G+ASC_LEWIS \
        --prefix "$OUT/invariant/trees/"$BAS \
        --threads 2 \
        --bs-trees 20 \
        --seed 13
        
        
done < $OUT/phy/windows.invariant.list

#!/usr/bin/env bash
#SBATCH --job-name=RAxML
#SBATCH --partition=standard
#SBATCH --constraint=hi_mem
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=32
#SBATCH --mem=3008G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/RAxML.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch RAxML.sh

BASE=/xdisk/mcnew/dannyjackson/sulidae
OUT=$BASE/analyses/raxml_windows_5kb_snp

mkdir -p $OUT/windows/by_scaffold $OUT/phy/by_scaffold

module load bcftools
module load samtools
module load parallel

module load raxml-ng

# Run this once, it takes a while -- comment out if you need to do multiple runs with different settings
find $OUT/phy/by_scaffold -name "*.phy" > $OUT/phy/windows.list

parallel --jobs 32 --bar '
  PHY={1}
  BAS=$(basename "$PHY" .phy)
  raxml-ng --all \
    --msa "$PHY" \
    --model GTR+G+ASC_LEWIS \
    --prefix "'$OUT'/trees/"$BAS \
    --threads 2 \
    --bs-trees 20 \
    --seed 13
' :::: $OUT/phy/windows.list



# check status
ls /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/trees/*bestTree | wc -l
wc -l /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/windows.list

### The above was taking too long, so I am rereunning it with a more parallelized approach:
# trying with ParGenes
find . -name "*.phy" -type f -print0 | xargs -0 -I{} cp {} ../all_phy/



cp /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/by_scaffold/*/*phy /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/all_phy/
