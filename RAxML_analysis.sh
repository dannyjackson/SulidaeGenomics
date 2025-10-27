#!/usr/bin/env bash
#SBATCH --job-name=IndexVCFS
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/IndexVCFS.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-29 IndexVCFS.sh

# Set paths and inputs
BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
LIST=${BASE}/referencelists/allsamplecodes.txt
IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')

### optional: mask (BED) of callable regions per scaffold, if you want to restrict windows
### e.g. a combined mask per scaffold
MASK_DIR=$BASE/datafiles/snpable_masks/mask

OUT=$BASE/analyses/raxml_windows_5kb_snp
mkdir -p $OUT/logs $OUT/merged_scaffold_vcfs $OUT/windows $OUT/phy $OUT/trees

# Load modules 
module load bcftools

# Build multisample VCFs per scaffold
## list scaffolds by parsing filenames
CONTIG_FILE=${BASE}/referencelists/CONTIGS.txt
### Explanation: ".../IND.SCAFFOLD.whatshap.samtools.vcf.gz" -> SCAFFOLD is field NF-3

## Index all VCFs
echo "Indexing individual $IND"    
while read -r SCAF; do
    echo "Indexing scaffold $SCAF"
    bcftools index "$VCF_DIR"/"$IND"."$SCAF".whatshap.samtools.vcf.gz 
done < $CONTIG_FILE


#!/usr/bin/env bash
#SBATCH --job-name=MergeVCFs
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/MergeVCFs.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-33 MergeVCFs.sh

# Set paths and inputs
BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

### optional: mask (BED) of callable regions per scaffold, if you want to restrict windows
### e.g. a combined mask per scaffold
MASK_DIR=$BASE/datafiles/snpable_masks/mask

OUT=$BASE/analyses/raxml_windows_5kb_snp
mkdir -p $OUT/logs $OUT/merged_scaffold_vcfs $OUT/windows $OUT/phy $OUT/trees

# Load modules 
module load bcftools
module load samtools

# Build multisample VCFs per scaffold
## list scaffolds by parsing filenames
CONTIG_FILE=${BASE}/referencelists/AUTOSOMES.txt

SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_FILE" | tr -d '\r')

## Merge per scaffold

# gather all individuals’ VCFs for this scaffold
ls "$VCF_DIR"/*."$SCAF".whatshap.samtools.vcf.gz > $OUT/tmp.list.$SCAF

# normalize/merge into a multi-sample VCF
bcftools merge -m none -Oz -o $OUT/merged_scaffold_vcfs/$SCAF.merged.vcf.gz \
-l $OUT/tmp.list.$SCAF

tabix -p vcf $OUT/merged_scaffold_vcfs/$SCAF.merged.vcf.gz

rm $OUT/tmp.list.$SCAF

### Optional filtering step


#!/usr/bin/env bash
#SBATCH --job-name=FilterVCFs
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/FilterVCFs.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-33 FilterVCFs.sh

# Set paths and inputs
BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
CONTIG_FILE=${BASE}/referencelists/AUTOSOMES.txt

### optional: mask (BED) of callable regions per scaffold, if you want to restrict windows
### e.g. a combined mask per scaffold
MASK_DIR=$BASE/datafiles/snpable_masks/mask

OUT=$BASE/analyses/raxml_windows_5kb_snp
mkdir -p $OUT/logs $OUT/merged_scaffold_vcfs $OUT/windows $OUT/phy $OUT/trees

# Load modules 
module load bcftools
module load samtools
module load plink

SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONTIG_FILE" | tr -d '\r')

bcftools view -m2 -M2 -v snps -i 'QUAL>30' \
-Oz -o $OUT/merged_scaffold_vcfs/$SCAF.filtered.vcf.gz \
$OUT/merged_scaffold_vcfs/$SCAF.merged.vcf.gz

tabix -p vcf $OUT/merged_scaffold_vcfs/$SCAF.filtered.vcf.gz





ls merged_scaffold_vcfs/*.filtered.vcf.gz > filtered_scaff_vcf_list.txt
bcftools concat -f filtered_scaff_vcf_list.txt -Oz -o all_chromosomes.filtered.vcf.gz
bcftools index -t all_chromosomes.filtered.vcf.gz


# filter vcfs 

plink \
  --vcf all_chromosomes.filtered.vcf.gz \
  --allow-extra-chr \
  --snps-only just-acgt \
  --geno 0.4 \
  --recode vcf-iid \
  --out all_chromosomes.filtered.geno0.4

bgzip -c all_chromosomes.filtered.geno0.4.min4.vcf > all_chromosomes.filtered.geno0.4.min4.vcf.gz

tabix -p vcf all_chromosomes.filtered.geno0.4.min4.vcf.gz

/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.merged.vcf.gz

/home/u15/dannyjackson/programs/raxml-ng --all \
    --msa /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.merged.vcf.gz \
    --model GTR+G \
    --prefix /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/test/CM062567 \
    --threads 1 \
    --bs-trees 20 \
    --seed 13

raxml-ng --all \
    --msa /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/all_chromosomes.filtered.mac2.geno0.2.min4.phy \
    --model GTR+G \
    --prefix /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/wgs_tree/all_chromosomes.filtered.mac2.geno0.2.min4 \
    --threads 1 \
    --bs-trees 20 \
    --seed 13



#!/usr/bin/env bash
#SBATCH --job-name=cellphy_test
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/cellphy_test.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch cellphy_test.sh
# compiled on puma

module load plink

plink2 /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.filtered.vcf.gz --indep-pairwise 50 5 0 
module load R 


# Prune vcf on ld


plink --vcf /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.filtered.vcf.gz \
--set-all-var-ids @:#  --recode vcf-iid --out tmp_varID

plink --vcf tmp_varID.vcf \
 --rm-dup force-first   --recode vcf-iid --out tmp_nodup

plink \
  --vcf tmp_nodup.vcf \
  --indep-pairwise 50 5 0.6 \
  --allow-extra-chr \
  --snps-only just-acgt \
  --recode vcf-iid \
  --out CM062567.LDpruned \
  --bad-ld 




#!/usr/bin/env bash
#SBATCH --job-name=cellphy_test
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=470G
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/cellphy_test.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch cellphy_test.sh
# compiled on puma

module purge
module load R/4.4.0

/home/u15/dannyjackson/programs/cellphy/cellphy.sh RAXML \
    --bootstrap \
    --msa /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.filtered.vcf.gz \
    --prefix /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/test/CM062567_cellphy.interactive \
    --model GTR+G \
    --bs-trees 20 \
    --seed 13 

# Convert to phylip (relaxed or strict) — SNPs only
python3 ~/programs/vcf2phylip/vcf2phylip.py --input /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.merged.vcf.gz \
  --output-folder /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy \
  --output-prefix CM062567

SNPVCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/CM062567.1.merged.vcf.gz
python3 ~/programs/vcf2phylip/vcf2phylip.py \
  --input "$SNPVCF" \
  --output-folder /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy \
  --output-prefix CM062567 \
  -f -r
  
  --output-fasta \
  --snps-only


raxml-ng --all \
    --msa /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/phy/all_chromosomes.filtered.mac2.geno0.2.min4.phy \
    --model GTR+G \
    --prefix /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/wgs_tree/all_chromosomes.filtered.mac2.geno0.2.min4 \
    --threads 1 \
    --bs-trees 20 \
    --seed 13


raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n sula_flightless_b -q /data5/sulidae/final/to_flightless/raxml/partition_file.txt -s /data5/sulidae/my_datasets/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V








# RAxML


#!/usr/bin/env bash
#SBATCH --job-name=vcf2phylip
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/vcf2phylip.%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp

~/programs/vcf2phylip/vcf2phylip.py -i all_chromosomes.filtered.geno0.4.vcf

#!/usr/bin/env bash
#SBATCH --job-name=generateconsensus
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/generateconsensus.%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/generate_consensus

python ~/programs/generate-consensus.py \
  --ref /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna \
  --vcf /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/generate_consensus/CM062567.1.filtered.vcf.gz \
  --format phylip \
  --variants

/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/generate_consensus/CM062567.1.filtered.vcf.gz
python ~/programs/filter_invariants_all.py all_chromosomes.filtered.geno0.4.phy
mv variantsites.phy all_chromosomes.filtered.geno0.4.invariant.phy


raxml-ng --check \
  --msa ../all_chromosomes.filtered.geno0.4.phy \
  --msa-format PHYLIP \
  --model \
  --prefix chk

#!/usr/bin/env bash
#SBATCH --job-name=RAxML
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/RAxML.%A_%a.out
#SBATCH --mail-type=ALL

module load raxml-ng

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/analysis

raxml-ng --all \
    --msa ../all_chromosomes.filtered.geno0.4.min4.phy \
    --model GTR+G \
    --prefix all_chr \
    --threads 32 \
    --bs-trees 20 \
    --seed 13
    

