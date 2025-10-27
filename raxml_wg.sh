# RAxML Analysis -- WG

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg



#!/usr/bin/env bash
#SBATCH --job-name=call
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/call.%A_%a.out
#SBATCH --mail-type=ALL

BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
LIST=$BASE/referencelists/allsamplecodes.txt
BAMDIR=$BASE/datafiles/finalbams
OUTDIR=$BASE/analyses/raxml_wg/

module load bcftools/1.19
module load vcftools
module load plink

ID="Sula_MorusBassanus_2"
bcftools mpileup -Ou -f "$REF" -a FORMAT/AD,DP,INFO/AD,SP "$BAMDIR"/*final.bam | bcftools call -mv -V indels > $OUTDIR/vcfs/"$ID"_snps_multiallelic.vcf


grep -v 'JAVH' Sula_MorusBassanus_snps_multiallelic.vcf > Sula_MorusBassanus_snps_multiallelic.contigs.vcf 

bcftools reheader -s /xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt \
    -o Sula_MorusBassanus_snps_multiallelic.contigs.renamed.vcf \
    Sula_MorusBassanus_snps_multiallelic.contigs.vcf 

bcftools view -i 'QUAL>30' Sula_MorusBassanus_snps_multiallelic.contigs.renamed.vcf  > Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf

#filters by depth and removes indels
#vcftools --vcf Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf \
#    --min-meanDP 2 --max-meanDP 8 --remove-indels --recode \
#    --out Sula_MorusBassanus.qualitysort_filtered

vcftools --vcf Sula_MorusBassanus_snps_multiallelic.contigs.qualitysort.vcf \
    --min-meanDP 2 --max-meanDP 20 --remove-indels --recode \
    --out Sula_MorusBassanus.qualitysort_filtered

#!/usr/bin/env bash
#SBATCH --job-name=plink_filter
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/plink_filter.%A_%a.out
#SBATCH --mail-type=ALL
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs

plink --vcf Sula_MorusBassanus.qualitysort_filtered.recode.vcf \
    --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid \
    --out Sula_MorusBassanus.qualitysort_filtered_mind2

bgzip Sula_MorusBassanus.qualitysort_filtered_mind2.vcf
bcftools index Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz

#!/usr/bin/env bash
#SBATCH --job-name=splitvcfchrom
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/splitvcfchrom.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-35 splitvcfchrom.sh 

set -euo pipefail

INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

bcftools view -r "${CHROM}" "${INDIR}"/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz > "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.vcf


# Convert to phylip

#!/usr/bin/env bash
#SBATCH --job-name=chromtophylip
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/chromtophylip.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-35 chromtophylip.sh 

set -euo pipefail

INDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

python ~/programs/vcf2phylip.py3 \
    -i "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.vcf \
    -o "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.phy \
    -r 


python ~/programs/ascbias.py -p "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.phy -o "${OUTDIR}"/"${CHROM}".qualitysort_filtered_mind2.noINV.phy


#!/usr/bin/env bash
#SBATCH --job-name=raxml_chroms
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=94
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/raxml_chroms.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-35 raxml_chroms.sh

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/
INDIR="${OUTDIR}/chroms"
TREEDIR="${OUTDIR}/trees"
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

mkdir -p "${TREEDIR}"

module load raxml-ng

raxml-ng --all \
    --msa "${INDIR}"/"${CHROM}".qualitysort_filtered_mind2.noINV.phy \
    --model GTR+G+ASC_LEWIS \
    --prefix "${TREEDIR}/${CHROM}" \
    --bs-trees 200 \
    --seed 13 \
    --redo



# 68, 80, 85, 600, 610
grep '68' /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt > didntrun.txt
grep '80' /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt >> didntrun.txt
grep '85' /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt >> didntrun.txt
grep '600' /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt >> didntrun.txt

#!/usr/bin/env bash
#SBATCH --job-name=raxml_chroms2
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=94
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/raxml_chroms2.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-4 raxml_chroms2.sh

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/
INDIR="${OUTDIR}/chroms"
TREEDIR="${OUTDIR}/trees"
LIST=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/scripts/didntrun.txt
CHROM="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"

module load raxml-ng

raxml-ng --all \
    --msa "${INDIR}"/"${CHROM}".qualitysort_filtered_mind2.noINV.phy \
    --model GTR+G+ASC_LEWIS \
    --prefix "${TREEDIR}/${CHROM}" \
    --bs-trees 200 \
    --seed 13 \
    --redo

