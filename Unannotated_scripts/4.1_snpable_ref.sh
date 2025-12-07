#!/bin/bash
#SBATCH --job-name=snpable_refmask
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12            # threads per sample
#SBATCH --mem=64G                     # ~5G/core is plenty for bwa-mem2
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/snpable_refmask.%A_%a.out

# script usage:
# sbatch run_snpable2.sh
snpable_script_path=~/programs/seqbility-20091110 # directory with snpable scripts
PATH=$PATH:$snpable_script_path

date

scriptdir=~/programs/msmc2_scripts/
module load gcc
module load msmc2
module load samtools
module load bcftools
module load vcftools
module load bwa
module load python/3.6.3

MSMCTOOLS=~/programs/msmc-tools # folder with msmc-tools binaries
PATH=$PATH:$MSMCTOOLS

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks # main directory for output files

# make directories for intermediate files-- will fail if these don't exist
mkdir -p ${OUTDIR}/vcf
mkdir -p ${OUTDIR}/mask
mkdir -p ${OUTDIR}/input
mkdir ${OUTDIR}/snpable
cd ${OUTDIR}/snpable

GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
prefix=${OUTDIR}/mask/GCA_031468815_ # prefix of genome masks
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams # directory with bamfiles

k=150

echo "Starting extraction of overlapping ${k}-mer subsequences"

~/programs/seqbility-20091110/splitfa $GENOME $k | split -l 20000000
cat x* >> ${prefix}_split.$k

# if it can't find splitfa, try adding seqbility to the path using 'PATH=$PATH:/project/WagnerLab/jrick/msmc_Sept2017/snpable/scripts'

echo "Aligning ${k}-mer reads to the genome with BWA, then converting to sam file"

# the genome needs to be indexed prior to this step-- if it has not already been indexed, run:
if [ -f "${GENOME}.bwt" ]; then
        echo "$GENOME already indexed"
else
        echo "indexing $GENOME"
        bwa index $GENOME
fi

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t ${SLURM_CPUS_PER_TASK} -R 1000000 -O 3 -E 3 ${GENOME} ${prefix}_split.${k} > ${prefix}_split.${k}.sai
bwa samse -f ${prefix}_split.${k}.sam $GENOME ${prefix}_split.${k}.sai ${prefix}_split.${k} 

echo "reads aligned, starting to generate rawMask"
~/programs/seqbility-20091110/gen_raw_mask.pl ${prefix}_split.${k}.sam > ${prefix}_rawMask.${k}.fa

~/programs/seqbility-20091110/gen_raw_mask.pl ${prefix}_split.${k}.sam > ${prefix}_rawMask.${k}.fa


echo "raw mask created as ${prefix}_rawMask.35.fa, now generating final mask with stringency r=50%"
~/programs/seqbility-20091110/gen_mask -l ${k} -r 0.5 ${prefix}_rawMask.${k}.fa > ${prefix}_mask.${k}.50.fa

echo "all done! final mask saved as ${prefix}_mask.${k}.50.fa"

date