# indel realignment
# Make dictionary for reference genome, necessary for doing indel realignment

REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
DICT=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.dict

module load picard/2.23.4 

picard CreateSequenceDictionary -R ${REF} -O ${DICT}


#!/bin/bash
#SBATCH --job-name=indexclip
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/indexclip.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
SORTBAM=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/
OUTDIR=${SORTBAM}/indelmaps

mkdir -p "$OUTDIR" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')

BAMFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.clip.bam

module load samtools
samtools index ${BAMFILE} -@ ${SLURM_CPUS_PER_TASK}

sbatch --array=1-29%10 indexclip.sh




# Make indel maps

#!/bin/bash
#SBATCH --job-name=indelrealignment
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/indelrealignment.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
SORTBAM=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/
OUTDIR=${SORTBAM}/indelmaps

mkdir -p "$OUTDIR" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
BAMFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.clip.bam
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R ${REF} \
-I ${BAMFILE} \
-o ${OUTDIR}/${SAMPLE}.intervals

sbatch --array=1-29%10 indelrealignment.sh

#!/bin/bash
#SBATCH --job-name=indelrealignment2
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/indelrealignment2.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
SORTBAM=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/
OUTDIR=${SORTBAM}/indelmaps

mkdir -p "$OUTDIR" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
BAMFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.clip.bam
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R ${REF} \
  --consensusDeterminationModel USE_READS \
  -I ${BAMFILE} \
  --targetIntervals ${OUTDIR}/${SAMPLE}.intervals \
  -o ${SORTBAM}/${SAMPLE}.final.bam

sbatch --array=1-29%10 indelrealignment2.sh

mv /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamfiles/*final* /xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams/