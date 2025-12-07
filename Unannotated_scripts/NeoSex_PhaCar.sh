# Neo Sex: Inference via alignment to PhaCar



#!/bin/bash
#SBATCH --job-name=align
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12            # threads per sample
#SBATCH --mem=64G                     # ~5G/core is plenty for bwa-mem2
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/align.%A_%a.out

set -euo pipefail

module load bwa
module load samtools
module list

LIST=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
TRIM=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq
OUT=/xdisk/mcnew/dannyjackson/sulidae/datafiles/PhaCar_alignment/bamfiles

SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$OUT" "$SCRATCH" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
echo "[$(date)] aligning $SAMPLE with ${SLURM_CPUS_PER_TASK} threads"

# Use node-local scratch for sort temp; then move result
TMPDIR="$SCRATCH"

bwa mem -t ${SLURM_CPUS_PER_TASK} ${REF} /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${SAMPLE}_trimmed_1P.fq.gz \
/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${SAMPLE}_trimmed_2P.fq.gz \
| samtools sort -@ $((SLURM_CPUS_PER_TASK-2)) -m 3G -T "$SCRATCH/${SAMPLE}.tmp" -o "$OUT/${SAMPLE}.bam" -

samtools index -@ 2 "$OUT/${SAMPLE}.bam"
echo "[$(date)] done $SAMPLE"


sbatch --array=1-34%10 align.sh


#!/usr/bin/env bash

conv="/xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt"

while read -r old new; do
  for f in ${old}*; do
    if [[ -e "$f" ]]; then
      newname="${f/$old/$new}"
      mv "$f" "$newname"
    fi
  done
done < "$conv"

# Make VCF
#!/usr/bin/env bash
#SBATCH --job-name=call
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=94
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/call.%A_%a.out
#SBATCH --mail-type=ALL

BASE=/xdisk/mcnew/dannyjackson/sulidae
REF=$BASE/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
LIST=$BASE/referencelists/allsamplecodes.txt
BAMDIR=$BASE/datafiles/PhaCar_alignment/bamfiles

OUTDIR=$BASE/datafiles/bamstats/sexing_by_heterozygosity_filtered/unfiltered_bams/

module load bcftools/1.19
module load vcftools
module load plink

THREADS=${SLURM_CPUS_PER_TASK:-1}
ID="Sula_PhaCar_unfilteredbams"


bcftools mpileup \
  -Ou \
  -f "$REF" \
  -a FORMAT/AD,DP,INFO/AD,SP \
  --threads "$THREADS" \
  "$BAMDIR"/*.bam \
| bcftools call \
  -mv -V indels \
  --threads "$THREADS" \
  > "$OUTDIR"/vcfs/"$ID".vcf

