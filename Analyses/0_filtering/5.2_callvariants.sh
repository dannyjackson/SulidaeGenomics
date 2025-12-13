#!/usr/bin/env bash
#SBATCH --job-name=callvariants
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_output/callvariants.%A_%a.out
#SBATCH --mail-type=ALL

# LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
# sbatch --array=1-29 callvariants.sh ${LIST}

set -euo pipefail

METHOD="samtools"
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/vcf
METHOD="samtools"
MSMCTOOLS=~/programs/msmc-tools/

IND="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"
GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams
BAMFILE=${BAMDIR}/${IND}.final.bam
CONTIGS=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt

module load samtools
module load bcftools
module load parallel
# Optional: module load mosdepth

mkdir -p "${OUTDIR}"/{mask,vcf,tmp,logs}
COVFILE="${OUTDIR}/coverage_samtoolsDepth_${IND}.txt"

echo -e "\n$(date)  IND=${IND}\nBAM=${BAMFILE}"

# ----- Stage data to fast local scratch -----
SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$SCRATCH"
cp -v "$BAMFILE" "${BAMFILE}.bai" "$SCRATCH/"
cp -v "$GENOME".fai "$SCRATCH/" 2>/dev/null || true
# faidx if needed (once)
[ -f "${GENOME}.fai" ] || samtools faidx "$GENOME"

LOCAL_BAM="$SCRATCH/$(basename "$BAMFILE")"
LOCAL_BAI="${LOCAL_BAM}.bai"

# ----- Threading strategy -----
# bcftools mpileup/call uses threads mostly for I/O; keep small to avoid contention.
PER_TASK_THREADS=4             # threads given to mpileup/call
MAX_JOBS=$(( SLURM_CPUS_PER_TASK / PER_TASK_THREADS ))
(( MAX_JOBS >= 1 )) || MAX_JOBS=1

export GENOME LOCAL_BAM IND OUTDIR METHOD PER_TASK_THREADS
export COVFILE
export PATH  # so parallel workers see loaded modules

# Pre-create coverage file header if you want
: > "$COVFILE"

scaffold_job() {
  s="$1"
  # 1) Mean coverage (fast but single-threaded; parallelization amortizes)
  MEANCOV=$(samtools depth -r "$s" "$LOCAL_BAM" | awk '{sum+=$3} END{if(NR==0)print 0;else printf("%.6f",sum/NR)}')
  printf "%s.%s\t%s\n" "$IND" "$s" "$MEANCOV" >> "$COVFILE"

  MASK_IND="${OUTDIR}/mask/ind_mask.${IND}.${s}.${METHOD}.bed.gz"
  VCF="${OUTDIR}/vcf/${IND}.${s}.${METHOD}.vcf"

  if [ "$METHOD" = "samtools" ]; then
    # mpileup/call with small thread count; keep output uncompressed for speed
    bcftools mpileup -Ou -r "$s" --threads "$PER_TASK_THREADS" -f "$GENOME" "$LOCAL_BAM" \
      | bcftools call -c --threads "$PER_TASK_THREADS" -V indels \
      | "$MSMCTOOLS"/bamCaller.py "$MEANCOV" "$MASK_IND" > "$VCF"
  else
    echo "Unsupported METHOD=$METHOD" >&2
    exit 1
  fi
}

export -f scaffold_job

echo "Starting parallel over scaffolds at $(date). Jobs: $MAX_JOBS, threads/job: $PER_TASK_THREADS"

# GNU parallel: one scaffold per job, up to MAX_JOBS concurrent
parallel --jobs "$MAX_JOBS" --halt soon,fail=1 --joblog "${OUTDIR}/logs/${IND}.parallel.log" \
  scaffold_job :::: "$CONTIGS"

echo "All scaffolds done for ${IND} at $(date). Outputs in ${OUTDIR}/{mask,vcf}"
