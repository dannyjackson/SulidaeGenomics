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

# sbatch --array=1-35 1_raxml_chroms.sh

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