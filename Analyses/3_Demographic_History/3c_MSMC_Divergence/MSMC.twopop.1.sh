#!/usr/bin/env bash
#SBATCH --job-name=MSMC_run
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/MSMC_run_%A_%a.out
#SBATCH --mail-type=ALL
# sbatch --array=1-N divergence.sh BFBO_PEBO

set -euo pipefail

# --- Command-line argument handling ---
if [ $# -lt 1 ]; then
  echo "Usage: sbatch divergence.sh <COMPARISON>"
  echo "Example: sbatch divergence.sh BFBO_PEBO"
  exit 1
fi

COMPARISON=$1

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
FILES_DIR="${BASE}/${COMPARISON}/files"      # where multihetsep .txt go
OUT_DIR="${BASE}/${COMPARISON}/results"      # where msmc2 outputs go
VCFLIST="${BASE}/vcflists"    
mkdir -p "$OUT_DIR" "$FILES_DIR" "$VCFLIST" slurm_output

LIST=${BASE}/${COMPARISON}/samples.txt
RUN_INFO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
# RUN_INFO format: SPECIES1,SPECIES2,IND1,IND2
IFS=',' read -r SPECIES1 SPECIES2 IND1 IND2 <<< "$RUN_INFO"

cd $BASE/$COMPARISON

MASK_GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062567.1.mask.150.50.bed.gz


# Limit to CM062567 + each target individual’s VCFs.
# NOTE: -F (fixed string) and -w (word match) help avoid partial matches like PE10 matching PE1.
# If your filenames do not contain space-delimited tokens, consider dropping -w and keep -F only.

# Multihetsep outputs
COMPARISON_TXT="${OUT_DIR}/${IND1}_${IND2}.CM062567.multihetsep.txt"


echo "Running MSMC2 for ${IND1} (${SPECIES1}) and ${IND2} (${SPECIES2})..."
# Run MSMC2 on combined
~/programs/msmc_2.0.0_linux64bit \
  -p 1*2+15*1+1*2 \
  -i 10 \
  -o ${OUT_DIR}/${IND1}_${IND2} \
  -I 0,1,2,3 \
  -P 0,0,1,1 \
  -s \
  ${COMPARISON_TXT}

echo "Combining cross-coalescence results for ${IND1} (${SPECIES1}) and ${IND2} (${SPECIES2})..."
# Analyze with individual runs
~/programs/msmc-tools/combineCrossCoal.py ${OUT_DIR}/${IND1}_${IND2}.final.txt "${BASE}/${SPECIES1}/results/${IND1}.CM062567.final.txt" \
    "${BASE}/${SPECIES2}/results/${IND2}.CM062567.final.txt" > "${OUT_DIR}/${IND1}_${IND2}.msmc2.final.txt"