#!/usr/bin/env bash
#SBATCH --job-name=divergence_pre
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/divergence_pre_%j.out
#SBATCH --mail-type=ALL


set -euo pipefail

# --- Command-line argument handling ---
if [ $# -lt 1 ]; then
  echo "Usage: sbatch divergence.sh <SPECIES>"
  echo "Example: sbatch divergence.pre.sh BFBO"
  exit 1
fi

SPECIES=$1

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
FILES_DIR="${BASE}/${SPECIES}/files"      # where multihetsep .txt go
OUT_DIR="${BASE}/${SPECIES}/results"      # where msmc2 outputs go
VCFLIST="${BASE}/vcflists"    
mkdir -p "$OUT_DIR" "$FILES_DIR" "$VCFLIST" slurm_output
REFLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists

# --- Input files ---
FILE=${REFLIST}/${SPECIES}_samplecodes.txt


cd $BASE/$SPECIES

MASK_GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062567.1.mask.150.50.bed.gz

# Limit to CM062567 + each target individual’s VCFs.
# NOTE: -F (fixed string) and -w (word match) help avoid partial matches like PE10 matching PE1.
# If your filenames do not contain space-delimited tokens, consider dropping -w and keep -F only.
for IND in $(cat "${FILE}"); do
echo $IND
grep "CM062567" "/xdisk/mcnew/dannyjackson/sulidae/referencelists/vcfprocessing/${SPECIES}_phasedvcfs.txt" \
  | grep -F -w "${IND}" > "${VCFLIST}/${IND}.txt"

# Multihetsep outputs
IND_TXT="${FILES_DIR}/${IND}.CM062567.multihetsep.txt"

# Generate per-species multihetsep
echo "Generating multihetsep for ${IND}..."
~/programs/msmc-tools/generate_multihetsep.py \
  --mask="${MASK_GENOME}" \
  $(cat "${VCFLIST}/${IND}.txt") \
  > "${IND_TXT}"
done 
