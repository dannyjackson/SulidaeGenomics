#!/usr/bin/env bash
#SBATCH --job-name=divergence_indvmsmc
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --output=slurm_output/divergence_indvmsmc_%j.out
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
REFLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists
FILE=${REFLIST}/${SPECIES}_samplecodes.txt

for IND in $(cat "${FILE}"); do
echo $IND

IND_MSMC="${OUT_DIR}/${IND}.CM062567"
IND_TXT="${FILES_DIR}/${IND}.CM062567.multihetsep.txt"

~/programs/msmc_2.0.0_linux64bit \
  -p 1*2+15*1+1*2 \
  -i 10 \
  -o ${IND_MSMC} \
  -I 0,1 \
  ${IND_TXT}

done 