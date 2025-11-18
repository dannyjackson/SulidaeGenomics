#!/usr/bin/env bash
# Usage: ./make_samples_list.sh BFBO PEBO
# Output: creates ${BASE}/BFBO_PEBO/samples.txt

set -euo pipefail

# --- Arguments ---
if [ $# -ne 2 ]; then
  echo "Usage: ./make_samples_list.sh <SPECIES1> <SPECIES2>"
  echo "Example: ./make_samples_list.sh BFBO PEBO"
  exit 1
fi

SPECIES1=$1
SPECIES2=$2

# --- Paths ---
BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence
REFLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists
OUTDIR=${BASE}/${SPECIES1}_${SPECIES2}
LIST=${OUTDIR}/samples.txt
mkdir -p "$OUTDIR"

# --- Input files ---
FILE1=${REFLIST}/${SPECIES1}_samplecodes.txt
FILE2=${REFLIST}/${SPECIES2}_samplecodes.txt

if [[ ! -f "$FILE1" || ! -f "$FILE2" ]]; then
  echo "Error: missing samplecode list(s):"
  echo "  $FILE1"
  echo "  $FILE2"
  exit 1
fi

# --- Generate combinations ---
echo "Creating $LIST ..."
> "$LIST"  # clear existing

while read -r IND1; do
  [[ -z "$IND1" ]] && continue
  while read -r IND2; do
    [[ -z "$IND2" ]] && continue
    echo "${SPECIES1},${SPECIES2},${IND1},${IND2}" >> "$LIST"
  done < "$FILE2"
done < "$FILE1"

echo "Done. Wrote $(wc -l < "$LIST") comparisons to $LIST"
