#!/usr/bin/env bash
#SBATCH --job-name=msmc_sets
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-100%10
#SBATCH --output=slurm_out/msmc_sets_%A_%a.out

set -euo pipefail

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/divergence/test
# --- inputs you already have ---
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/divergencepairs.txt
ALL=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

# --- dirs you likely use (edit if needed) ---
BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc
FILES_DIR=$BASE/files                      # where multihetsep / txt go
SPECIES_DIR=$BASE/files/species            # where per-set stems live
OUT_DIR=$BASE/files/divergence/test                          # where msmc2 final outputs go

mkdir -p "$OUT_DIR" "${OUT_DIR}/inputfiles"

# 1) Read the Nth line = divergence set (e.g., "BFBO PEBO MABO")
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r' | xargs)
# testing:
# line=$(sed -n 1p $LIST | tr -d '\r' | xargs)

if [[ -z "${line}" ]]; then
  echo "Empty line for task ${SLURM_ARRAY_TASK_ID}; nothing to do."
  exit 0
fi

# 2) Split into codes array
read -r -a codes <<< "$line"

# 3) Label like BFBO_PEBO_MABO
IND=$(IFS=_; printf '%s' "${codes[*]}")

# 4) Regex that matches any code at start of line
regex="^($(printf '%s|' "${codes[@]}" | sed 's/|$//'))"

# 5) Count and write sample list
NR_IND=$(grep -E -c "$regex" "$ALL")
grep -E "$regex" "$ALL" | awk '{print $1}' > "${OUT_DIR}/${IND}.samples"

# 6) Stems/paths you can plug into your existing MSMC steps
MSMC_STEM="${SPECIES_DIR}/${IND}_allchr"       # e.g., input stem you referenced
MSMC_OUT="${OUT_DIR}/${IND}.final"             # final msmc output prefix
MSMC_IN="${OUT_DIR}/inputfiles/${IND}.final"             # final msmc output prefix

echo "SET:         $line"
echo "IND:         $IND"
echo "NR_IND:      $NR_IND"
echo "SAMPLES:     ${OUT_DIR}/${IND}.samples"
echo "MSMC_STEM:   $MSMC_STEM"
echo "MSMC_OUT:    $MSMC_OUT"
echo "MSMC_IN:    $MSMC_IN"

# Make list of vcfs for input into each divergence set's multihetsep generation
VCFPROC=/xdisk/mcnew/dannyjackson/sulidae/referencelists/vcfprocessing
VCFLIST="${VCFPROC}/${IND}_phasedvcfs.txt"

# clear or create it fresh
: > "$VCFLIST"

echo "Building $VCFLIST ..."
for sp in "${codes[@]}"; do
    SRC="${VCFPROC}/${sp}_phasedvcfs.txt"
    if [[ -f "$SRC" ]]; then
        cat "$SRC" >> "$VCFLIST"
    else
        echo "Warning: missing $SRC" >&2
    fi
done


MASK_GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062587.1.mask.150.50.bed.gz
OUTFILE=${OUT_DIR}/${IND}.all.txt

PEBO_OUTFILE=PEBO.CM062587.txt
BFBO_OUTFILE=BFBO.CM062587.txt


# test 

grep "CM062587" /xdisk/mcnew/dannyjackson/sulidae/referencelists/vcfprocessing/BFBO_phasedvcfs.txt | head -n 2 > BFBO_test.txt
grep "CM062587" /xdisk/mcnew/dannyjackson/sulidae/referencelists/vcfprocessing/PEBO_phasedvcfs.txt | head -n 2 > PEBO_test.txt

~/programs/msmc-tools/generate_multihetsep.py --mask=$MASK_GENOME \
    `cat BFBO_test.txt PEBO_test.txt` > test.BFBOPEBO.txt
# 8 haplotypes


~/programs/msmc-tools/generate_multihetsep.py --mask=$MASK_GENOME \
    `cat BFBO_test.txt` > test.BFBO.txt

~/programs/msmc-tools/generate_multihetsep.py --mask=$MASK_GENOME \
    `cat PEBO_test.txt` > test.PEBO.txt


~/programs/msmc_2.0.0_linux64bit \
  -p 1*2+15*1+1*2 \
  -i 10 \
  -o BFBO_21 \
  -I 0,1,2,3 \
  test.BFBOPEBO.txt


~/programs/msmc_2.0.0_linux64bit \
  -p 1*2+15*1+1*2 \
  -i 10 \
  -o PEBO_21 \
  -I 4,5,6,7 \
  test.BFBOPEBO.txt


~/programs/msmc_2.0.0_linux64bit \
  -p 1*2+15*1+1*2 \
  -i 10 \
  -o BFBO_PEBO_21 \
  -I 0,1,4,5 \
  -P 0,0,1,1 \
  test.BFBOPEBO.txt

~/programs/msmc-tools/combineCrossCoal.py BFBO_PEBO_21.final.txt BFBO_21.final.txt \
    PEBO_21.final.txt > BFBO_PEBO.combined.msmc2.final.txt

mu <- 1.25e-8
gen <- 30
crossPopDat<-read.table("BFBO_PEBO.combined.msmc2.final.txt", header=TRUE)
pdf("crosspop.pdf", width = 6, height = 6)  # higher resolution and larger canvas

plot(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     xlim=c(1000,500000),ylim=c(0,1), type="n", xlab="Years ago", ylab="relative cross-coalescence rate")
lines(crossPopDat$left_time_boundary/mu*gen, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11), type="s")
dev.off()

