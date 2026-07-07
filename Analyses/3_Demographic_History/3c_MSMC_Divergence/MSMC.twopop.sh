

# COMPARISON=BFBO_PEBO  # specify comparison name here
COMPARISON=MABO_NABO  # specify comparison name here

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
# FILES_DIR="${BASE}/${COMPARISON}/files"      # where multihetsep .txt go
OUT_DIR="${BASE}/${COMPARISON}/results"      # where msmc2 outputs go

LIST=${BASE}/${COMPARISON}/samples.txt

while read -r LINE; do 
# RUN_INFO format: SPECIES1,SPECIES2,IND1,IND2

IFS=',' read -r SPECIES1 SPECIES2 IND1 IND2 <<< "$LINE"
echo SPECIES1: $SPECIES1
echo SPECIES2: $SPECIES2
echo IND1: $IND1
echo IND2: $IND2

~/programs/msmc-tools/combineCrossCoal.py ${OUT_DIR}/${IND1}_${IND2}.final.txt "${BASE}/${SPECIES1}/results/${IND1}.CM062567.final.txt" \
    "${BASE}/${SPECIES2}/results/${IND2}.CM062567.final.txt" > "${OUT_DIR}/${IND1}_${IND2}.msmc2.final.txt"

done < "$LIST"

