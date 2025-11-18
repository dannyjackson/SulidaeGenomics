cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence

chmod +x divergence.generatesamplefile.sh

# generate preliminary files
sbatch divergence.pre.sh BFBO
sbatch divergence.pre.sh PEBO
sbatch divergence.pre.sh MABO
sbatch divergence.pre.sh NABO
sbatch divergence.pre.sh BRBO
sbatch divergence.pre.sh RFBO

# generate individual msmc results
sbatch divergence.indvmsmc.sh BFBO
sbatch divergence.indvmsmc.sh PEBO
sbatch divergence.indvmsmc.sh MABO
sbatch divergence.indvmsmc.sh NABO
sbatch divergence.indvmsmc.sh BRBO
sbatch divergence.indvmsmc.sh RFBO


# BFBO PEBO
./divergence.generatesamplefile.sh BFBO PEBO
sbatch --array=1-25 divergence.2.sh BFBO_PEBO

# 4857399

# MABO NABO
./divergence.generatesamplefile.sh MABO NABO
sbatch --array=1-20 divergence.2.sh MABO_NABO

# 4857415

# BFBOPEBO MABONABO
## BFBO MABO
./divergence.generatesamplefile.sh BFBO MABO
sbatch --array=1-20 divergence.2.sh BFBO_MABO

## BFBO NABO
./divergence.generatesamplefile.sh BFBO NABO
sbatch --array=1-25 divergence.2.sh BFBO_NABO

## PEBO MABO
./divergence.generatesamplefile.sh PEBO MABO
sbatch --array=1-20 divergence.2.sh PEBO_MABO

## PEBO NABO
./divergence.generatesamplefile.sh PEBO NABO
sbatch --array=1-25 divergence.2.sh PEBO_NABO

# BRBO COBO
./divergence.generatesamplefile.sh BRBO COBO
sbatch --array=1-3 divergence.2.sh BRBO_COBO



# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence//BFBO_PEBO/results/BFBO504_PEBO606.final.txt # finished
# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence//BFBO_PEBO/results/BFBO503_PEBO603.final.txt # finished
# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence//BFBO_PEBO/results/BFBO503_PEBO601.final.txt # finished
# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence//MABO_NABO/results/MABO305_NABO402.final.txt # finished

COMPARISON=MABO_NABO

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
FILES_DIR="${BASE}/${COMPARISON}/files"      # where multihetsep .txt go
OUT_DIR="${BASE}/${COMPARISON}/results"      # where msmc2 outputs go
VCFLIST="${BASE}/vcflists"    


MASK_GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/mask/GCA_031468815_CM062567.1.mask.150.50.bed.gz

# Limit to CM062567 + each target individual’s VCFs.
# NOTE: -F (fixed string) and -w (word match) help avoid partial matches like PE10 matching PE1.
# If your filenames do not contain space-delimited tokens, consider dropping -w and keep -F only.
IND1=MABO305
IND2=NABO402
SPECIES1=BFBO
SPECIES2=PEBO
# Multihetsep outputs
COMPARISON_TXT="${FILES_DIR}/${IND1}_${IND2}.CM062567.multihetsep.txt"

# Generate combined multihetsep
~/programs/msmc-tools/generate_multihetsep.py \
  --mask="${MASK_GENOME}" \
  $(cat "${VCFLIST}/${IND1}.txt" "${VCFLIST}/${IND2}.txt") \
  > "${COMPARISON_TXT}"

echo "Running MSMC2 for ${IND1} (${SPECIES1}) and ${IND2} (${SPECIES2})..."

# BRBO_COBO
PEBO_NABO
PEBO_MABO
BFBO_NABO 
BFBO_MABO
BFBO_PEBO
MABO_NABO 

COMPARISON=BRBO_COBO  # specify comparison name here

cd $COMPARISON


BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
# FILES_DIR="${BASE}/${COMPARISON}/files"      # where multihetsep .txt go
OUT_DIR="${BASE}/${COMPARISON}/results"      # where msmc2 outputs go

LIST=${BASE}/${COMPARISON}/samples.txt

while read -r LINE; do 
# RUN_INFO format: SPECIES1,SPECIES2,IND1,IND2

IFS=',' read -r SPECIES1 SPECIES2 IND1 IND2 <<< "$LINE"

# define this for COBO
SPECIES2=BRBO

echo SPECIES1: $SPECIES1
echo SPECIES2: $SPECIES2
echo IND1: $IND1
echo IND2: $IND2

~/programs/msmc-tools/combineCrossCoal.py "${OUT_DIR}/${IND1}_${IND2}.final.txt" "${BASE}/${SPECIES1}/results/${IND1}.CM062567.final.txt" \
    "${BASE}/${SPECIES2}/results/${IND2}.CM062567.final.txt" > "${OUT_DIR}/${IND1}_${IND2}.msmc2.final.txt"

python3 ~/programs/MSMC-IM/MSMC_IM.py -o ${BASE}/${COMPARISON}/MSMC-IM/${IND1}_${IND2} -p 1*2+15*1+1*2 -mu 1.913e-9 "${OUT_DIR}/${IND1}_${IND2}.msmc2.final.txt" --printfittingdetails --plotfittingdetails

done < "$LIST"


# Computing generation time:
G = α + (S/(1-S))
α = age at first breeding, ~4
S = survival rate = 0.91
G = 4 + (0.91/0.09)
G = 14.11
Schreiber, E. A., & Burger, J. (Eds.). (2001). Biology of marine birds. CRC press.

# We used a mutation rate of 4.6 × 10− 9 mutations/site/generation [97], which has been used in other avian systems for PSMC analysis (e.g., [66, 98,99,100]). 
# 1.82E-09
# 3.00E-09
# We assumed a nucleotide substitution rate of m = 1.91*10−9 (Penguin paper)
# Comparative Genomics Supports Ecologically Induced Selection as a Putative Driver of Banded Penguin Diversification 

########################################################
# MSMC-IM
git clone https://github.com/wangke16/MSMC-IM

# BRBO_COBO # in queue
PEBO_NABO # doing
PEBO_MABO # done
BFBO_NABO # done
BFBO_MABO # done
BFBO_PEBO # done
# MABO_NABO # done


# 
COMPARISON=MABO_NABO  # specify comparison name here

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
# FILES_DIR="${BASE}/${COMPARISON}/files"      # where multihetsep .txt go
OUT_DIR="${BASE}/${COMPARISON}/results"      # where msmc2 outputs go

LIST=${BASE}/${COMPARISON}/samples.txt
mkdir -p ${BASE}/${COMPARISON}/MSMC-IM
while read -r LINE; do 
# RUN_INFO format: SPECIES1,SPECIES2,IND1,IND2

IFS=',' read -r SPECIES1 SPECIES2 IND1 IND2 <<< "$LINE"
echo SPECIES1: $SPECIES1
echo SPECIES2: $SPECIES2
echo IND1: $IND1
echo IND2: $IND2


python3 ~/programs/MSMC-IM/MSMC_IM.py -o ${BASE}/${COMPARISON}/MSMC-IM/${IND1}_${IND2} -p 1*2+15*1+1*2 -mu 1.913e-9 "${OUT_DIR}/${IND1}_${IND2}.msmc2.final.txt" --printfittingdetails --plotfittingdetails

done < "$LIST"

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence

python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison BRBO_COBO
python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison PEBO_NABO
python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison PEBO_MABO
python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison BFBO_NABO
python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison BFBO_MABO
python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison BFBO_PEBO
python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/extract_split_quantiles.py --base ${BASE} --comparison MABO_NABO

# BRBO_COBO 1675177.04
# PEBO_NABO 1068119.72
# PEBO_MABO 1065428.08
# BFBO_NABO 1437685.88
# BFBO_MABO 1249922.16
# BFBO_PEBO 1352416.16
# MABO_NABO 1412187.46

COMPARISON=MABO_NABO  # specify comparison name here

tail -n 1 ${COMPARISON}/MSMC-IM/split_times_quantiles.csv


python3 ~/programs/msmc-tools/plot_utils.py results/BFBO501_PEBO601.msmc2.final.txt


# BRBO_COBO 
# PEBO_NABO 
# PEBO_MABO 
# BFBO_NABO 
# BFBO_MABO 
# BFBO_PEBO 
# MABO_NABO 

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
# FILES_DIR="${BASE}/${COMPARISON}/files"      # where multihetsep .txt go
COMPARISON=BRBO_COBO

LIST=${BASE}/${COMPARISON}/samples.txt

cd ${COMPARISON}

while read -r LINE; do 
# RUN_INFO format: SPECIES1,SPECIES2,IND1,IND2

IFS=',' read -r SPECIES1 SPECIES2 IND1 IND2 <<< "$LINE"
echo SPECIES1: $SPECIES1
echo SPECIES2: $SPECIES2
echo IND1: $IND1
echo IND2: $IND2


python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/plot_rccr.py ${IND1}_${IND2}

done < "$LIST"

python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/msmcim_split.py BFBO501_PEBO601 \
  --base /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/BFBO_PEBO \
  --g 14.11



BRBO_COBO 
# PEBO_NABO 
# PEBO_MABO 
# BFBO_NABO 
# BFBO_MABO 
# BFBO_PEBO 
# MABO_NABO 

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/
COMPARISON=BRBO_COBO
LIST=${BASE}/${COMPARISON}/samples.txt

while read -r LINE; do 
IFS=',' read -r SPECIES1 SPECIES2 IND1 IND2 <<< "$LINE"
echo SPECIES1: $SPECIES1
echo SPECIES2: $SPECIES2
echo IND1: $IND1
echo IND2: $IND2

python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/msmcim_split.py ${IND1}_${IND2} \
  --base /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/${SPECIES1}_${SPECIES2} \
  --g 14.11
  
done < "$LIST"


python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/msmcim_split.py BFBO501_PEBO601 \
  --base /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/divergence/BFBO_PEBO \
  --g 14.11