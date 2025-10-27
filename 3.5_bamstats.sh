
## Compute statistics on bam files 

# compute alignment stats


# compute depth per sample 

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/

#!/bin/bash
#SBATCH --job-name=bamstats
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/bamstats.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
BAMFILE=${BAMDIR}/${SAMPLE}.final.bam
OUTFILE=${OUTDIR}/${SAMPLE}_depthstats.txt
STATSFILE=${OUTDIR}/alldepthstats.txt
module load samtools
samtools depth "${BAMDIR}/${SAMPLE}.final.bam" >> "$OUTFILE"

echo $SAMPLE >> ${STATSFILE}

# average and standard deviation
awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' ${OUTFILE} >> ${STATSFILE}

sbatch --array=1-29%10 bamstats.sh




while read -r bird;
do 
  echo $bird >> depthstats.txt

  # average and standard deviation
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$bird"_depthstats.txt >> depthstats.txt

done < ${LIST}

sbatch --account=mcnew \
--job-name=calcdepth \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.calcdepth.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=48:00:00 \
calcdepth.sh

# Submitted batch job 11131350

#!/usr/bin/env bash
# Parallel mean/sd over column 3 of *_depthstats.txt files
# Output: depthstats.tsv with columns: bird<TAB>depth<TAB>stdev

set -euo pipefail

module load parallel

# ---- Config (edit if needed) ----
LIST=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt
STATDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth
OUT=${STATDIR}/depthstats.parallel.tsv

# Use all allocated CPUs if under SLURM, otherwise all local cores
JOBS="${SLURM_CPUS_ON_NODE:-$(nproc)}"

cd "$STATDIR"

# Header
printf "bird\tdepth\tstdev\n" > "$OUT"

# Speed tips: C locale for faster awk; guard against small negative var from FP errors
export LC_ALL=C

# Run in parallel; keep output in the same order as LIST (-k)
parallel -k --jobs "$JOBS" --no-run-if-empty '
  bird={}
  f="'"$STATDIR"'/${bird}_depthstats.txt"

  if [[ -s "$f" ]]; then
    awk -v bird="$bird" '"'"'
      { sum += $3; sumsq += $3*$3 }
      END {
        n = NR
        if (n > 0) {
          mean = sum / n
          var  = (sumsq / n) - (mean * mean)
          if (var < 0) var = 0  # FP guard
          sd = sqrt(var)
          printf "%s\t%.6f\t%.6f\n", bird, mean, sd
        } else {
          printf "%s\tNA\tNA\n", bird
        }
      }
    '"'"' "$f"
  else
    # Missing or empty file → NA row (still emit a line to keep order)
    printf "%s\tNA\tNA\n" "$bird"
  fi
' :::: "$LIST" >> "$OUT"

echo "Wrote $OUT"

sbatch --account=mcnew \
--job-name=calcdepthparallel \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.calcdepthparallel.%j \
--nodes=1 \
--ntasks-per-node=1 \
--cpus-per-task=12 \
--time=6:00:00 \
calcdepth.parallel.sh


# less than 3 average depth
MABO301 NABO401 BFBO506, PEBO602
# greater than 50 average Stdev
BRBO204

cat depthstats.txt 
RFBO101
Average =  6.44044
Stdev =  20.5312
RFBO102
Average =  7.42855
Stdev =  22.5939
RFBO103
Average =  13.8523
Stdev =  40.189
RFBO104
Average =  4.93989
Stdev =  14.66
RFBO105
Average =  13.2342
Stdev =  37.292
RFBO106
Average =  5.26899
Stdev =  20.1239
BRBO201
Average =  5.72509
Stdev =  26.0884
BRBO202
Average =  5.23748
Stdev =  32.5151
BRBO203
Average =  4.00899
Stdev =  28.1569
BRBO204
Average =  3.19737
Stdev =  162.598
BRBO205
Average =  3.72565
Stdev =  37.3861
MABO301 
Average =  2.90332
Stdev =  19.414
MABO302
Average =  5.42676
Stdev =  18.7655
MABO304
Average =  4.45272
Stdev =  32.3557
MABO305
Average =  4.12955
Stdev =  32.6736
MABO306
Average =  4.11271
Stdev =  25.1543
NABO401 
Average =  2.57944
Stdev =  26.273
NABO402
Average =  3.9151
Stdev =  29.6072
NABO403
Average =  5.51654
Stdev =  25.4479
NABO404
Average =  5.53378
Stdev =  22.5798
NABO405
Average =  5.43561
Stdev =  29.529
NABO406
Average =  5.11753
Stdev =  29.8089
BFBO501
Average =  5.12754
Stdev =  34.996
BFBO502
Average =  5.39171
Stdev =  31.2432
BFBO503
Average =  4.6679
Stdev =  26.2793
BFBO504
Average =  3.89637
Stdev =  27.9247
BFBO505
Average =  5.39507
Stdev =  31.5529
BFBO506
Average =  2.41157
Stdev =  42.5753
PEBO601
Average =  4.15249
Stdev =  18.0643
PEBO602
Average =  2.58575
Stdev =  25.1816
PEBO603
Average =  4.40899
Stdev =  36.3871
PEBO604
Average =  5.19002
Stdev =  36.9476
PEBO605
Average =  4.38301
Stdev =  35.3789
PEBO606
Average =  4.456
Stdev =  31.5736


IND=BFBO501
DEPTHFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/${IND}_depthstats.txt
OUTFILE=${IND}_coverage_by_scaffold.txt

awk '
{sum[$1]+=$3; n[$1]++}
END {for (scaf in sum) print scaf, sum[scaf]/n[scaf]}
' "$DEPTHFILE" > "$OUTFILE"


LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

while read -r IND;
do 
  DEPTHFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/${IND}_depthstats.txt
  OUTFILE=${IND}_coverage_by_scaffold.txt

  awk '
  {sum[$1]+=$3; n[$1]++}
  END {for (scaf in sum) print scaf, sum[scaf]/n[scaf]}
  ' "$DEPTHFILE" > "$OUTFILE"
done < ${LIST}

grep 'CM062600.1' /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/BFBO501_depthstats.txt | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` # calculate mean coverage



LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
BAMSTATS=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats
CONTIGS=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt

while read -r IND; do
  IN="$BAMSTATS/${IND}_depthstats.txt"
  while read -r SCAF; do
    OUT="${IND}_${SCAF}_mean_coverage.txt"
    LC_ALL=C awk -v sc="$SCAF" '$1==sc{sum+=$3; n++} END{if(n) printf("%.6f\n", sum/n); else print 0}' "$IN" > "$OUT"
  done < $CONTIGS
done < "$LIST"


#!/usr/bin/env bash
#SBATCH --job-name=contigdepth
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/contigdepth.%A_%a.out
#SBATCH --mail-type=ALL

set -euo pipefail

module load parallel

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
BAMSTATS=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats
CONTIGS=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
JOBS=${SLURM_CPUS_PER_TASK:-8}   # adjust or set via SLURM

export BAMSTATS
LC_ALL=C parallel --bar --jobs "$JOBS" '
  IN="'$BAMSTATS'"/{1}_depthstats.txt
  OUT="'$BAMSTATS'"/mean_by_scaffold/{1}_{2}_mean_coverage.txt
  awk -v sc="{2}" '"'"'$1==sc{sum+=$3; n++} END{if(n) printf("%.6f\n", sum/n); else print 0}'"'"' "$IN" > "$OUT"
' :::: "$LIST" :::: "$CONTIGS"


while read -r IND;
 do 
  mv $IND.final.bai $IND.final.bam.bai 
done < ${LIST}


# run once per IND before the parallel loop
IND=BFBO501
DEPTHFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/${IND}_depthstats.txt
COVMAP=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/mean_by_scaffold/${IND}.mean_by_scaffold.tsv

# scaffold\tmean_depth
awk '{sum[$1]+=$3; n[$1]++} END{for(s in sum) printf "%s\t%.6f\n", s, sum[s]/n[s]}' "$DEPTHFILE" \
  > "$COVMAP"




