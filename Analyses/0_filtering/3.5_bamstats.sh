
## Compute statistics on bam files 



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




