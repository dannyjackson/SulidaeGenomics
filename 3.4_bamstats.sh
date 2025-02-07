
## Compute statistics on bam files 

# compute alignment stats

#!/bin/bash
module load samtools
while read -r bird; do
echo $bird >> /xdisk/mcnew/dannyjackson/sulidae/bamstats/"$bird"_depthstats.txt 
samtools depth /xdisk/mcnew/dannyjackson/sulidae/indelrealignment/"$bird".final.bam >> /xdisk/mcnew/dannyjackson/sulidae/bamstats/"$bird"_depthstats.txt 
done <  /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt

sbatch --account=mcnew \
--job-name=bamstats \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.bamstats.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/sulidae/bamstats.sh 
# Submitted batch job 12031592


# compute depth per sample 

#!/bin/bash

cd /xdisk/mcnew/dannyjackson/sulidae/bamstats

while read -r bird;
do 
  echo $bird >> depthstats.txt

  # average and standard deviation
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$bird"_depthstats.txt >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt

sbatch --account=mcnew \
--job-name=calcdepth \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.calcdepth.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/sulidae/calcdepth.sh

# Submitted batch job 11131350



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