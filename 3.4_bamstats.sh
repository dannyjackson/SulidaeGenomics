
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