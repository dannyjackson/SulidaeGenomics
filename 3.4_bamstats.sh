
## Compute statistics on bam files 

# compute alignment stats

#!/bin/bash
module load samtools
while read -r bird; do
samtools depth /xdisk/mcnew/dannyjackson/cardinals_dfinch/indelrealignment/$bird.realigned.bam >> /xdisk/mcnew/dannyjackson/cardinals_dfinch/bamstats/$bird_depthstats.txt 
done <  /xdisk/mcnew/dannyjackson/cardinals/sampleids.txt

sbatch --account=mcnew \
--job-name=bamstats \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.bamstats.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/cardinals_dfinch/bamstats.sh 
# Submitted batch job 12031584