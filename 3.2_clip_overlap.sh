# clip overlap
# Clip overlap

# run as slurm array

for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=clip_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.clip_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=240:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/clip.sh $i
done

#!/bin/bash
IND=$1

~/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.sorted_RGadded_dupmarked.bam --out /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/${IND}.sorted_RGadded_dupmarked.clipped.bam --stats --params --poolSize 10000000

echo ${IND} >> /xdisk/mcnew/dannyjackson/sulidae/clippingdone.txt

# 11125903 - 11125936



sbatch --account=mcnew \
--job-name=clip_SRR19149562 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.clip_SRR19149562.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=240:00:00 \
/xdisk/mcnew/dannyjackson/sulidae/clip_alt.sh SRR19149562

# Submitted batch job 11126846

#!/bin/bash
IND=$1

~/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.sorted_RGadded_dupmarked.bam --out /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/${IND}.alt.sorted_RGadded_dupmarked.clipped.bam --stats --params --poolSize 10000000

echo ${IND} "alt" >> /xdisk/mcnew/dannyjackson/sulidae/clippingdone.txt

