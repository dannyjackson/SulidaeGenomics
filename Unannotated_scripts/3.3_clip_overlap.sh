# clip overlap
# Clip overlap

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/

BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/sortedbams

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

~/programs/bamUtil-master/bin/bam clipOverlap \
--in /xdisk/mcnew/dannyjackson/sulidae/datafiles/sortedbams/${IND}.sorted_RGadded_dupmarked.bam \
--out /xdisk/mcnew/dannyjackson/sulidae/datafiles/clipoverlap/${IND}.sorted_RGadded_dupmarked.clipped.bam \
--stats --params --poolSize 10000000

echo ${IND} >> /xdisk/mcnew/dannyjackson/sulidae/clippingdone.txt


#!/bin/bash
#SBATCH --job-name=clipoverlap
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/clipoverlap.%A_%a.out

set -euo pipefail

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
SORTBAM=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/
SCRATCH="${TMPDIR:-/scratch/$USER/$SLURM_JOB_ID}"
mkdir -p "$OUTDIR" "$SCRATCH" slurm_output

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')
INFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.bam
OUTFILE=${SORTBAM}/${SAMPLE}.sorted_RGadded_dupmarked.clip.bam
echo "[$(date)] aligning $SAMPLE with ${SLURM_CPUS_PER_TASK} threads"

~/programs/bamUtil-master/bin/bam clipOverlap \
--in ${INFILE} \
--out ${OUTFILE} \
--stats --params --poolSize 10000000

echo ${SAMPLE} >> /xdisk/mcnew/dannyjackson/sulidae/clippingdone.txt

sbatch --array=1-29%10 clip.sh

INFILE=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/BFBO501.sorted_RGadded_dupmarked.bam
OUTFILE=/xdisk/mcnew/dannyjackson/sulidae/bamfiles/BFBO501.sorted_RGadded_dupmarked.bam

~/programs/bamUtil-master/bin/bam clipOverlap \
--in ${INFILE} \
--out ${OUTFILE} \
--stats --params --poolSize 10000000



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

