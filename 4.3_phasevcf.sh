# Phase data

cp ~/programs/msmc2_scripts/SCAFFOLDS.txt .

#!/bin/sh 
module load python

scriptdir=~/programs/msmc2_scripts
source ${scriptdir}/msmc_params.sh

IND=$1
BAMFILE=${BAMDIR}/${IND}.final.bam

echo "working with individual $IND"


for s in `cat ~/programs/msmc2_scripts/SCAFFOLDS.txt`;
    do echo "working with scaffold $s";
    if [ -f ${OUTDIR}/vcf/${IND}.${s}.${phasing}.samtools.vcf.gz ]; then
            echo "phased VCF already exists; moving onto next scaffold";
    else
            echo "phased VCF does not exist; phasing VCF for scaffold $s"
            sed -i 's/^ //g' ${OUTDIR}/vcf2/${IND}.${s}.samtools.vcf

            whatshap phase --reference ${GENOME} --ignore-read-groups -o ${OUTDIR}/vcf3/${IND}.${s}.${phasing}.samtools.vcf.gz ${OUTDIR}/vcf2/${IND}.${s}.samtools.vcf ${BAMFILE}
            whatshap stats --tsv=${OUTDIR}/stats/${IND}.${s}.${prefix}.minDP10.${phasing}.stats.tsv ${OUTDIR}/vcf/${IND}.${s}.${phasing}.samtools.vcf.gz ;
    fi;
done


echo "finished with individual $IND"
date


for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=phasing${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.phasing${i}.%j \
	--nodes=1 \
	--ntasks-per-node=4 \
	--time=25:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/msmc/phasing.sh $i
done

12043660..12043693



# 11152196 - 11152229
for job in {3583203..3583233}; do
    scancel $job
done

