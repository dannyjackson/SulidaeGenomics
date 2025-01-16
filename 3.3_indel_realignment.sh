# indel realignment
# Make dictionary for reference genome, necessary for doing indel realignment

cd ..
mkdir indelmaps
cd indelmaps 


java -jar picard.jar Create

~/programs/jdk1.8.0_411/bin/java -jar ~/programs/picard.jar CreateSequenceDictionary -R GCF_963921805.1_bPhaCar2.1_genomic.fna -O GCF_963921805.1_bPhaCar2.1_genomic.dict

for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=indexclip_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.indexclip_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/index_clip.sh $i
done

# 11127239 - 11127272
#!/bin/bash
IND=$1
module load samtools
samtools index /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/$IND.sorted_RGadded_dupmarked.clipped.bam
echo "done " ${IND} >>/xdisk/mcnew/dannyjackson/sulidae/index_clippedstats.txt 

# lagging indv
sbatch --account=mcnew \
--job-name=indexclip_SRR19149562 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.indexclip_SRR19149562.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=5:00:00 \
/xdisk/mcnew/dannyjackson/sulidae/index_clip.sh SRR19149562
# Submitted batch job 11127430

# Make indel maps
# make indel maps slurm array
11127362
for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=indelrealign_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.indelrealign_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=48:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/indelmap.sh $i
done

# lagging indv
sbatch --account=mcnew \
--job-name=indelrealign_SRR19149562 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.indelrealign_SRR19149562.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/sulidae/indelmap.sh SRR19149562
# Submitted batch job 12031204

#!/bin/bash
IND=$1

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
-I /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/${IND}.sorted_RGadded_dupmarked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/sulidae/indelmaps/${IND}.intervals


for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=indelrealign_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.indelrealign2_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=48:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/indelrealign.sh $i
done
12031173
#!/bin/bash
IND=$1

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R /xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
  --consensusDeterminationModel USE_READS \
  -I /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/${IND}.sorted_RGadded_dupmarked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/sulidae/indelmaps/${IND}.intervals \
  -o /xdisk/mcnew/dannyjackson/sulidae/indelrealignment/${IND}.final.bam



for job in {11127476..11127509}; do
    scancel $job
done

# lingering individual

sbatch --account=mcnew \
--job-name=indelrealign_SRR19149562 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.indelrealign2_SRR19149562.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/sulidae/indelrealign.sh SRR19149562

# Submitted batch job 12031275