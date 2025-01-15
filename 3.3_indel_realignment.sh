# indel realignment
# Make dictionary for reference genome, necessary for doing indel realignment

cd ..
mkdir indelmaps
cd indelmaps 


java -jar picard.jar Create

~/programs/jdk1.8.0_411/bin/java -jar ~/programs/picard.jar CreateSequenceDictionary -R GCF_963921805.1_bPhaCar2.1_genomic.fna -O GCF_963921805.1_bPhaCar2.1_genomic.dict


# Make indel maps
# make indel maps slurm array

for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=sort_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=output.maps_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=48:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/indelmap.sh $i
done

#!/bin/bash
IND=$1

module load samtools
module load parallel

apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
-I /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/${IND}_sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/sulidae/indelmaps/${IND}.intervals


sbatch indelmap.sh 
Submitted batch job 11118922

# test
apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
-I /xdisk/mcnew/dannyjackson/sulidae/clipoverlap/BFBO501_sorted.marked.clipped.bam \
-o /xdisk/mcnew/dannyjackson/sulidae/indelmaps/BFBO501.intervals

# Realign indels using indel maps
#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j


mkdir /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/indelrealignment


parallel -j 12 -a /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt 'C={}; ~/programs/jdk1.8.0_411/bin/java -jar ~/programs/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
  --consensusDeterminationModel USE_READS \
  -I /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/clipoverlap/${C}_sorted.marked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/indelmaps/${C}.intervals \
  -o /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/indelrealignment/${C}.realigned.bam' 




# didn't do these, not sure what's up with them yet

# some alignment stats 
echo "Computing stats on sorted clipped and indel realigned bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt 

echo "Computing stats on unsorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/alignmentstats.txt 

echo "Job ended on `date`"


sbatch align_finaltrim.round1.sh 
Submitted batch job 9481231




#!/bin/bash

#SBATCH --job-name=align_sequencess
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.align.%j

module load bowtie2
module load picard
module load samtools
module list
 
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

echo "Job started on `date`"

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507)

echo "Aligning fastas"


for u in "${finch[@]}"; do
echo ${u}
bowtie2 -p 24 --phred33 --very-sensitive-local -x /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -I 149 -X 900 --rg-id ${u} --rg SM:${u} -1 ${u}_trimmed2_1P.fq.gz -2 ${u}_trimmed2_2P.fq.gz -U ${u}_trimmed2_1U.fq.gz, ${u}_trimmed2_2U.fq.gz -S /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam 2>> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/alignment.output.txt
done

# convert to bams and sort
echo "Converting sams to bams"

for u in "${finch[@]}"; do 
samtools view -S -b /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam

echo "Sorting bams"

samtools sort /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam 
done


# some alignment stats 
echo "Computing stats on sorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt 

echo "Computing stats on unsorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/alignmentstats.txt 

echo "Job ended on `date`"

sbatch align_finaltrim.sh 
Submitted batch job 9483129


#!/bin/bash

#SBATCH --job-name=align_sequencess
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.align.%j

module load bowtie2
module load picard
module load samtools
module list

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507)


echo "Computing stats on sorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt 

sbatch sortedbam_stats.sh 
Submitted batch job 9515763



ls /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/

ls *all*bam | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt
