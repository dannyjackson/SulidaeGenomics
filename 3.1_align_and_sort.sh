# align and sort 
#!/bin/bash

#SBATCH --job-name=aligning
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.aligning.%j

module load bwa
module load picard
module load samtools
module load parallel
module list
 
cd /xdisk/mcnew/dannyjackson/sulidae/

echo "Job started on `date`"

booby=(SRR19149590 SRR19149589 SRR19149578 SRR19149567 SRR19149562 SRR19149561 SRR19149560 SRR19149559 SRR19149558 SRR19149557 SRR19149588 SRR19149587 SRR19149586 SRR19149585 SRR19149584 SRR19149583 SRR19149582 SRR19149581 SRR19149580 SRR19149579 SRR19149577 SRR19149576 SRR19149575 SRR19149574 SRR19149573 SRR19149572 SRR19149571 SRR19149570 SRR19149569 SRR19149568 SRR19149566 SRR19149565 SRR19149564 SRR19149563)

# echo "indexing reference">>bwa_alignment_log.txt

# bwa index /xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna

echo "Aligning fastas"
parallel -j 12 -a /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt 'C={}; bwa mem -t 12 /xdisk/mcnew/dannyjackson/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${C}_trimmed_1P.fq.gz \
/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${C}_trimmed_2P.fq.gz | \
samtools view -b -o /xdisk/mcnew/dannyjackson/sulidae/bamfiles/${C}.bam -S' 

# sbatch align.sh 


#!/bin/bash

#SBATCH --job-name=sorting
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sorting.%j

module load bwa
module load picard
module load samtools
module load parallel
module list


echo "Sorting bams"

parallel -j 12 -a /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt 'C={}; samtools sort /xdisk/mcnew/dannyjackson/sulidae/bamfiles/${C}.bam -o /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${C}.sorted.bam' 

for job in {3418931..3418964}; do
    scancel $job
done

for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=sort_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=output.sort_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=12 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/sort.sh $i
done

#!/bin/bash
IND=$1
module load samtools
samtools sort /xdisk/mcnew/dannyjackson/sulidae/bamfiles/${IND}.bam -o /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}.sorted.bam

for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=RG_MD_index_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=output.RG_MD_index_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=12 \
	--time=50:00:00 \
	/xdisk/mcnew/dannyjackson/sulidae/RG_MD_index.sh $i
done

#!/bin/bash
IND=$1
module load picard
module load samtools
module load parallel

mkdir /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}

picard AddOrReplaceReadGroups \
I=/xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}.sorted.bam \
O=/xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.sorted_RGadded.bam \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${IND}

picard MarkDuplicates \
I=/xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.sorted_RGadded.bam \
O=/xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.sorted_RGadded_dupmarked.bam \
M=/xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.duplicate.metrics.txt

samtools index \
/xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/${IND}/${IND}.sorted_RGadded_dupmarked.bam

# sbatch sort.sh 
# 11124640 - 11124673



# rename sorted bams using species codes and sample numbers
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149590.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/RFBO101.sorted.bam

mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149589.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/RFBO102.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149578.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/RFBO103.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149567.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/RFBO104.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149562.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/RFBO105.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149561.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/RFBO106.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149560.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BRBO201.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149559.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BRBO202.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149558.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BRBO203.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149557.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BRBO204.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149588.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BRBO205.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149587.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/MABO301.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149586.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/MABO302.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149585.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/MABO304.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149584.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/MABO305.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149583.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/MABO306.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149582.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/NABO401.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149581.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/NABO402.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149580.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/NABO403.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149579.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/NABO404.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149577.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/NABO405.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149576.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/NABO406.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149575.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BFBO501.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149574.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BFBO502.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149573.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BFBO503.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149572.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BFBO504.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149571.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BFBO505.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149570.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/BFBO506.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149569.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/PEBO601.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149568.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/PEBO602.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149566.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/PEBO603.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149565.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/PEBO604.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149564.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/PEBO605.sorted.bam
mv /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/SRR19149563.sorted.bam /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/PEBO606.sorted.bam



ls /xdisk/mcnew/dannyjackson/sulidae/sortedbamfiles/*bam | sed 's/\.sorted\.bam//g' | awk 'BEGIN {FS = "/"} {print $7}' > /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt