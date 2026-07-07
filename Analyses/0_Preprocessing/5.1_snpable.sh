
datasets download genome accession GCF_963921805.1 --include gff3,rna,cds,protein,genome,seq-report

# msmc

cd ~/programs
git clone https://github.com/stschiff/msmc2
git clone https://github.com/stschiff/msmc-tools
git clone https://github.com/jessicarick/msmc2_scripts

cd ~/programs/msmc2_scripts/
# edit msmc_params.sh
pip3 install --user whatshap

# following Jessi Rick's tutorial, found here: https://github.com/jessicarick/msmc2_scripts?tab=readme-ov-file
# Step 0: Create mappability mask

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc

sbatch --account=mcnew \
--job-name=run_snpable2 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.run_snpable2.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
~/programs/msmc2_scripts/run_snpable2.sh

# Submitted batch job 12032071

awk '{print $1}' /xdisk/mcnew/dannyjackson/sulidae/old/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna.fai > ~/programs/msmc2_scripts/SCAFFOLDS.txt

# then edit the makeMappabilityMask.py script (lines 26 & 30) for your genome, and run the script
# on line 30, use curly braces {} to indicate where in the name the scaffold name should go

# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc

splitfa $GENOME $k | split -l 20000000
cat snpable/x* >> GCF_963921805__split.150

gzip -dc xx??.sam.gz | gen_raw_mask.pl > rawMask_35.fa  

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 8 -R 1000000 -O 3 -E 3 ${GENOME} ${prefix}_split.${k} > ${prefix}_split.${k}.sai
bwa samse -f ${prefix}_split.${k}.sam $GENOME ${prefix}_split.${k}.sai ${prefix}_split.${k}


sbatch --account=mcnew \
--job-name=snpable_2 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.snpable_2.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
snpable_2.sh

Submitted batch job 3526455

11143339

sbatch --account=mcnew \
--job-name=snpable_1 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.snpable_1.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=48:00:00 \
snpable_1.sh

Submitted batch job 12044740 # only successfully ran the splitfa command because my module load bwa line didn't work
Submitted batch job 12044756 # all but the splitfa command

#!/bin/bash

module load bwa

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1

echo "status: indexing reference genome"

bwa index /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc

echo "status: splitting genome"

# ~/programs/seqbility-20091110/splitfa /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna 150 | split -l 20000000

cat snpable/x* >> GCF_963921805_split.150

# gzip -dc xx??.sam.gz | gen_raw_mask.pl > rawMask_35.fa  

echo "status: performing bwa aln"

bwa aln -t 8 -R 1000000 -O 3 -E 3 /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna GCF_963921805_split.150 > GCF_963921805_split.150.sai

echo "status: performing bwa samse"

bwa samse /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna GCF_963921805_split.150.sai GCF_963921805_split.150 > GCF_963921805_split.150.sam

mv GCF_963921805_split.150 GCF_963921805_split.150.tmp

#!/bin/bash
module load bwa/0.7.17
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module spider samtools/1.19.2

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc

sed 's/>//g' GCF_963921805_split.150.tmp > GCF_963921805_split.150

echo "status: starting to generate rawMask"
~/programs/seqbility-20091110/gen_raw_mask.pl GCF_963921805_split.150.sam >  GCF_963921805_rawMask.150.fa
# gen_raw_mask.pl ${prefix}_split.${k}.sam > ${prefix}_rawMask.${k}.fa

echo "status: raw mask created as  GCF_963921805__rawMask.35.fa, now generating final mask with stringency r=50%"
~/programs/seqbility-20091110/gen_mask -l 150 -r 0.5  GCF_963921805_rawMask.150.fa >  GCF_963921805_mask.150.50.fa

echo "status: all done! final mask saved as GCF_963921805_mask.150.50.fa"

sbatch --account=mcnew \
--job-name=snpable_3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.snpable_3.%j \
--nodes=1 \
--ntasks-per-node=12 \
--time=24:00:00 \
snpable_3.sh

Submitted batch job 11162246 # elgato

# Next, you'll have to turn this into a "mappability mask" to be used with MSMC. To do this, edit the makeMappabilityMask.py script. You'll need to edit the paths indicated on lines 26 and 30. On line 26, direct the script to the masked fasta file created by the run_snpable2 script (note: I recommend using the full path).
# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/GCF_963921805__mask.150.50.fa
# On line 30, you need to specify the output location for the individual scaffold masks. In the path, use curly braces {} to indicate where the name of the scaffold should go. For example, my line 30 reads as follows:

#/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask/GCF_963921805_{}.mask.150.50.bed.gz
# Once those two lines are edited, run the script. Note: this script requires python2 (NOT python3) to run.
sed 's/>>/>/g' /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/GCF_963921805_mask.150.50.fa > /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/GCF_963921805_revised_mask.150.50.fa

#!/bin/bash
module load python/2.7.15
python2 ~/programs/msmc-tools/makeMappabilityMask.py

sbatch --account=mcnew \
--job-name=makeMappabilityMask \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.makeMappabilityMask.%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=24:00:00 \
makeMappabilityMask.sh

Submitted batch job 11162993



# Done with snppable!
Step 1: Calling variants
# check if any files in the mask directory are empty
find mask/* -size 0 -delete
# edit msmc_params
nano ~/programs/msmc-tools/msmc_params 
source ~/programs/msmc2_scripts/submit_1.txt


for i in `cat INDS.txt`;
        do echo $i
        IND=$i
        sbatch --account=latesgenomics \
        --job-name=msmc_${i} \
        --mail-type=ALL \
        --output=outs/stdout_${i}_msmc1 \
        --error=outs/stderr_${i}_msmc1 \
        --nodes=1 \
        --ntasks-per-node=16 \
        --time=1-00:00:00 \
        $scriptdir/msmc_1_call.sh $i
done

awk '{print $1}' /xdisk/mcnew/dannyjackson/sulidae/old/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna.fai > SCAFFOLDS.txt


for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=msmc_1_call_${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.msmc_1_call_${i}.%j \
	--nodes=1 \
	--ntasks-per-node=16 \
	--time=25:00:00 \
	~/programs/msmc2_scripts/msmc_1_call.sh $i
done

head /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt

~/programs/msmc2_scripts/msmc_1_call.sh BRBO201

# This created masks with '>' in the filename before each chromosome. Use this script to rename them without that character
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask
ls > masknames.txt

while IFS= read -r file; do
    new_name="${file//>/}"  # Remove all '>' characters from the filename
    if [[ "$file" != "$new_name" ]]; then
        echo "$file" "$new_name"
    fi
done < masknames.txt
rm masknames.txtnode


