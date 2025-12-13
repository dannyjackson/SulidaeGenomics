#!/bin/bash

#SBATCH --job-name=post_fastqc
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.post_fastqc.%j

module load fastqc/0.11.9
mkdir  /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/post_fastqcs/
cd /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/post_fastqcs/

while read -r ID;
do
fastqc -t 12 /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/"$ID"_trimmed_1P.fq.gz

fastqc -t 12 /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/"$ID"_trimmed_2P.fq.gz

done < /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt

while read -r file;
do
bgzip $file
done <- fastq_filenames.txt

mkdir  /xdisk/mcnew/dannyjackson/sulidae/trimming/