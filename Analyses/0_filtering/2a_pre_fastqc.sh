#!/bin/bash

#SBATCH --job-name=fastqc_pre
#SBATCH --ntasks=8
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.fastqc_pre.%j

module load fastqc/0.11.9
module load samtools/1.19.2
module load htslib/1.19.1
module load parallel/20241122

cd /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/
ls *q > fastq_filenames.txt

parallel -j 8 -a fastq_filenames.txt 'C={}; bgzip ${C}' 

mkdir  /xdisk/mcnew/dannyjackson/sulidae/trimming/
cd /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/
mkdir /xdisk/mcnew/dannyjackson/sulidae/trimming/pre_fastqcs/

cd /xdisk/mcnew/dannyjackson/sulidae/trimming/pre_fastqcs/

fastqc -t 8 /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/*fastq.gz
