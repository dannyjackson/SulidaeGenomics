#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.trimming.%j

module load parallel/20241122
module load trimmomatic/0.39

cd /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/

echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/sulidae/trimming/trimmedseq//trim_log.txt

parallel -j 12 -a /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt 'C={}; trimmomatic PE -threads 12 \
/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/polyg/${C}_1_polyg.fastq.gz  /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/polyg/${C}_2_polyg.fastq.gz  \
-baseout /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/${C}_trimmed.fq.gz \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/trim_log.txt' 
