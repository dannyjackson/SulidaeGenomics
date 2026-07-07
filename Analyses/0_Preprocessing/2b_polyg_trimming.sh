#!/bin/bash

#SBATCH --job-name=polygtrimming
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.polyg.%j

while read -r ID;
do

echo "Beginning polyg trimming for "$ID>>/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/polyg/polyg_trim_log.txt

  ~/programs/fastp --trim_poly_g -Q -L -A -w 24 --poly_g_min_len 10 -i /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/"$ID"_1.fastq.gz -o /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/polyg/"$ID"_1_polyg.fastq.gz -I /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/"$ID"_2.fastq.gz -O /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/polyg/"$ID"_2_polyg.fq.gz

done < /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt
