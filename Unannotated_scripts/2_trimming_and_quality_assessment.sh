# trimming and quality assessment

## fastqc

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

# sbatch pre_fastqc.sh


## polyg trimming

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


# sbatch polyg_trimming.sh 


## Trimmomatic 


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

# Submitted batch job 3215923


## Post-trimming fastqc 

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

sbatch post_fastqc.sh 
Submitted batch job 3296593



# xdisk/mcnew/dannyjackson/sulidae/raw_sequences/trimming/trimmedseq/_trimmed_2P.fq.gz'





while read -r file;
do
bgzip $file
done <- fastq_filenames.txt

mkdir  /xdisk/mcnew/dannyjackson/sulidae/trimming/


#!/bin/bash

#SBATCH --job-name=fastqc_trimmed
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=fastqc_trimmed.%j

module load fastqc/0.11.9

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/

fastqc -t 12 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/*.fastq.gz

sbatch fastqc_polyg.sh 
Submitted batch job 9469286

while read file;
 do unzip $file 
done < filenames_SRA.txt

grep 'warn' */fastqc_data.txt > warnings.txt
grep 'fail' */fastqc_data.txt > failings.txt






grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_adaptertrimmed_fastas.txt

grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_adaptertrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' failings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_polygtrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' warnings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_polygtrimmed_fastas.txt

grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_finaltrim_fastas.txt

grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_finaltrim_fastas.txt