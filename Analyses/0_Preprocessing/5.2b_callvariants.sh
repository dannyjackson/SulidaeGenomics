#!/usr/bin/env bash
#SBATCH --job-name=callvariants
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_output/callvariants.%A_%a.out
#SBATCH --mail-type=ALL

set -euo pipefail


# all parameters come from the msmc_param control file
# make edits there before using this script!

## VARIABLES:
LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

IND=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')

GENOME=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks
prefix=${OUTDIR}/mask/GCA_031468815_ # prefix of genome masks
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams # directory with bamfiles
BAMFILE=${BAMDIR}/${IND}.final.bam
DEPTHFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/${IND}_depthstats.txt 
DEPTHOUTFILE=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/mean_by_scaffold/${IND}_depthstats.txt 
METHOD="samtools" # method for variant calling (e.g., samtools)
MSMCTOOLS=~/programs/msmc-tools/ # path to msmc-tools

module load samtools
module load bcftools
module load python/3.7.4


printf "\n \n \n \n"
date
echo "Individual: $IND"
echo "Bamfile: $BAMFILE"

if [ -f "${BAMFILE}.bai" ]
        then
                echo "Bamfile index already exists, moving on!"
        else
                echo "Bamfile index does not exist, creating index"
                samtools index $BAMFILE ${BAMFILE}.bai
fi


for s in `cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt`; \
        do echo "Handling scaffold $s"; \

        ### Calculate mean coverage (to be used as input for bamCaller.py):
        MEANCOV=`grep $s $DEPTHFILE | awk '{sum += $3} END {if (NR==0) print NR; else print sum / NR}' | tr ',' '.'` # calculate mean coverage
        echo ${IND} ${s} $MEANCOV >> $DEPTHOUTFILE # save mean coverage in separate file
        echo "Mean coverage for this individual, scaffold ${s}: $MEANCOV"

        ### Generate a single-sample VCF and a mask-file:
        MASK_IND=${OUTDIR}/mask/ind_mask.${IND}.${s}.bed.gz # Individual mask file to be created
        VCF=${OUTDIR}/vcf/${IND}.${s}.${METHOD}.vcf # VCF file to be created

        #If genome isn't indexed, add:
        #samtools faidx $GENOME
        if [ "$METHOD" == "samtools" ]; then
                echo "starting samtools alignment"
                bcftools mpileup -Ou -r ${s} --threads 16 -f $GENOME $BAMFILE \
                | bcftools call -c --threads 16 -V indels \
                | $MSMCTOOLS/bamCaller.py $MEANCOV $MASK_IND > ${VCF}

                ## Only DP > 9:
                #samtools mpileup -q 20 -Q 20 -C 50 -u -r $SCAFFOLD -f $GENOME $BAMFILE | bcftools call -c -V indels | bcftools view -i 'INFO/DP>9' | $MSMCTOOLS/bamCaller.py $MEANCOV $MASK_IND | gzip -c > $VCF.gz

                # -q = min. mapping qual; -Q = min. base qual; -C = coefficient for downgrading mapping qual for reads w/ excessive mismatches; -u = generate uncompressed VCF/BCF; -r = only specified region; -f = fasta.
                # bcftools: "call" will call SNPs/indels; "-V indels" skips indels; -c = consensus caller.
        fi
        echo "done with scaffold ${s}; moving on to next scaffold"
done

### Report:
echo "Filtered VCF and mask created for all scaffolds for ${IND}."
echo "Done with script."
date

# sbatch --array=1-29 4.3_callvariants.sh


