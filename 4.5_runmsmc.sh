4.5_runmsmc.sh
# for a single individual
sbatch msmc_3_runMSMC_slurm.sh
# or without using slurm
./msmc_3_runMSMC.sh
~/programs/msmc2_scripts/msmc_3_runMSMC.sh 


sbatch --account=mcnew \
--job-name=msmcinput_BFBO501 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.msmcinput_BFBO501.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=5:00:00 \
~/programs/msmc2_scripts/msmc_3_runMSMC.sh 

Submitted batch job 3690843


for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=msmc_run.${i} \
        --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/msmc_run.${i}.%j \
	--nodes=1 \
	--ntasks-per-node=1 \
	--time=10:00:00 \
	~/programs/msmc2_scripts/msmc_3_runMSMC.sh $i
done

3690918 - 3690951

# RFBO103 and RFBO105 failed due to mem
sbatch --account=mcnew \
--job-name=msmc_run.RFBO103 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/msmc_run.RFBO103.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=10:00:00 \
--mem=100G \
~/programs/msmc2_scripts/msmc_3_runMSMC.sh RFBO103
#Submitted batch job 11164669

sbatch --account=mcnew \
--job-name=msmc_run.RFBO105 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/msmc_run.RFBO105.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=10:00:00 \
--mem=100G \
~/programs/msmc2_scripts/msmc_3_runMSMC.sh RFBO105
# Submitted batch job 11164670


# first, remove scaffolds and sex chromosomes from input directory
mkdir input_scaffolds
mv input/*NW_* input_scaffolds/
mkdir input_sexchroms
mv input/*NC_087547* input_sexchroms/ # W
mv input/*NC_087548* input_sexchroms/ # Z


#!/bin/sh

# to run this script:
# ./msmc_3_runMSMC.sh

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=~/programs/msmc2_scripts
source ${scriptdir}/msmc_params.sh

POP_OR_IND=$1
echo $1

if [ $NR_IND == 1 ]; then

        find ${OUTDIR}/input/msmc_input.${POP_OR_IND}.*.txt -size 0 -delete
        ls ${OUTDIR}/input/msmc_input.${POP_OR_IND}.*.txt > ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}
else
        for i in `cat ${POP_OR_IND}_IND.txt`
       do echo $i
       IND=$i

       for s in `cat SCAFFOLDS.txt`
               do echo $s
               ls ${OUTDIR}/input/msmc_input.${IND}.${s}.txt >> ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}
      done
  done
fi

module load msmc2


### Report settings/parameters:
date
echo "Script: msmc_3_onepop.sh"
echo "Run name: $RUN_NAME"
echo "SNP calling method: $METHOD"
echo "Period setting: $P_PAR"
echo "Nr of individuals (1 or 2+): $NR_IND"
echo "Population or individuals ID: $POP_OR_IND"
echo "Individual: "
echo "Scaffolds: SCAFS_INPUT_${POP_OR_IND}"
echo "Iterations: 100"

if [ $NR_IND == 1 ]
        then
        echo "Running MSMC for one individual"
        MSMC_INPUT=`cat ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${OUTDIR}/output/msmc_output.${RUN_NAME}

        if [ -f "${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}" ]
                then
                        echo "MSMC_INPUTS: SCAFS_INPUT_${POP_OR_IND}_noLG9"
                        echo "MSMC_OUTPUT: $MSMC_OUTPUT"
                else
                        echo "MSMC_INPUT does not exist! Exiting now"
                        exit 1
        fi

        ~/programs/msmc_2.0.0_linux64bit -t 16 -p $P_PAR -i 100 -o $MSMC_OUTPUT -I 0,1 $MSMC_INPUT

        mv $MSMC_OUTPUT*loop.txt ${OUTDIR}/output/log_and_loop/
        mv $MSMC_OUTPUT*log ${OUTDIR}/output/log_and_loop/
else
        echo "Running MSMC for $NR_IND individuals"
        MSMC_INPUT=`cat ${OUTDIR}/input/SCAFS_INPUT_${POP_OR_IND}`
        MSMC_OUTPUT=${OUTDIR}/output/msmc_output.${POP_OR_IND}.${RUN_NAME}
        n=$(expr ${NR_IND} - 2)
        INDEX=$(for num in `seq 0 ${n}`; do echo -n "${num},"; done; echo ${NR_IND})
        
        ~/programs/msmc_2.0.0_linux64bit -t 16 -p $P_PAR -i 100 -o ${MSMC_OUTPUT} -I `echo $INDEX` $MSMC_INPUT

        mv $MSMC_OUTPUT*loop.txt ${OUTDIR}/output/log_and_loop/
        mv $MSMC_OUTPUT*log ${OUTDIR}/output/log_and_loop/
fi


echo "done running msmc2 for ${POP_OR_IND}"
date



# P_PAR=1*2+25*1+1*2+1*3 

# now make plots in R
for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt`;
	do echo $i
	IND=$i
	Rscript msmc_plot.r /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc $i
done


scp -r dannyjackson@filexfer.hpc.arizona.edu:/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots /Users/danjack/Documents/SulidaePaper/Summer2024/RawFigures




11164669




# The -p 1*2+15*1+1*2 option defines the time segment patterning. By default, MSMC uses 32 time segments, grouped as 1*2+25*1+1*2+1*3, which means that the first 2 segments are joined (forcing the coalescence rate to be the same in both segments), then 25 segments each with their own rate, and then again two groups of 2 and 3, respectively. MSMC2 run time and memory usage scales quadratically with the number of time segments. Here, since we are only analysing a single chromosome, you should reduce the number of segments to avoid overfitting. That's why I set 18 segments, with two groups in the front and back. Grouping helps avoiding overfitting, as it reduces the number of free parameters.












# MSMC on Multiple Genomes 
# step 2 -- generate input
grep 'BFBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO_samplecodes.txt
grep 'PEBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt
grep 'MABO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_samplecodes.txt
grep 'NABO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO_samplecodes.txt
grep 'BRBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO_samplecodes.txt
grep 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_samplecodes.txt

echo -e 'BFBO\nPEBO\nMABO\nNABO\nBRBO\nRFBO' > /xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt

To combine multiple individuals, the steps generally follow what is detailed above. However, in step 2, you will instead use the script msmc_2_generateInput_multiInd.sh, which will combine all of the individuals from a given population. For this script, which is submitted using submit_2_multi.txt, you will need a file with the names of the individuals you want combined in your MSMC run (named POP_IND.txt, where "POP" is your population name).

nano ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh
nano ~/programs/msmc2_scripts/submit_2_multi.txt

for POP in `cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/speciescodes.txt`;
	do echo $POP
	sbatch --account=mcnew \
	--job-name=msmc_run.${POP} \
        --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/msmc_run.${POP}.%j \
	--nodes=1 \
	--ntasks-per-node=1 \
	--time=10:00:00 \
                ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh `echo ${POP}_IND.txt` `echo $POP`
done

sbatch --account=mcnew \
--job-name=msmc_run.BFBO \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/msmc_run.BFBO.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=10:00:00 \
        ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh `echo BFBO_IND.txt` `echo BFBO`
# 3690978 - 3690983

chmod +x ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh 
~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh `echo BFBO_IND.txt` `echo BFBO`

cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/BFBO_samplecodes.txt ./BFBO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt ./PEBO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/MABO_samplecodes.txt ./MABO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/NABO_samplecodes.txt ./NABO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/BRBO_samplecodes.txt ./BRBO_IND.txt
cp /xdisk/mcnew/dannyjackson/sulidae/referencelists/RFBO_samplecodes.txt ./RFBO_IND.txt

ls vcf/ > vcf_filenames.txt

for job in {3690995..3691382}; do
    scancel $job
done


#!/bin/bash
module load parallel
parallel -j 12 -a vcf_filenames.txt 'C={}; bgzip -c vcf/$file > gzvcf/$file.gz'

sbatch --account=mcnew \
--job-name=bgzipvcfs \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/bgzipvcfs.%j \
--nodes=1 \
--ntasks-per-node=12 \
--time=6:00:00 \
./bgzipvcfs.sh

# Submitted batch job 3691397

# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask/ind_mask.BFBO501.NC_087518.1.bed.gz:


# Do I really need a gzipped vcf?
less ~/programs/msmc2_scripts/msmc_2_generateInput_multiInd.sh
 #${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
generate_multihetsep.py mask/BFBO.mask_file.NC_087513.1 --mask=mask/GCF_963921805_NC_087513.1.mask.150.50.bed.gz vcf/BFBO.vcf_file.NC_087513.1 > input/BFBO_testinput.txt # wi
thout repeat mask
MASK_GENOME= mask/GCF_963921805_${SCAFFOLD}.mask.${k}.50.bed.gz

In step 3, you will then need to change -I, which specifies which haplotypes to analyze. For one individual, this is 0,1, for two individuals 0,1,2,3, and so on. This will be done automatically in the “multiInd" script, which can simply be run as follows:

# for multiple individuals in a population
sbatch msmc_3_runMSMC_slurm.sh
# or without using slurm
~/programs/msmc2_scripts/msmc_3_runMSMC.sh









#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=~/programs/msmc2_scripts
source ${scriptdir}/msmc_params.sh

### Variables:
IND=`cat $1`
POP=$2

for s in `cat SCAFFOLDS.txt`
        do SCAFFOLD=$s

        MSMC_INPUT=${OUTDIR}/input/msmc_input.${POP}.${SCAFFOLD}.txt

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_multiInd"
        echo "Individuals: ${IND}"
        echo "Population: $POP"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"

        for ind in $IND
                do INDMASK=`ls ${OUTDIR}/mask/ind_mask.${ind}.${SCAFFOLD}.bed.gz`
                echo "--mask=$INDMASK " >> ${OUTDIR}/mask/${POP}.mask_file.$SCAFFOLD
                INDVCF=`ls ${OUTDIR}/vcf/${ind}.${SCAFFOLD}.msmc.vcf.gz`
                echo $INDVCF >> ${OUTDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}
        done

### Generate MSMC input files:
        if [ $METHOD == samtools ]
                then
                MASK_GENOME=${OUTDIR}/mask/${prefix}_${SCAFFOLD}.mask.${k}.50.bed.gz

                echo "MAPPABILITY MASK: ${MASK_GENOME}"
                echo "Creating MSMC input file WITH individual mask (samtools)"
                #${MSMCTOOLS}/generate_multihetsep.py --negative_mask=$MASK_REPEATS --mask=$MASK_INDIV $VCF > $MSMC_INPUT # with repeat mask
                generate_multihetsep.py `cat ${OUTDIR}/mask/${POP}.mask_file.${SCAFFOLD}` --mask=$MASK_GENOME `cat ${OUTDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}` > ${MSMC_INPUT} # without repeat mask


        elif [ $METHOD == gatk ]
                then
                echo "Creating MSMC input file WITHOUT individual mask (gatk)"
                MASK_GENOME=`ls ${OUTDIR}/mask/${prefix}_${SCAFFOLD}.mask.${k}.50.bed.gz`
                #msmc-tools/generate_multihetsep.py --negative_mask=$MASK_REPEATS $VCF > $MSMC_INPUT # with repeat mask
                generate_multihetsep.py --mask=$MASK_GENOME `cat ${OUTDIR}/vcf/${POP}.vcf_file.${SCAFFOLD}` > $MSMC_INPUT # without repeat mask

                # NOTE THAT THIS WAS CHANGED 10 FEB 2024
                # AND HAS NOT YET BEEN TESTED TO MAKE SURE IT WORKS
                                
                echo "Creating individual mask. Note that your input VCF should include ALL sites (variant & invariant)."
                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz
                                
                VCF_OUT=${VCF}.parsed.vcf
                vcfAllSiteParser.py $SCAFFOLD $MASK_INDIV $VCF_OUT
                                
                echo "Creating MSMC input file with new individual mask"
                                
                generate_multihetsep.py --mask=$MASK_INDIV --mask=$MASK_GENOME $VCF > $MSMC_INPUT # with new repeat mask
        fi

done

echo "Done with script."
date

####