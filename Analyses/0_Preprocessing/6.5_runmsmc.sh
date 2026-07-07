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










