# 4.4 Generate input for msmc
Now that you have vcf files and mask files for each individual and each chromosome/scaffold, you can use these to generate the input for MSMC2 using msmc_2_generateInput_singleInd.sh, which is run using the submit_2.txt script. This script is for analyzing a single genome/individual separately; if you have multiple individuals per population/species that you want to combine in analyses, then follow the instructions in "MSMC on Multiple Individuals" below. NOTE: I would recommend doing a first-pass analysis with all individuals separate to make sure that they’re all reasonably similar, and then combining multiple individuals in a later run for greater statistical power.

The workhorse of this script is the line

generate_multihetsep.py --mask=${MASK_INDIV} --mask=${MASK_GENOME} $VCF > $MSMC_INPUT
BFBO501

~/programs/msmc-tools/generate_multihetsep.py --mask=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask/ind_mask.BFBO501.NC_087532.1.samtools.bed.gz --mask=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask/GCF_963921805_NC_087532.1.mask.150.50.bed.gz /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/BFBO501.NC_087532.1..samtools.vcf.gz > /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/input/BFBO501.NC_087532.1.input.txt

The output of this script, the MSMC input, has a list of heterozygous sites in the genome, and the distance between these heterozygous sites.

# check if any files in the input directory are empty
find input/* -size 0 -delete
# 
#!/bin/bash

# all parameters come from the msmc_param control file
# make edits there before using this script!
scriptdir=$(dirname "$0")
source ${scriptdir}/msmc_params.sh

### Variables:
IND=$1

module list

for s in `cat SCAFFOLDS.txt`
        do echo "working on scaffold $s"
        SCAFFOLD=$s
        VCF=`ls ${OUTDIR}/vcf/${IND}.${SCAFFOLD}.${METHOD}.vcf.gz`

        #MASK_REPEATS=repeats.bed.gz # Needs to be gzipped
        MSMC_INPUT=${OUTDIR}/input/msmc_input.${IND}.${SCAFFOLD}.txt

        printf "\n \n \n \n"
        date
        echo "Script: msmc_2_generateInput_singleInd"
        echo "Individual: ${IND}"
        echo "Scaffold: ${SCAFFOLD}"
        echo "Phasing: ${PHASING}"
        echo "Method: ${METHOD}"
        echo "MSMC input file: ${MSMC_INPUT}"
        echo "VCF: ${VCF}"

        if [ -f "$VCF" ]
                then
                        echo "VCF exists, starting creation of input for MSMC2!"

### Generate MSMC input files:
                if [ $METHOD == samtools ]
                        then
                                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz # store indiv.mask file path
                                MASK_GENOME=${OUTDIR}/mask/prefix_${SCAFFOLD}.mask.${k}.50.bed.gz
                                echo "MASK: ${MASK_INDIV}"
                                echo "MAPPABILITY MASK: ${MASK_GENOME}"
                                echo "Creating MSMC input file WITH individual mask (samtools)"

                                generate_multihetsep.py --mask=${MASK_INDIV} --mask=${MASK_GENOME} $VCF > $MSMC_INPUT # without repeat mask

                elif [ $METHOD != samtools ]
                        then
                                # NOTE THAT THIS WAS CHANGED 10 FEB 2024
                                # AND HAS NOT YET BEEN TESTED TO MAKE SURE IT WORKS
                                
                                echo "Creating individual mask. Note that your input VCF should include ALL sites (variant & invariant)."
                                MASK_INDIV=${OUTDIR}/mask/ind_mask.${IND}.${SCAFFOLD}.${METHOD}.bed.gz
                                MASK_GENOME=${OUTDIR}/mask/prefix_${SCAFFOLD}.mask.${k}.50.bed.gz
                                
                                VCF_OUT=${VCF}.parsed.vcf
                                vcfAllSiteParser.py $SCAFFOLD $MASK_INDIV $VCF_OUT
                                
                                echo "Creating MSMC input file with new individual mask"
                                
                                generate_multihetsep.py --mask=$MASK_INDIV --mask=$MASK_GENOME $VCF > $MSMC_INPUT # with new repeat mask
                fi

        else
                echo "VCF does not exist, moving on to next scaffold"
        fi

        echo "Done with ${SCAFFOLD}; moving on to next scaffold"
done;

echo "Done with script."
date

####


# test script
submit_2.txt

sbatch --account=mcnew \
--job-name=msmcinput_BFBO501 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.msmcinput_BFBO501.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=5:00:00 \
~/programs/msmc2_scripts/msmc_2_generateInput_singleInd.sh BFBO501

Submitted batch job 3690788


~/programs/msmc-tools/generate_multihetsep.py --mask=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask/ind_mask.BFBO501.NC_087532.1.samtools.bed.gz --mask=/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/mask/GCF_963921805_NC_087532.1.mask.150.50.bed.gz  /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/BFBO501.NC_087532.1..samtools.vcf.gz > /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/input/msmc_input.BFBO501.NC_087532.1.txt # without repeat mask

# remove missing genotypes didn't work -- the vcf doesn't have any missing files

grep -n 'self.end' generate_multihetsep.py

adding mask: GCF_963921805_mask.150.50.fa
Traceback (most recent call last):
  File "/home/u15/dannyjackson/programs/msmc-tools/generate_multihetsep.py", line 207, in <module>
    maskIterators.append(MaskIterator(f))
                         ^^^^^^^^^^^^^^^
  File "/home/u15/dannyjackson/programs/msmc-tools/generate_multihetsep.py", line 19, in __init__
    self.readLine()
  File "/home/u15/dannyjackson/programs/msmc-tools/generate_multihetsep.py", line 30, in readLine
    self.end = int(fields[2])
               ^^^^^^^^^^^^^^
ValueError: invalid literal for int() with base 10: '0.500'


for i in `cat /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_samplecodes.txt`;
	do echo $i
	IND=$i
	sbatch --account=mcnew \
	--job-name=msmcinput.${i} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/msmcinput.${i}.%j \
	--nodes=1 \
	--ntasks-per-node=4 \
	--time=25:00:00 \
	~/programs/msmc2_scripts/msmc_2_generateInput_singleInd.sh $i
done

3690791 - 3690824