#!/usr/bin/env bash
#SBATCH --job-name=dsuite_rfbo
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/dsuite_rfbo%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/RFBO

VCF=RFBO_BRBO.autosomes.vcf.gz

# ~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt --ABBAclustering -g
~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t rfbo.nwk --ABBAclustering -g -n RFBOtree