#!/usr/bin/env bash
#SBATCH --job-name=BRBO4taxa_dsuite
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBO4taxa_dsuite.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch 1a_BRBO4taxa_dsuite.sh

# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO_4taxa_likelihoods

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t BRBO_4taxa.nwk --ABBAclustering -g
~/programs/Dsuite/Build/Dsuite Fbranch --Pb-matrix BRBO_4taxa.nwk SETS_tree.txt > BRBO_4taxa_Fbranch.txt
