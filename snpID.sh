#!/usr/bin/env bash
#SBATCH --job-name=snpID
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=5:00:00
#SBATCH --output=slurm_output/snpID.out
#SBATCH --mail-type=ALL

OUTBASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen

~/programs/angsd/angsd -GL 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 \
      -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 15 -minQ 30 -minMapQ 30 \
      -minMaf 0.05 \
      -bam /xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt \
      -out $OUTBASE \
      -P 20