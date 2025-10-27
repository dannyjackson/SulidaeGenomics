#!/bin/bash
#SBATCH --job-name=snpable_mapabilitymask
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1          # threads per sample
#SBATCH --mem=5G                     # ~5G/core is plenty for bwa-mem2
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/snpable_mapabilitymask.%A_%a.out


module load python/2.7.15

python2 /xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/makeMappabilityMask.py

echo "all done! final mask saved as ${prefix}_mask.${k}.50.fa"

date