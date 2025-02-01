# 4.X_snpID 

cd /xdisk/mcnew/dannyjackson/sulidae/referencelists

ls /xdisk/mcnew/dannyjackson/cardinals_dfinch/indelrealignment/*bam > allsamplebams.txt


#!/bin/bash

#SBATCH --job-name=snps
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.snps.%j

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists

~/programs/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/allsamplebams.txt -out /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/allsnps -nThreads 10 