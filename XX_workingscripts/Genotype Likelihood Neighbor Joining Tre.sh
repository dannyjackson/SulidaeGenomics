Genotype Likelihood Neighbor Joining Tree

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/NJ_tree

ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/indelrealignment/*bam > /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams.txt

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams.txt -GL 1 -anc /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -out sula -setMinDepthInd 4 -nThreads 8 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6



# filter sample bams to just individuals kept
for line in $(cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt);
do
grep $line /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt
done 

#!/bin/bash

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/NJ_tree

~/programs/angsd/angsd \
  -b /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt \
  -GL 2 \
  -anc /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
  -out sula_filtered_mind24 \
  -doGlf 2 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -minInd 24 \
  -minMapQ 30 \
  -minQ 20 \
  -SNP_pval 1e-6 \
  -nThreads 32

sbatch --account=mcnew \
    --job-name=makebeagle \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.makebeagle.%j \
    --nodes=1 \
    --ntasks-per-node=32 \
    --time=5:00:00 \
    makebeagle.sh

#!/bin/bash

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/NJ_tree

~/programs/angsd/angsd \
  -b /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt \
  -GL 2 \
  -anc /xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna \
  -out sula_filtered_minmaf \
  -doGlf 2 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -minInd 4 \
  -minMapQ 30 \
  -minQ 20 \
  -SNP_pval 1e-6 \
  -nThreads 32 \
  -minMaf 0.05


sbatch --account=mcnew \
    --job-name=makebeagle_minmaf \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.makebeagle_minmaf.%j \
    --nodes=1 \
    --ntasks-per-node=32 \
    --time=5:00:00 \
    makebeagle_minmaf.sh

# rename individuals in beagle file

# 1. Extract real names
sample_names=($(basename -s .final.bam -a $(cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt)))

# 2. Replace header line
zcat sula_filtered.beagle.gz | awk -v OFS="\t" -v names="${sample_names[*]}" '
  BEGIN {
    split(names, n)
  }
  NR==1 {
    printf "%s\t%s\t%s", $1, $2, $3
    for (i=1; i<=length(n); i++) {
      printf "\t%s\t%s\t%s", n[i], n[i], n[i]
    }
    printf "\n"
    next
  }
  { print }
' | gzip > sula_filtered.labeled.beagle.gz


# angsd -bam bamlist.txt -GL 1 -out output -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6
# ngsDist -glf sula.glf.gz -n_ind N -n_sites S -out output.dist
# ngsTools installed with elgato, not puma

#!/bin/bash

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/NJ_tree

# ~/programs/ngsTools/ngsDist/ngsDist -geno sula_filtered.labeled.beagle.gz -n_ind 29 -n_sites 35216786 -out sula_filtered.dist -n_threads 16 
# ~/programs/ngsTools/ngsDist/ngsDist -geno sula_filtered_mind24.beagle.gz -n_ind 29 -n_sites 30054477 -out sula_filtered.mind24.dist -n_threads 16 
~/programs/ngsTools/ngsDist/ngsDist -geno sula_filtered_minmaf.beagle.gz -n_ind 29 -n_sites 27269944 -out sula_filtered.minmaf.dist -n_threads 16 

# -tot_sites 1129395008

sbatch --account=mcnew \
    --job-name=makedist \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.makedist.%j \
    --nodes=1 \
    --ntasks-per-node=16 \
    --time=1:00:00 \
    makedist.sh
sbatch --account=mcnew \
    --job-name=makedist \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.makedist.%j \
    --nodes=1 \
    --ntasks-per-node=16 \
    --time=1:00:00 \
    makedist_mind24.sh
sbatch --account=mcnew \
    --job-name=makedist \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.makedist.%j \
    --nodes=1 \
    --ntasks-per-node=16 \
    --time=1:00:00 \
    makedist_minmaf.sh

# installed fastme on elgato
fastme -i output.dist -o output.tree
# bootstapping
# If you want branch support values on your tree, you can use ngsDist with the option --n_boot_rep and --boot_block_size to bootstrap the input data. ngsDist will output one distance matrix (the first) for the input full dataset, plus --n_boot_rep matrices for each of the bootstrap replicates. After, infer a tree for each of the matrices using the program of your choice and plot them. For example, using FastME on a dataset with 5 bootstrap replicates:


basename -s .final.bam -a $(cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams.txt)

BFBO501
BFBO502
BFBO503
BFBO504
BFBO505
BFBO506
BRBO201
BRBO202
BRBO203
BRBO204
BRBO205
MABO301
MABO302
MABO304
MABO305
MABO306
NABO401
NABO402
NABO403
NABO404
NABO405
NABO406
PEBO601
PEBO602
PEBO603
PEBO604
PEBO605
PEBO606
RFBO101
RFBO102
RFBO103
RFBO104
RFBO105
RFBO106

head -n 2 sula_filtered.dist > sula_filtered.dist.names
paste <(basename -s .final.bam -a $(cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt)) <(grep '^Ind_' sula_filtered.dist | cut -f2-) | awk '{printf "%s\t", $1; for (i=2; i<=NF; i++) printf "%s%s", $i, (i==NF?"\n":"\t")}' >> sula_filtered.dist.names

head -n 2 sula_filtered.minmaf.dist > sula_filtered.minmaf.dist.names
paste <(basename -s .final.bam -a $(cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt)) <(grep '^Ind_' sula_filtered.minmaf.dist | cut -f2-) | awk '{printf "%s\t", $1; for (i=2; i<=NF; i++) printf "%s%s", $i, (i==NF?"\n":"\t")}' >> sula_filtered.minmaf.dist.names

head -n 2 sula_filtered.mind24.dist > sula_filtered.mind24.dist.names
paste <(basename -s .final.bam -a $(cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplebams_filtered.txt)) <(grep '^Ind_' sula_filtered.mind24.dist | cut -f2-) | awk '{printf "%s\t", $1; for (i=2; i<=NF; i++) printf "%s%s", $i, (i==NF?"\n":"\t")}' >> sula_filtered.mind24.dist.names

/home/u15/dannyjackson/programs/FastME-master/src/fastme -i sula_filtered.dist.names -o sula_filtered.tree -q
/home/u15/dannyjackson/programs/FastME-master/src/fastme -i sula_filtered.minmaf.dist.names -o sula_filtered.minmaf.tree -q
/home/u15/dannyjackson/programs/FastME-master/src/fastme -i sula_filtered.mind24.dist.names -o sula_filtered.mind24.tree -q














fastme -T 20 -i testA_8B.dist -s -D 6 -o testA_8B.nwk
# split the input dataset tree from the bootstraped ones:

head -n 1 testA_8B.nwk > testA_8B.main.nwk
tail -n +2 testA_8B.nwk | awk 'NF' > testA_8B.boot.nwk
and, to place supports on the main tree, use RAxML:

raxmlHPC -f b -t testA_8B.main.nwk -z testA_8B.boot.nwk -m GTRCAT -n testA_8B
or RAxML-NG:

raxml-ng --support --tree testA_8B.main.nwk --bs-trees testA_8B.boot.nwk --prefix testA_8B


# OR in R
library(ape)
d <- as.dist(read.table("output.dist"))
tree <- nj(d)
plot(tree)

Fri Jun 13 16:29:16 MST 2025
NC_087521.1_24014288

NC_087521.1_37303840       
Fri Jun 13 16:32:07 MST 2025
265,791,060


26938933 + 24782708 + 23960227 + 18845095 + 18068808 + 17445492 + 15300844 + 13918507 + 13167548 + 10041279 + 9979391 + 8926215 + 8777388 + 8193927 + 8090667 + 7619428 + 5560761 + 4397208 + 3969748 + 1555303 + 1526071 + 898613 + 1101926 + 996000 + 960721 + 5556839

260,579,647

85,784,825
346364472