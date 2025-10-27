# abba baba
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

grep -Ev 'CM062595\.1|CM062599\.1|CM062600\.1|CM062610\.1' \
    /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.txt \
    > /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

grep 'CM' bam.abbababa2 | grep -Ev 'CM062595\.1|CM062599\.1|CM062600\.1|CM062610\.1' > autosomes.abbababa2

# TEST 1 MABO VS NABO
## MABO_ATLCAR MABO_INDOPA     NABO  RFBO

#!/usr/bin/env bash
#SBATCH --job-name=MABA
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/MABA.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch MABA.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p MANA
cd MANA
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

ls ${BAMDIR}/MABO304*bam > angsd.MANA.pops
ls ${BAMDIR}/MABO305*bam >> angsd.MANA.pops

ls ${BAMDIR}/MABO302*bam >> angsd.MANA.pops
ls ${BAMDIR}/MABO306*bam >> angsd.MANA.pops

ls ${BAMDIR}/NABO402*bam >> angsd.MANA.pops
ls ${BAMDIR}/NABO403*bam >> angsd.MANA.pops
ls ${BAMDIR}/NABO404*bam >> angsd.MANA.pops
ls ${BAMDIR}/NABO405*bam >> angsd.MANA.pops
ls ${BAMDIR}/NABO406*bam >> angsd.MANA.pops

ls ${BAMDIR}/RFBO101*bam >> angsd.MANA.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.MANA.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.MANA.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.MANA.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.MANA.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.MANA.pops

echo "2" > angsd.MANA.abba
echo "2" >> angsd.MANA.abba
echo "5" >> angsd.MANA.abba
echo "6" >> angsd.MANA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.MANA.pops \
    -sizeFile angsd.MANA.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.MANA.abba 



# TEST 3 BFBO vs PEBO
# BFBO_GofCA BFBO_sym     PEBO  RF


#!/usr/bin/env bash
#SBATCH --job-name=BFPEa
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BFPEa.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BFPE.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BFPEa
cd BFPEa
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

ls ${BAMDIR}/BFBO504*bam > angsd.BFPE.pops

ls ${BAMDIR}/BFBO502*bam >> angsd.BFPE.pops
ls ${BAMDIR}/BFBO505*bam >> angsd.BFPE.pops

ls ${BAMDIR}/PEBO601*bam >> angsd.BFPE.pops
ls ${BAMDIR}/PEBO603*bam >> angsd.BFPE.pops
ls ${BAMDIR}/PEBO604*bam >> angsd.BFPE.pops
ls ${BAMDIR}/PEBO605*bam >> angsd.BFPE.pops
ls ${BAMDIR}/PEBO606*bam >> angsd.BFPE.pops

ls ${BAMDIR}/RFBO101*bam >> angsd.BFPE.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BFPE.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BFPE.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BFPE.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BFPE.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BFPE.pops


echo "1" > angsd.BFPE.abba
echo "2" >> angsd.BFPE.abba
echo "5" >> angsd.BFPE.abba
echo "6" >> angsd.BFPE.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFPE.pops \
    -sizeFile angsd.BFPE.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BFPE.abba 

# TEST 2a BRBO BFBO

#!/usr/bin/env bash
#SBATCH --job-name=BRBF
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BRBF.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBF.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRBF
cd BRBF
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 BRBO 3,5
ls ${BAMDIR}/BRBO203*bam > angsd.BRBF.pops
ls ${BAMDIR}/BRBO205*bam >> angsd.BRBF.pops

# p2 BRBO 1,2
ls ${BAMDIR}/BRBO201*bam >> angsd.BRBF.pops
ls ${BAMDIR}/BRBO202*bam >> angsd.BRBF.pops

# p3 BFBO
ls ${BAMDIR}/BFBO504*bam >> angsd.BRBF.pops
ls ${BAMDIR}/BFBO501*bam >> angsd.BRBF.pops
ls ${BAMDIR}/BFBO502*bam >> angsd.BRBF.pops
ls ${BAMDIR}/BFBO503*bam >> angsd.BRBF.pops
ls ${BAMDIR}/BFBO505*bam >> angsd.BRBF.pops

# p4

ls ${BAMDIR}/RFBO101*bam >> angsd.BRBF.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRBF.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRBF.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRBF.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRBF.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRBF.pops


echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops \
    -sizeFile angsd.BRBF.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 

# TEST 2b

#!/usr/bin/env bash
#SBATCH --job-name=BRBFb
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BRBFb.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBFb.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRBFb
cd BRBFb
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 BFBO southern
ls ${BAMDIR}/BFBO501*bam > angsd.BRBFb.pops
ls ${BAMDIR}/BFBO502*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/BFBO503*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/BFBO505*bam >> angsd.BRBFb.pops

# p2 BFBO GofCA
ls ${BAMDIR}/BFBO504*bam >> angsd.BRBFb.pops


# p3 BRBO GofCA
ls ${BAMDIR}/BRBO201*bam >> angsd.BRBFb.pops


# p4
ls ${BAMDIR}/RFBO101*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRBFb.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRBFb.pops


echo "4" > angsd.BRBFb.abba
echo "1" >> angsd.BRBFb.abba
echo "1" >> angsd.BRBFb.abba
echo "6" >> angsd.BRBFb.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBFb.pops \
    -sizeFile angsd.BRBFb.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBFb.abba 



# TEST 2 d
#!/usr/bin/env bash
#SBATCH --job-name=BRBFd
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BRBFd.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBFd.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRBFd
cd BRBFd
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 BRBO 3,5
ls ${BAMDIR}/BRBO203*bam > angsd.BRBFd.pops
ls ${BAMDIR}/BRBO205*bam >> angsd.BRBFd.pops

# p2 BRBO 1,2
ls ${BAMDIR}/BRBO201*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/BRBO202*bam >> angsd.BRBFd.pops

# p3 MABO
ls ${BAMDIR}/MABO302*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/MABO306*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/MABO304*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/MABO305*bam >> angsd.BRBFd.pops

# p4

ls ${BAMDIR}/RFBO101*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRBFd.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRBFd.pops


echo "2" > angsd.BRBFd.abba
echo "2" >> angsd.BRBFd.abba
echo "4" >> angsd.BRBFd.abba
echo "6" >> angsd.BRBFd.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBFd.pops \
    -sizeFile angsd.BRBFd.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBFd.abba 

# TEST 2 e
#!/usr/bin/env bash
#SBATCH --job-name=BRBFe
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BRBFe.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBFe.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRBFe
cd BRBFe
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 BRBO 3,5
ls ${BAMDIR}/BRBO203*bam > angsd.BRBFe.pops
ls ${BAMDIR}/BRBO205*bam >> angsd.BRBFe.pops

# p2 BRBO 1,2
ls ${BAMDIR}/BRBO201*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/BRBO202*bam >> angsd.BRBFe.pops

# p3 NABO
ls ${BAMDIR}/NABO402*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/NABO403*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/NABO404*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/NABO405*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/NABO406*bam >> angsd.BRBFe.pops

# p4

ls ${BAMDIR}/RFBO101*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRBFe.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRBFe.pops


echo "2" > angsd.BRBFe.abba
echo "2" >> angsd.BRBFe.abba
echo "5" >> angsd.BRBFe.abba
echo "6" >> angsd.BRBFe.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBFe.pops \
    -sizeFile angsd.BRBFe.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBFe.abba 

# TEST 2 f
#!/usr/bin/env bash
#SBATCH --job-name=BRBFf
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BRBFf.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBFf.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRBFf
cd BRBFf
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 PEBO
ls ${BAMDIR}/PEBO601*bam > angsd.BRBFf.pops
ls ${BAMDIR}/PEBO603*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/PEBO604*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/PEBO605*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/PEBO606*bam >> angsd.BRBFf.pops

# p2 BFBO
ls ${BAMDIR}/BFBO501*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/BFBO502*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/BFBO503*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/BFBO504*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/BFBO505*bam >> angsd.BRBFf.pops


# p3 BRBO Pacific
ls ${BAMDIR}/BRBO201*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/BRBO202*bam >> angsd.BRBFf.pops

# p4

ls ${BAMDIR}/RFBO101*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRBFf.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRBFf.pops


echo "5" > angsd.BRBFf.abba
echo "5" >> angsd.BRBFf.abba
echo "2" >> angsd.BRBFf.abba
echo "6" >> angsd.BRBFf.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBFf.pops \
    -sizeFile angsd.BRBFf.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBFf.abba 

# TEST 2g
#!/usr/bin/env bash
#SBATCH --job-name=BRBFg
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBFg.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBFg.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRBFg
cd BRBFg
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 BRBO 3,5
ls ${BAMDIR}/BRBO203*bam > angsd.BRBFg.pops
ls ${BAMDIR}/BRBO205*bam >> angsd.BRBFg.pops

# p2 BRBO 1,2
ls ${BAMDIR}/BRBO201*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/BRBO202*bam >> angsd.BRBFg.pops

# p3 4taxon
ls ${BAMDIR}/PEBO601*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/PEBO603*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/PEBO604*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/PEBO605*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/PEBO606*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/BFBO501*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/BFBO502*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/BFBO503*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/BFBO504*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/BFBO505*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/NABO402*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/NABO403*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/NABO404*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/NABO405*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/NABO406*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/MABO302*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/MABO306*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/MABO304*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/MABO305*bam >> angsd.BRBFg.pops

# p4

ls ${BAMDIR}/RFBO101*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRBFg.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRBFg.pops


echo "2" > angsd.BRBFg.abba
echo "2" >> angsd.BRBFg.abba
echo "19" >> angsd.BRBFg.abba
echo "6" >> angsd.BRBFg.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBFg.pops \
    -sizeFile angsd.BRBFg.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBFg.abba 


# TEST 4 a
# NAZCA vs BFBO

#!/usr/bin/env bash
#SBATCH --job-name=NABF
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/NABF.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch NABF.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p NABF
cd NABF
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 PEBO 
ls ${BAMDIR}/PEBO601*bam > angsd.NABF.pops
ls ${BAMDIR}/PEBO603*bam >> angsd.NABF.pops
ls ${BAMDIR}/PEBO604*bam >> angsd.NABF.pops
ls ${BAMDIR}/PEBO605*bam >> angsd.NABF.pops
ls ${BAMDIR}/PEBO606*bam >> angsd.NABF.pops


# p2 BFBO
ls ${BAMDIR}/BFBO504*bam >> angsd.NABF.pops
ls ${BAMDIR}/BFBO501*bam >> angsd.NABF.pops
ls ${BAMDIR}/BFBO502*bam >> angsd.NABF.pops
ls ${BAMDIR}/BFBO503*bam >> angsd.NABF.pops
ls ${BAMDIR}/BFBO505*bam >> angsd.NABF.pops

# p3 NABO
ls ${BAMDIR}/NABO402*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO403*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO404*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO405*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO406*bam >> angsd.NABF.pops


# p4
ls ${BAMDIR}/RFBO101*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.NABF.pops


echo "5" > angsd.NABF.abba
echo "5" >> angsd.NABF.abba
echo "5" >> angsd.NABF.abba
echo "6" >> angsd.NABF.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.NABF.pops \
    -sizeFile angsd.NABF.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.NABF.abba 

# Test 4b

#!/usr/bin/env bash
#SBATCH --job-name=NABFb
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/NABFb.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch NABFb.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p NABFb
cd NABFb
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
# p1 BFBO  GofCA
ls ${BAMDIR}/BFBO504*bam > angsd.NABF.pops

# p2 BFBO Galapagos
ls ${BAMDIR}/BFBO505*bam >> angsd.NABF.pops

# p3 NABO
ls ${BAMDIR}/NABO402*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO403*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO404*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO405*bam >> angsd.NABF.pops
ls ${BAMDIR}/NABO406*bam >> angsd.NABF.pops


# p4
ls ${BAMDIR}/RFBO101*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.NABF.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.NABF.pops


echo "1" > angsd.NABF.abba
echo "1" >> angsd.NABF.abba
echo "5" >> angsd.NABF.abba
echo "6" >> angsd.NABF.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.NABF.pops \
    -sizeFile angsd.NABF.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.NABF.abba 


# Test 5
# Brown Masked

#!/usr/bin/env bash
#SBATCH --job-name=BRMA
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/BRMA.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRMA.sh
# compiled on puma

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dancingqueen

mkdir -p BRMA
cd BRMA
BASE=/xdisk/mcnew/dannyjackson/sulidae
BAMDIR=$BASE/datafiles/finalbams
REF=$BASE/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# p1 MABO INDOPAC
ls ${BAMDIR}/MABO302*bam > angsd.BRMA.pops
ls ${BAMDIR}/MABO306*bam >> angsd.BRMA.pops


# p2 MABO atlcar
ls ${BAMDIR}/MABO304*bam >> angsd.BRMA.pops
ls ${BAMDIR}/MABO305*bam >> angsd.BRMA.pops

# p3 BRBO atlcar
ls ${BAMDIR}/BRBO203*bam >> angsd.BRMA.pops
ls ${BAMDIR}/BRBO205*bam >> angsd.BRMA.pops

# p4
ls ${BAMDIR}/RFBO101*bam >> angsd.BRMA.pops
ls ${BAMDIR}/RFBO102*bam >> angsd.BRMA.pops
ls ${BAMDIR}/RFBO103*bam >> angsd.BRMA.pops
ls ${BAMDIR}/RFBO104*bam >> angsd.BRMA.pops
ls ${BAMDIR}/RFBO105*bam >> angsd.BRMA.pops
ls ${BAMDIR}/RFBO106*bam >> angsd.BRMA.pops


echo "2" > angsd.BRMA.abba
echo "2" >> angsd.BRMA.abba
echo "2" >> angsd.BRMA.abba
echo "6" >> angsd.BRMA.abba

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRMA.pops \
    -sizeFile angsd.BRMA.abba -doCounts 1 -out bam \
    -anc ${REF} -useLast 1 -setMinDepthInd 2 -nThreads 10 -blocksize 500000 -enhance 1 -rf /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

module load R

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRMA.abba 