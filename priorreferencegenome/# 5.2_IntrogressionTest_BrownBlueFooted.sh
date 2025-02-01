# Test for introgression between Brown and Blue-footed boobies


# Test 2a	Brown (3,5)	Brown (1,2)	Blue-footed (All)	Red-footed (all)
cd /data5/sulidae/angsd/BRBF

mkdir testa testb testc testd teste testf
cd testa 

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO501*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO503*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO505*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO506*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 


# Test 2b Blue-footed (1,2,3,5)	Blue-footed (4)	Brown (1)	Red-footed (all)
cd /data5/sulidae/angsd/BRBF/testb


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO501*bam > angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO503*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO505*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BRBF.pops

echo "3" > angsd.BRBF.abba
echo "1" >> angsd.BRBF.abba
echo "1" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 


# Test 2c	Brown (3,5)	Brown (1,2)	Peruvian (All)	Red-footed (all)
cd /data5/sulidae/angsd/BRBF/testc


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO601*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO602*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO603*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO604*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO605*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO606*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "6" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 


# Test 2d	Brown (3,5)	Brown (1,2)	Masked (All)	Red-footed (all)
cd /data5/sulidae/angsd/BRBF/testd


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO302*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO304*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO305*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO306*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "4" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 


# Test 2e	Brown (3,5)	Brown (1,2)	Nazca (All)	Red-footed (all)
cd /data5/sulidae/angsd/BRBF/teste


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam > angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO402*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO403*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO404*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO405*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO406*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BRBF.pops

echo "2" > angsd.BRBF.abba
echo "2" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 



# Test 2f	Peruvian (all)	Blue-footed (all)	Brown (1,2)	Red-footed (all)
cd /data5/sulidae/angsd/BRBF/testf


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO601*bam > angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO602*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO603*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO604*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO605*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO606*bam >> angsd.BRBF.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO501*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO503*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO505*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO506*bam >> angsd.BRBF.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO202*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.BRBF.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BRBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BRBF.pops

echo "6" > angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba
echo "4" >> angsd.BRBF.abba
echo "5" >> angsd.BRBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BRBF.pops -sizeFile angsd.BRBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BRBF.abba 
