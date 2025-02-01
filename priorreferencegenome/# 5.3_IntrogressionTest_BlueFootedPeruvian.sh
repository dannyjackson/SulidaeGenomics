
# Test for introgression between Blue-footed and Peruvian boobies

# Test 3	Blue-footed (4)	Blue-footed (2,3)	Peruvian (all)	Red-footed (all)
cd /data5/sulidae/angsd/BFPE

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam > angsd.BFPE.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO502*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO503*bam >> angsd.BFPE.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO601*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO603*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO604*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO605*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO606*bam >> angsd.BFPE.pops



ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BFPE.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BFPE.pops

echo "1" > angsd.BFPE.abba
echo "1" >> angsd.BFPE.abba
echo "5" >> angsd.BFPE.abba
echo "5" >> angsd.BFPE.abba

echo "BFBO_CA" > angsd.BFPE.abba.name
echo "BFBO_peru" >> angsd.BFPE.abba.name
echo "PEBO" >> angsd.BFPE.abba.name
echo "RFBO" >> angsd.BFPE.abba.name

~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFPE.pops -sizeFile angsd.BFPE.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BFPE.abba