# Test for introgression between Blue-footed and Nazca boobies

# Test 4a	Peruvian (all)	Blue-footed (all)	Nazca (all)	Red-footed (all)
cd /data5/sulidae/angsd/BFNA/testa

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO601*bam > angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO603*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO604*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO605*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO606*bam >> angsd.BFNA.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO501*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO503*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO505*bam >> angsd.BFNA.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO402*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO403*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO404*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO405*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO406*bam >> angsd.BFNA.pops



ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BFNA.pops

echo "5" > angsd.BFNA.abba
echo "4" >> angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFNA.pops -sizeFile angsd.BFNA.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BFNA.abba 



# Test 4b	Blue-footed (4)	Blue-footed (5)	Nazca (all)	Red-footed (all)
cd /data5/sulidae/angsd/BFNA/testb

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam > angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO505*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO402*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO403*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO404*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO405*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO406*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.BFNA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.BFNA.pops

echo "1" > angsd.BFNA.abba
echo "1" >> angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba
echo "5" >> angsd.BFNA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.BFNA.pops -sizeFile angsd.BFNA.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.BFNA.abba 
