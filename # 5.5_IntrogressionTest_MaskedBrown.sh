# Test for introgression between Masked and Brown boobies

# Test 5	Masked (1, 2)	Masked (4,5)	Brown (3,5)	Red-footed (all)
cd /data5/sulidae/angsd/MABR

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO302*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO304*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO305*bam >> angsd.MABR.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.MABR.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.MABR.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.MABR.pops

echo "1" > angsd.MABR.abba
echo "2" >> angsd.MABR.abba
echo "2" >> angsd.MABR.abba
echo "5" >> angsd.MABR.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.MABR.pops -sizeFile angsd.MABR.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.MABR.abba 
