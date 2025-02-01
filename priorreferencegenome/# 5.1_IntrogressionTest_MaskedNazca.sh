# Test for introgression between Masked and Nazca boobies

# MA_ATLCAR MA_INDOPA     NA  RF
cd /data5/sulidae/angsd/MANA


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO304*bam > angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO305*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO301*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO302*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/MABO306*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO402*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO403*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO404*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO405*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/NABO406*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.MANA.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.MANA.pops

echo "2" > angsd.MANA.abba
echo "3" >> angsd.MANA.abba
echo "5" >> angsd.MANA.abba
echo "5" >> angsd.MANA.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.MANA.pops -sizeFile angsd.MANA.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.MANA.abba 