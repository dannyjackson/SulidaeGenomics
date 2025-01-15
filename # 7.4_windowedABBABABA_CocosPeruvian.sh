14_windowedABBABABA_CocosPeruvian.sh

cd /data5/sulidae/angsd/dsuite/COPE

cd /windowed

~/programs/angsd/angsd -bam /data5/sulidae/angsd/COBF/autosomes/testc/angsd.COBF.pops -out COPE -setMinDepthInd 2 -nThreads 8 -doBcf 1 -doGeno 4 -doPost 1 -gl 1 -doMajorMinor 1 -doMaf 2 # -r NC_087513.1


bcftools convert -O z -o COPE.vcf.gz COPE.bcf

cat ../../COBF/autosomes/testc/angsd.COBF.pops  > SETS.txt 

sed -i 's/\/data5\/sulidae\/angsd\/bamfiles\/sorted\_bam\_files\/indelrealignment\///g' SETS.txt
sed -i 's/\.realigned\.bam//g' SETS.txt


sed -i 's/BRBO201/BRBO201\tBRBO/g' SETS.txt
sed -i 's/BRBO203/BRBO203\tBRBO/g' SETS.txt
sed -i 's/BRBO205/BRBO205\tBRBO/g' SETS.txt
sed -i 's/BRBO202/BRBO202\tCOBO/g' SETS.txt
sed -i 's/PEBO601/PEBO601\tPEBO/g' SETS.txt
sed -i 's/PEBO602/PEBO602\tPEBO/g' SETS.txt
sed -i 's/PEBO603/PEBO603\tPEBO/g' SETS.txt
sed -i 's/PEBO604/PEBO604\tPEBO/g' SETS.txt
sed -i 's/PEBO605/PEBO605\tPEBO/g' SETS.txt
sed -i 's/PEBO606/PEBO606\tPEBO/g' SETS.txt
sed -i 's/RFBO101/RFBO101\tOutgroup/g' SETS.txt
sed -i 's/RFBO102/RFBO102\tOutgroup/g' SETS.txt
sed -i 's/RFBO103/RFBO103\tOutgroup/g' SETS.txt
sed -i 's/RFBO104/RFBO104\tOutgroup/g' SETS.txt
sed -i 's/RFBO105/RFBO105\tOutgroup/g' SETS.txt

echo -e 'BRBO\tCOBO\tPEBO' > COPE_trios.txt

~/programs/Dsuite/Build/Dsuite Dtrios COPE.vcf.gz SETS.txt -o COPE --use-genotype-probabilities 

~/programs/Dsuite/Build/Dsuite Dinvestigate COPE.vcf.gz SETS.txt COPE_trios.txt --use-genotype-probabilities -w 1000,1000








cd /data5/sulidae/angsd/COBF/autosomes/testc

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO201*bam > angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO203*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO205*bam >> angsd.COBF.pops

ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BRBO202*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO601*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO602*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO603*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO604*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO605*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO606*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO101*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO102*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO103*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO104*bam >> angsd.COBF.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/RFBO105*bam >> angsd.COBF.pops

echo "3" > angsd.COBF.abba
echo "1" >> angsd.COBF.abba
echo "6" >> angsd.COBF.abba
echo "5" >> angsd.COBF.abba


~/programs/angsd/angsd -doAbbababa2 1 -bam angsd.COBF.pops -sizeFile angsd.COBF.abba -doCounts 1 -out bam -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt


Rscript ~/programs/angsd/R/estAvgError.R angsdFile="bam" out="result" sizeFile=angsd.COBF.abba 

