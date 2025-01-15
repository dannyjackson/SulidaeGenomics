# Windowed ABBA BABA Blue-footed Nazca


# plot abba baba 


while read -r line;
do

samtools index $line

done < bamlist.txt

cd /data5/sulidae/angsd/BFNA/testa/regions

~/programs/angsd/angsd -doAbbababa2 1 -bam ../angsd.BFNA.pops -sizeFile ../angsd.BFNA.abba -doCounts 1 -out NC_087538 -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -useLast 1 -setMinDepthInd 2 -nThreads 8 -blocksize 20000 -enhance 1 -rf /data5/sulidae/angsd/reference_lists/autosomes.txt

Rscript ~/programs/angsd/R/estAvgError.R angsdFile="NC_087538" out="result" sizeFile=../angsd.BFNA.abba 

