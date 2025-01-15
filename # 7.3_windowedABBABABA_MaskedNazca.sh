# Windowed ABBA BABA Masked Nazca

cd /data5/sulidae/angsd/MANA/all/windowed

~/programs/angsd/angsd -bam ../angsd.MANA.pops -out MANA -setMinDepthInd 2 -nThreads 8 -doBcf 1 -doGeno 4 -doPost 1 -gl 1 -doMajorMinor 1 -doMaf 2 -rf /data5/sulidae/angsd/reference_lists/autosomes_filtered.txt


bcftools convert -O z -o MANA.vcf.gz MANA.bcf

rm MANA.bcf

# run on mabo nabo
cd /data5/sulidae/angsd/dsuite/MANA

# echo -e 'SAMPLE_ID\tSPECIES_ID' > SETS.txt
cat /data5/sulidae/angsd/MANA/windowed/pops.txt > SETS.txt 
sed -i 's/P1/MABO_allo/g' SETS.txt 
sed -i 's/P2/MABO_sym/g' SETS.txt 
sed -i 's/P3/NABO/g' SETS.txt 
sed -i 's/P4/Outgroup/g' SETS.txt 

mv /data5/sulidae/angsd/MANA/all/windowed/MANA.vcf.gz .

# gunzip MANA.vcf
cd /data5/sulidae/angsd/dsuite/MANA/parallel

~/programs/Dsuite/utils/DtriosParallel SETS.txt /data5/sulidae/angsd/MANA/all/windowed/MANA.vcf.gz -n MANA --use-genotype-probabilities --cores 12

# this one worked !

~/programs/Dsuite/Build/Dsuite Dtrios MANA.vcf.gz SETS.txt -o MANA --use-genotype-probabilities 

~/programs/Dsuite/Build/Dsuite Dinvestigate MANA.vcf.gz SETS.txt --use-genotype-probabilities -w 1000,1000

# but with different d stat than angsd? idk why... i will need to redo with the proper filters including autosomes.filter.txt

MABO_allo       MABO_sym        NABO
D=0.0298134
f_d=0.0403082   397976/9.87332e+06
f_dM=0.0187862  397976/2.11844e+07
ABBA_KSpval = 0
BABA_KSpval = 0




sed -i 's/NC_087513.1/1/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087514.1/2/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087515.1/3/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087516.1/4/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087517.1/5/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087518.1/6/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087519.1/7/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087520.1/8/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087521.1/9/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087522.1/10/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087523.1/11/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087524.1/12/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087525.1/13/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087526.1/14/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087527.1/15/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087528.1/16/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087529.1/17/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087530.1/18/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087531.1/19/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087532.1/20/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087533.1/21/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087534.1/22/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087535.1/23/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087536.1/24/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087537.1/25/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087538.1/26/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087539.1/27/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087540.1/28/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087541.1/29/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087542.1/30/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087543.1/31/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087544.1/32/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087545.1/33/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087546.1/34/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087547.1/90/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt
sed -i 's/NC_087548.1/91/g' MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt



library(qqman)
library(dplyr)
df = read.csv("MABO_allo_MABO_sym_NABO_localFstats__1000_1000.txt", sep='\t')

df$f_d = ifelse(df$D < 0, 0, df$f_d)

df$f_d = ifelse(df$f_d > 100, NA, df$f_d)

# plot histogram of fd
pdf(file = "MANA_fdstats_hist.pdf", width = 5, height = 5, useDingbats=FALSE)
hist(df$f_d, xlim=c(0,1))
  dev.off()

# plot manhattan of fd
dfsubset<-df[complete.cases(df$f_d),]
SNP<-c(1: (nrow(dfsubset)))
mydf<-data.frame(SNP,dfsubset)

pdf(file = "MANA.slidingwindow_fdstats_plot_1.pdf", width = 20, height = 7, useDingbats=FALSE)

manhattan(mydf,chr="chr",bp="windowStart",p="f_d", logp=FALSE,ylab="fd stat",cex = 0.5)

dev.off()

/data5/sulidae/angsd/dsuite/MANA/MANA.slidingwindow_fdstats_plot_1.pdf 













echo -e 'PEBO\tBFBO\tNABO' > BFNA_trios.txt


~/programs/Dsuite/Build/Dsuite Dinvestigate MANA.vcf.gz SETS.txt BFNA_trios.txt --use-genotype-probabilities -w 1000,1000

~/programs/Dsuite/Build/Dsuite Dtrios MANA.vcf.gz SETS.txt -o MANA --use-genotype-probabilities 

echo -e 'MABO_allo\tMABO_sym\tNABO' > MANA_trios.txt

~/programs/Dsuite/Build/Dsuite Dinvestigate MANA.vcf.gz SETS.txt MANA_trios.txt --use-genotype-probabilities -w 1000,1000


PEBO    BFBO    NABO
D=-0.0334056
f_d=-0.1691     -2908.77/17201.5
f_dM=-0.0698076 -2908.77/41668.4
ABBA_KSpval = 6.82526e-08
BABA_KSpval = 0.0145576

head PEBO_BFBO_NABO_localFstats__1000_1000.txt 
chr     windowStart     windowEnd       D       f_d     f_dM    d_f
NC_087538.1     7805    15551   -0.056252       -0.827096       -0.241965       -0.025704
NC_087538.1     15552   21625   -0.008419       -0.051355       -0.019739       -0.002217
NC_087538.1     21626   31676   0.068752        0.481858        0.269870        0.029971
NC_087538.1     31865   38626   0.009821        0.062844        0.031278        0.005345
NC_087538.1     38627   41463   0.218262        0.596644        0.368097        0.089747

sed -i 's/NC_087538.1/26/g' PEBO_BFBO_NABO_localFstats__1000_1000.txt 

library(qqman)
library(dplyr)
df = read.csv("PEBO_BFBO_NABO_localFstats__1000_1000.txt", sep='\t')

df$f_d = ifelse(df$D < 0, 0, df$f_d)

df$f_d = ifelse(df$f_d > 100, NA, df$f_d)

# plot histogram of fd
pdf(file = "NC_087538_fdstats_hist.pdf", width = 5, height = 5, useDingbats=FALSE)
hist(df$f_d, xlim=c(0,1))
  dev.off()

# plot manhattan of fd
dfsubset<-df[complete.cases(df$f_d),]
SNP<-c(1: (nrow(dfsubset)))
mydf<-data.frame(SNP,dfsubset)

pdf(file = "NC_087538.slidingwindow_fdstats_plot_1.pdf", width = 20, height = 7, useDingbats=FALSE)

manhattan(mydf,chr="chr",bp="windowStart",p="f_d", logp=FALSE,ylab="fd stat",cex = 0.5)

dev.off()
