# FST Scans between Blue-footed and Peruvian


# BFBO PEBO
cd /data5/sulidae/angsd/fst/BFPE


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO501*bam > angsd.BFBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO503*bam >> angsd.BFBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO504*bam >> angsd.BFBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO505*bam >> angsd.BFBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/BFBO506*bam >> angsd.BFBO.pops


ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO601*bam > angsd.PEBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO602*bam >> angsd.PEBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO603*bam >> angsd.PEBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO604*bam >> angsd.PEBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO605*bam >> angsd.PEBO.pops
ls /data5/sulidae/angsd/bamfiles/sorted_bam_files/indelrealignment/PEBO606*bam >> angsd.PEBO.pops

#first calculate per pop saf for each populatoin
~/programs/angsd/angsd -b angsd.BFBO.pops  -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -out BFBO -dosaf 1 -gl 1 -setMinDepthInd 2 -nThreads 8 -r NC_087518.1

~/programs/angsd/angsd -b angsd.PEBO.pops  -anc /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna -out PEBO -dosaf 1 -gl 1 -setMinDepthInd 2 -nThreads 8 -r NC_087518.1

#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS BFBO.saf.idx PEBO.saf.idx >BFBO.PEBO.ml

#prepare the fst for easy window analysis etc
~/programs/angsd/misc/realSFS fst index BFBO.saf.idx PEBO.saf.idx -sfs BFBO.PEBO.ml -fstout BFBOPEBO.out
#get the global estimate
~/programs/angsd/misc/realSFS fst stats BFBOPEBO.out.fst.idx 

# 50000kb window
~/programs/angsd/misc/realSFS fst stats2 BFBOPEBO.out.fst.idx -win 50000 -step 10000 >slidingwindow

# each snp
~/programs/angsd/misc/realSFS fst stats2 BFBOPEBO.out.fst.idx -win 1 -step 1 >snps




# plotting 


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
grep 'NC_' slidingwindow >> slidingwindow_fst.txt

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst_foroutlier.txt
grep 'NC_' slidingwindow >> slidingwindow_fst_foroutlier.txt

sed -i 's/NC_087513.1/1/g' slidingwindow_fst.txt
sed -i 's/NC_087514.1/2/g' slidingwindow_fst.txt
sed -i 's/NC_087515.1/3/g' slidingwindow_fst.txt
sed -i 's/NC_087516.1/4/g' slidingwindow_fst.txt
sed -i 's/NC_087517.1/5/g' slidingwindow_fst.txt
sed -i 's/NC_087518.1/6/g' slidingwindow_fst.txt
sed -i 's/NC_087519.1/7/g' slidingwindow_fst.txt
sed -i 's/NC_087520.1/8/g' slidingwindow_fst.txt
sed -i 's/NC_087521.1/9/g' slidingwindow_fst.txt
sed -i 's/NC_087522.1/10/g' slidingwindow_fst.txt
sed -i 's/NC_087523.1/11/g' slidingwindow_fst.txt
sed -i 's/NC_087524.1/12/g' slidingwindow_fst.txt
sed -i 's/NC_087525.1/13/g' slidingwindow_fst.txt
sed -i 's/NC_087526.1/14/g' slidingwindow_fst.txt
sed -i 's/NC_087527.1/15/g' slidingwindow_fst.txt
sed -i 's/NC_087528.1/16/g' slidingwindow_fst.txt
sed -i 's/NC_087529.1/17/g' slidingwindow_fst.txt
sed -i 's/NC_087530.1/18/g' slidingwindow_fst.txt
sed -i 's/NC_087531.1/19/g' slidingwindow_fst.txt
sed -i 's/NC_087532.1/20/g' slidingwindow_fst.txt
sed -i 's/NC_087533.1/21/g' slidingwindow_fst.txt
sed -i 's/NC_087534.1/22/g' slidingwindow_fst.txt
sed -i 's/NC_087535.1/23/g' slidingwindow_fst.txt
sed -i 's/NC_087536.1/24/g' slidingwindow_fst.txt
sed -i 's/NC_087537.1/25/g' slidingwindow_fst.txt
sed -i 's/NC_087538.1/26/g' slidingwindow_fst.txt
sed -i 's/NC_087539.1/27/g' slidingwindow_fst.txt
sed -i 's/NC_087540.1/28/g' slidingwindow_fst.txt
sed -i 's/NC_087541.1/29/g' slidingwindow_fst.txt
sed -i 's/NC_087542.1/30/g' slidingwindow_fst.txt
sed -i 's/NC_087543.1/31/g' slidingwindow_fst.txt
sed -i 's/NC_087544.1/32/g' slidingwindow_fst.txt
sed -i 's/NC_087545.1/33/g' slidingwindow_fst.txt
sed -i 's/NC_087546.1/34/g' slidingwindow_fst.txt
sed -i 's/NC_087547.1/90/g' slidingwindow_fst.txt
sed -i 's/NC_087548.1/91/g' slidingwindow_fst.txt

# 90 = W, 91 = Z

# outliers
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
fst <- read.csv('slidingwindow_fst_foroutlier.txt', sep ='\t')
max(fst$fst)
# 0.868438
min(fst$fst)
#  -0.030072

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

# 121678 snps, so top 0.1% would be 122

outlier_fst_disorder <- ordered_fst[1:122,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.579191
max(outlier_fst_disorder$fst)
# 0.868438
write_tsv(outlier_fst, "outlierfst.tsv")
# write_tsv(outlier_fst, "outlierfst_12.tsv")



df <- read.csv('slidingwindow_fst.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("BFBOPEBO.fst.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.579191, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()




# Get list of significant genes

tail -n +2 outlierfst.tsv > outlierfst.headless.tsv



Consider a line in the output file.
Search the gff for genes from the same chromosome.
Only keep genes that are within 50kb of the midpoint of the significant window.
Append those genes to a text file.
Move to the next line in the output file.
# gff
/data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/genomic.gff

# code to find relevant genes from a gff given a list of outlier snps

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 50000))
maxpos=$((midpos + 50000))

grep "$chr" /data5/sulidae/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_50k.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' relevantgenes_50k.txt | sed 's/ID\=gene\-//g' | sort -u > relevantgenenames_50k.txt

done < outlierfst.headless.tsv 