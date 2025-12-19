# Sex Inference Plot
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/plot_sex_inference

DEPTH=/xdisk/mcnew/dannyjackson/sulidae/analyses/sexing_by_depth/normalized_depth_all_chromosomes.csv

head -n +1 $DEPTH > chrZ1_depth.csv
grep 'CM062600.1' $DEPTH >> chrZ1_depth.csv

head -n +1 $DEPTH > chrZ2_depth.csv
grep 'CM062595.1' $DEPTH >> chrZ2_depth.csv

DEPTH_Z1=chrZ1_depth.csv
DEPTH_Z2=chrZ2_depth.csv
HET_Z1=/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/hetplot_allchr/CM062600.1.with_HO.tsv
HET_Z2=/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/hetplot_allchr/CM062595.1.with_HO.tsv


head $HET_Z1 
head $DEPTH_Z1

# Extract relevant columns
awk 'NR>1 {print $1, $6}' OFS="\t" $HET_Z1 > HET_Z1.clean
awk 'NR>1 {print $1, $6}' OFS="\t" $HET_Z2 > HET_Z2.clean
awk -F',' 'NR>1 {print $1, $8}' OFS="\t" $DEPTH_Z1 > DEPTH_Z1.clean
awk -F',' 'NR>1 {print $1, $8}' OFS="\t" $DEPTH_Z2 > DEPTH_Z2.clean

# Make matrix for each chromosome
join -t $'\t' HET_Z1.clean DEPTH_Z1.clean \
  | awk '{print "Z1",$1,$2,$3}' OFS="\t" > Z1.matrix

join -t $'\t' HET_Z2.clean DEPTH_Z2.clean \
  | awk '{print "Z2",$1,$2,$3}' OFS="\t" > Z2.matrix

# Combine matrices
echo -e "Chromosome\tInd\tH_O\tnorm_depth" > Z_matrix.csv
cat Z1.matrix Z2.matrix >> Z_matrix.csv
