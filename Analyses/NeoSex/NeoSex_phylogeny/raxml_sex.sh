# RAxML Analysis -- WG

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_sex
############################################
# Z1
############################################
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_sex/Z1


VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062600.1.qualitysort_filtered_mind2.vcf

awk '$2 == "M" {print $1}' /xdisk/mcnew/dannyjackson/sulidae/referencelists/samples.sex > /xdisk/mcnew/dannyjackson/sulidae/referencelists/males.txt
awk '$2 == "F" {print $1}' /xdisk/mcnew/dannyjackson/sulidae/referencelists/samples.sex > /xdisk/mcnew/dannyjackson/sulidae/referencelists/females.txt

vcftools --vcf $VCF --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/males.txt --out CM062600.1.qualitysort_filtered_mind2.males --recode
vcftools --vcf $VCF --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/females.txt --out CM062600.1.qualitysort_filtered_mind2.females --recode

VCF=CM062600.1.qualitysort_filtered_mind2.males.recode.vcf
PHY=CM062600.1.qualitysort_filtered_mind2.males.phy

python ~/programs/vcf2phylip.py3 \
    -i $VCF \
    -o $PHY \
    -r 


python ~/programs/ascbias.py -p $PHY -o CM062600.1.qualitysort_filtered_mind2.males.noINV.phy

module load raxml-ng

raxml-ng --all \
    --msa CM062600.1.qualitysort_filtered_mind2.males.noINV.phy \
    --model GTR+G+ASC_LEWIS \
    --prefix trees.CM062600.1 \
    --bs-trees 200 \
    --seed 13 \
    --redo


############################################
# Z2
############################################
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_sex/Z2


VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062595.1.qualitysort_filtered_mind2.vcf

/xdisk/mcnew/dannyjackson/sulidae/referencelists/samples.sex 

awk '$2 == "M" {print $1}' /xdisk/mcnew/dannyjackson/sulidae/referencelists/samples.sex > /xdisk/mcnew/dannyjackson/sulidae/referencelists/males.txt
awk '$2 == "F" {print $1}' /xdisk/mcnew/dannyjackson/sulidae/referencelists/samples.sex > /xdisk/mcnew/dannyjackson/sulidae/referencelists/females.txt


vcftools --vcf $VCF --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/males.txt --out CM062595.1.qualitysort_filtered_mind2.males --recode
vcftools --vcf $VCF --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/females.txt --out CM062595.1.qualitysort_filtered_mind2.females --recode

VCF=CM062595.1.qualitysort_filtered_mind2.males.recode.vcf
PHY=CM062595.1.qualitysort_filtered_mind2.males.phy

python ~/programs/vcf2phylip.py3 \
    -i $VCF \
    -o $PHY \
    -r 


python ~/programs/ascbias.py -p $PHY -o CM062595.1.qualitysort_filtered_mind2.males.noINV.phy

module load raxml-ng

raxml-ng --all \
    --msa CM062595.1.qualitysort_filtered_mind2.males.noINV.phy \
    --model GTR+G+ASC_LEWIS \
    --prefix trees \
    --bs-trees 200 \
    --seed 13 \
    --redo
