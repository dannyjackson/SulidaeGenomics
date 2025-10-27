# do it with chr 4

	CM062570.1

#!/bin/bash
BASE=/xdisk/mcnew/dannyjackson/sulidae
OUT=$BASE/analyses/raxml_windows_5kb_snp
WINBED=$OUT/Chr4phylogeny/windows/Chr4.5kb.callable80.bed
VCF_DIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/merged_scaffold_vcfs/
OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/Chr4phylogeny/windows
mkdir -p "$OUTDIR" "$OUTDIR/phy" "$OUTDIR/vcf" "$OUTDIR/trees" "$OUTDIR/phy_noINV"

module load bcftools
module load htslib/1.19.1
module load python/3.11/3.11.4

while read -r SCAF START END _; do
  echo "Processing $SCAF:$START-$END"

  bcftools view -r "$SCAF:$START-$END" -Oz -o "$OUTDIR/vcf/${SCAF}.${START}.${END}.vcf.gz" --threads "$THREADS" "$VCF_DIR/${SCAF}.filtered.vcf.gz"

done < "$WINBED" 

while read -r SCAF START END _; do

  python ~/programs/vcf2phylip.py3 \
    -i $OUTDIR/vcf/${SCAF}.${START}.${END}.vcf.gz \
    -o $OUTDIR/phy/${SCAF}.${START}.${END}.phy \
    -r 
done < "$WINBED" 




while read -r SCAF START END _; do
  mkdir -p $OUTDIR/phy_noINV/

  sites=`head -n 1 $OUTDIR/phy/${SCAF}.${START}.${END}.phy | awk '{print $2}'` 

  if [ "$sites" -eq "0" ]; then 
    echo "no sites left"
    #cat gene${gene}_miss${miss}_mac${mac}.REF.phy > gene${gene}_miss${miss}_mac${mac}.REF.noInv.phy
  else 
    python ~/programs/ascbias.py -p $OUTDIR/phy/${SCAF}.${START}.${END}.phy -o $OUTDIR/phy_noINV/${SCAF}.${START}.${END}.noINV.phy
    rm -f $OUTDIR/phy_noINV/${SCAF}.${START}.${END}.noINV.phy.felsenstein
    rm -f $OUTDIR/phy_noINV/${SCAF}.${START}.${END}.noINV.phy.stamatakis
  fi 


  if [ "$sites" -eq "0" ]; then
    nsnps=0
    rm -f $OUTDIR/phy_noINV/${SCAF}.${START}.${END}.noINV.phy
  else
    nsnps=`cat $OUTDIR/phy_noINV/${SCAF}.${START}.${END}.noINV.phy | head -n 1 | awk '{print $2}'`
  fi
  echo "${SCAF},${START},${END},${nsnps}" >> ${OUTDIR}/SNP_counts_per_window.csv

done < "$WINBED"  


awk -F',' '$4 >= 10' SNP_counts_per_window.csv > SNP_counts_per_window.filtered.csv

awk -F',' '{print $1"."$2"."$3".noINV.phy"}' SNP_counts_per_window.filtered.csv > 10SNPphy.txt


while read -r SCAF START END _; do
  echo "Analyzing $SCAF:$START-$END"

PHY=$OUTDIR/phy_noINV/${SCAF}.${START}.${END}.noINV.phy
OUT=$OUTDIR/trees/${SCAF}.${START}.${END}.noINV

apptainer exec ~/programs/raxml-ng.v2/raxml-ng_2.sif raxml-ng --all \
    --msa "$PHY" \
    --model GTR+G+ASC_LEWIS \
    --prefix "$OUT" \
    --threads 1 \
    --seed 9301 \
    --bs-metric rbs

done < "$WINBED"  



                                                                                  
#!/usr/bin/env bash
#SBATCH --job-name=PreparePhParGene10SNPChr4ylips
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=36
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/ParGene10SNPChr4.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch ParGene10SNPChr4.sh

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate pargenes 

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate pargenes 
# micromamba install bioconda::newick_utils

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/Chr4phylogeny/windows/ParGene_10SNP/
msa_dir=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/Chr4phylogeny/windows/phy_noINV

pargenes.py -a ${msa_dir} -o ${OUTDIR} -c 36 -d nt -R "--model GTR+G+ASC_LEWIS" --use-astral --msa-filter /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/Chr4phylogeny/windows/10SNPphy.txt




# trim 

find /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/Chr4phylogeny/windows/ParGene_10SNP/mlsearch_run/results -name "*.bestTree" \
  -exec awk '{printf "%s\n", $0}' {} \; > allbesttrees.10SNP.tre

nw_ed  allbesttrees.10SNP.tre 'i & b<=0.25' o > allbesttrees.bs-25.10SNP.tre


java -jar ~/programs/ASTRAL/Astral/astral.5.7.8.jar \
  -i allbesttrees.bs-25.10SNP.tre \
  -o astral_species_tree.tre

cat astral_species_tree.tre