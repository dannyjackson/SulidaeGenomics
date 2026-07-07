# PCA

# B1_PCA.sh with trans, subset

#!/usr/bin/env bash
#SBATCH --job-name=genolike_all
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/genolike_all%j.out
#SBATCH --mail-type=ALL

# cd /xdisk/mcnew/dannyjackson/sulidae/analyses/pca
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/pca
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
BAMLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allbams.txt
CHROMLIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/sulidae/analyses/angsd_processing/allsnps_popgen.sites_headless.mafs -bam ${BAMLIST} -anc ${REF} -doPlink 2 -out all -doGeno -1 -dopost 1  -nThreads 8 -rf $CHROMLIST


#!/usr/bin/env bash
#SBATCH --job-name=genolike_filtering
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --output=slurm_output/genolike_filtering%j.out
#SBATCH --mail-type=ALL

module load plink

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/

plink --tped /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/all.tped --tfam /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/all.tfam --allow-extra-chr --snps-only 'just-acgt' --indep-pairwise 50kb 1 0.5 --out genolike_filtered --bad-ld

plink --tped /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/all.tped --tfam /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/all.tfam --allow-extra-chr --snps-only 'just-acgt' --extract genolike_filtered.prune.in --out genolike_pruned --make-bed 



# Submitted batch job 3637791



# save the following as: /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt

sample  species
BFBO501 BFBO
BFBO502 BFBO
BFBO503 BFBO
BFBO504 BFBO
BFBO505 BFBO
BRBO201 BRBO
BRBO202 BRBO
BRBO203 BRBO
BRBO205 BRBO
MABO302 MABO
MABO304 MABO
MABO305 MABO
MABO306 MABO
NABO402 NABO
NABO403 NABO
NABO404 NABO
NABO405 NABO
NABO406 NABO
PEBO601 PEBO
PEBO603 PEBO
PEBO604 PEBO
PEBO605 PEBO
PEBO606 PEBO
RFBO101 RFBO
RFBO102 RFBO
RFBO103 RFBO
RFBO104 RFBO
RFBO105 RFBO
RFBO106 RFBO





# Produce pcangsd files

## Analyze all samples

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/pca

module load micromamba
micromamba activate pcangsd

# ================
# Loop through K = 2–12
# ================

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/pca
GENOLIKE_DIR=$BASE/genolike_pruned
OUTDIR=$BASE/all
THREADS=12

for K in {2..12}; do
    echo ">>> Running NGSadmix for K=$K"

    PREFIX=${OUTDIR}/all.k${K}

    pcangsd \
        -p "$GENOLIKE_DIR" \
        -o "$PREFIX" \
        -t $THREADS \
        --admix \
        --admix-K $K \
        > ${PREFIX}.${K}.likelihood.log 2>&1

    echo ">>> Finished K=$K"
done

echo "All K completed."


# ================
# extract likelihoods
# ================

OUTDIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/pca/all

echo -e "K\tlogL" > ${OUTDIR}/admix_likelihoods.tsv

for K in {2..12}; do
    LOG=${OUTDIR}/all.${K}.likelihood.log

    L=$(grep -o "likelihood[^ ]* [^ ]*" "$LOG" | awk '{print $2}')

    echo -e "${K}\t${L}" >> ${OUTDIR}/admix_likelihoods.tsv
done

micromamba deactivate
micromamba activate r_ocelote
Rscript ../plot_admix_likelihoods.R

## First analyze species
module load micromamba
micromamba activate pcangsd
module load plink
groups=(BFBO BRBO MABO NABO PEBO RFBO)


BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/pca
THREADS=2

for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  OUTDIR=$BASE/${g}

  cd ${OUTDIR}

  tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

  awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/${g}_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/${g}.pca.txt 

  plink --bed ${BASE}/genolike_pruned.bed --bim ${BASE}/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/${g}.pca.txt --out genolike_pruned_${g} --make-bed 

  GENOLIKE_DIR=$BASE/${g}/genolike_pruned_${g}

  echo " pcangsd for $g"

  for K in {2..6}; do
      echo ">>> Running NGSadmix for K=$K"

      PREFIX=${OUTDIR}/${g}.k${K}

      pcangsd \
          -p "$GENOLIKE_DIR" \
          -o "$PREFIX" \
          -t $THREADS \
          --admix \
          --admix-K $K \
          > ${PREFIX}.${K}.likelihood.log 2>&1

      echo ">>> Finished K=$K"
  done

  echo "All K completed."

    
  echo -e "K\tlogL" > ${OUTDIR}/admix_likelihoods.tsv

  for K in {2..6}; do
      LOG=${OUTDIR}/${g}.k${K}.${K}.likelihood.log

      L=$(grep -o "likelihood[^ ]* [^ ]*" "$LOG" | awk '{print $2}')

      echo -e "${K}\t${L}" >> ${OUTDIR}/admix_likelihoods.tsv
  done

done



micromamba deactivate
micromamba activate r_ocelote

groups=(BFBO BRBO MABO NABO PEBO RFBO)

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/pca

for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  OUTDIR=$BASE/${g}

  cd ${OUTDIR}

  Rscript ../plot_admix_likelihoods.R
done

## Next analyze species pairs
micromamba deactivate
micromamba activate pcangsd
module load plink
groups=( MABO_NABO PEBO_BFBO RFBO_BRBO)


BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/pca
THREADS=2


for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  
  OUTDIR=$BASE/${g}

  cd ${OUTDIR}

  # split into sp1 (before _) and sp2 (after _)
  sp1=${g%_*}
  sp2=${g#*_}
  echo "    sp1 = $sp1, sp2 = $sp2"


  tail -n +2 /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

  awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp1}_samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/${g}.pca.txt 
  awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/dannyjackson/sulidae/referencelists/${sp2}_samplecodes.txt >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/${g}.pca.txt 


  plink --bed /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/sulidae/analyses/pca/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/dannyjackson/sulidae/referencelists/${g}.pca.txt --out genolike_pruned_${g} --make-bed 

  GENOLIKE_DIR=$BASE/${g}/genolike_pruned_${g}


  echo " pcangsd for $g"

  for K in {2..12}; do
      echo ">>> Running NGSadmix for K=$K"

      PREFIX=${OUTDIR}/${g}.k${K}

      pcangsd \
          -p "$GENOLIKE_DIR" \
          -o "$PREFIX" \
          -t $THREADS \
          --admix \
          --admix-K $K \
          > ${PREFIX}.${K}.likelihood.log 2>&1

      echo ">>> Finished K=$K"
  done

  echo "All K completed."
done



micromamba deactivate
micromamba activate r_ocelote

BASE=/xdisk/mcnew/dannyjackson/sulidae/analyses/pca

for g in "${groups[@]}"; do

  echo ">>> Processing $g"
  OUTDIR=$BASE/${g}

  cd ${OUTDIR}
  echo -e "K\tlogL" > ${OUTDIR}/admix_likelihoods.tsv

  for K in {2..12}; do
      LOG=${OUTDIR}/${g}.k${K}.${K}.likelihood.log

      L=$(grep -o "likelihood[^ ]* [^ ]*" "$LOG" | awk '{print $2}')

      echo -e "${K}\t${L}" >> ${OUTDIR}/admix_likelihoods.tsv
  done


  Rscript ../plot_admix_likelihoods.R
done


# Plot all PCAs
# Establish ks
all 4
RFBO  3
BRBO  4
MABO  3
NABO  3
PEBO  3
BFBO  3
# These are all 2k, but the plots look messy
MABO_NABO 2, 5, 7, 10, 12
PEBO_BFBO 2, 5, 8
RFBO_BRBO 2, 3, 7
# run_all_pca_admix.sh
micromamba deactivate
micromamba activate r_ocelote
BASE="/xdisk/mcnew/dannyjackson/sulidae/analyses/pca"

R_SCRIPT="/xdisk/mcnew/dannyjackson/sulidae/analyses/pca/pca_admix_plot.nonumbers.R"
Rscript "$R_SCRIPT" all "$BASE" 

R_SCRIPT="pca_admix_plot.onesp.R"
Rscript "$R_SCRIPT" RFBO "$BASE" 3
Rscript "$R_SCRIPT" BRBO "$BASE" 4
Rscript "$R_SCRIPT" MABO "$BASE" 3
Rscript "$R_SCRIPT" NABO "$BASE" 3
Rscript "$R_SCRIPT" PEBO "$BASE" 3
Rscript "$R_SCRIPT" BFBO "$BASE" 3

R_SCRIPT="/xdisk/mcnew/dannyjackson/sulidae/analyses/pca/pca_admix_plot.multisp.R"
cd MABO_NABO
mkdir -p plots_k2 plots_k5 plots_k7 plots_k10 plots_k12
Rscript "$R_SCRIPT" MABO_NABO "$BASE" 2
mv *pruned*pdf plots_k2
Rscript "$R_SCRIPT" MABO_NABO "$BASE" 5
mv *pruned*pdf plots_k5
Rscript "$R_SCRIPT" MABO_NABO "$BASE" 7
mv *pruned*pdf plots_k7
Rscript "$R_SCRIPT" MABO_NABO "$BASE" 10
mv *pruned*pdf plots_k10
Rscript "$R_SCRIPT" MABO_NABO "$BASE" 12
mv *pruned*pdf plots_k12

cd PEBO_BFBO
mkdir -p plots_k2 plots_k5 plots_k8 
Rscript "$R_SCRIPT" PEBO_BFBO "$BASE" 2
mv *pruned*pdf plots_k2
Rscript "$R_SCRIPT" PEBO_BFBO "$BASE" 5
mv *pruned*pdf plots_k5
Rscript "$R_SCRIPT" PEBO_BFBO "$BASE" 8
mv *pruned*pdf plots_k8


cd RFBO_BRBO
mkdir -p plots_k2 plots_k3 plots_k7
Rscript "$R_SCRIPT" RFBO_BRBO "$BASE" 2
mv *pruned*pdf plots_k2
Rscript "$R_SCRIPT" RFBO_BRBO "$BASE" 3
mv *pruned*pdf plots_k3
Rscript "$R_SCRIPT" RFBO_BRBO "$BASE" 7
mv *pruned*pdf plots_k7



R_SCRIPT="pca_admix_plot.multisp.R"
groups=( MABO_NABO PEBO_BFBO RFBO_BRBO)
for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  Rscript "$R_SCRIPT" "$g" "$BASE"
done


groups=(BFBO BRBO MABO NABO PEBO RFBO)
R_SCRIPT="pca_admix_plot.onesp.R"
for g in "${groups[@]}"; do
  echo ">>> Processing $g"
  Rscript "$R_SCRIPT" "$g" "$BASE"
done