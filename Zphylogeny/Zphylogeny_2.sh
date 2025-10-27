# Z Phylogeny

BASE=/xdisk/mcnew/dannyjackson/sulidae
VCF_DIR=$BASE/datafiles/snpable_masks/phasedvcfs/renamedvcfs
OUT=$BASE/analyses/raxml_windows_5kb_snp/Zphylogeny

module load bcftools
module load htslib

SCAF=CM062600.1 
# gather all individuals’ VCFs for this scaffold
ls "$VCF_DIR"/*.CM062600.1.whatshap.samtools.vcf.gz > $OUT/tmp.list.$SCAF

# normalize/merge into a multi-sample VCF
bcftools merge -m none -Oz -o $OUT/$SCAF.merged.vcf.gz \
-l $OUT/tmp.list.$SCAF

tabix -p vcf $OUT/$SCAF.merged.vcf.gz

rm $OUT/tmp.list.$SCAF

# filter vcfs 

bcftools view -m2 -M2 -v snps -i 'MAF>0.04 && QUAL>40' \
-Oz -o $OUT/$SCAF.filtered.vcf.gz \
$OUT/$SCAF.merged.vcf.gz

tabix -p vcf $OUT/$SCAF.filtered.vcf.gz



#!/usr/bin/env bash
#SBATCH --job-name=cellphy_Z
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=470G
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/cellphy_Z.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch cellphy_Z.sh
# compiled on puma

/home/u15/dannyjackson/programs/cellphy/cellphy.sh RAXML \
    --bootstrap \
    --msa /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/Zphylogeny/vcfs/CM062600.1.filtered.maf0.4.vcf.gz \
    --prefix /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/test/CM062600_cellphy \
    --model GTR+G \
    --bs-trees 20 \
    --seed 13 

apptainer exec ~/programs/raxml-ng.v2/raxml-ng_2.sif raxml-ng --consense \
  --tree CM062600_cellphy.raxml.bootstraps \
  --prefix CM062600_cellphy.consensus


raxml-ng --support \
  --tree CM062600_cellphy.raxml.bestTree \
  --bs-trees /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_windows_5kb_snp/test/CM062600_cellphy.raxml.bootstraps \
  --prefix CM062600_cellphy.support
