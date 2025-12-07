 

#!/usr/bin/env bash
#SBATCH --job-name=dsuite_brbo
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/dsuite_brbo%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/RFBO/RFBO_BRBO.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt --ABBAclustering -g


(((CPacific,Cocos),(Caribbean,Atlantic)),Outgroup);

#!/usr/bin/env bash
#SBATCH --job-name=dsuite_brbo
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/dsuite_brbo%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/RFBO/RFBO_BRBO.autosomes.vcf.gz
~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t BRBO.nwk --ABBAclustering -g -n BRBOtree
~/programs/Dsuite/Build/Dsuite Fbranch BRBO.nwk SETS_BRBOtree_tree.txt > BRBO_Fbranch.txt


module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py BRBO_Fbranch.txt BRBO.nwk --tree-label-size 3 --ladderize --use_distances

# pvalue
pvalues <- c(2.30E-16,2.30E-16,2.30E-16,2.30E-16)
2.3e-16
2.3e-16
2.3e-16
2.3e-16

# clustering sensitive
pvalues <- c(7.94E-05,3.21E-11,3.99E-08,5.75E-05)
7.940000e-05
1.284000e-10
7.980000e-08
7.666667e-05
# clustering robust
pvalues <- c(0.000135023,0.154974,0.00159699,0.652178)
0.000540092
0.206632000
0.003193980
0.652178000

p.adjust(pvalues, method = "fdr")