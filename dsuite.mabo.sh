# # Dsuite MABO
# filter vcf to just autosomes

#!/usr/bin/env bash
#SBATCH --job-name=filtervcf
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/filtervcf%A_%a.out
#SBATCH --mail-type=ALL


module load bcftools

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/MABO

bcftools view \
  -r $(paste -sd, /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt) \
  -s MABO302,MABO304,MABO305,MABO306,RFBO101,RFBO102,RFBO103,RFBO104,RFBO105,RFBO106 \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered.recode.vcf.gz \
  -Oz -o RFBO_MABO.autosomes.vcf.gz


# Dsuite
# Save this as SETS.txt

MABO302 CPacific
MABO304 Caribbean
MABO305 Atlantic
MABO306 Indian
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup


#!/usr/bin/env bash
#SBATCH --job-name=dsuite_mabo
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/dsuite_mabo%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/MABO

VCF=RFBO_MABO.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt --ABBAclustering -g


# pvalue
pvalues <- c(2.82E-13,2.30E-16,2.30E-16,2.30E-16)
2.820000e-13
3.066667e-16
3.066667e-16
3.066667e-16
# clustering sensitive
pvalues <- c(2.30E-16,3.56E-06,2.30E-16,0.00585484)
4.600000e-16
4.746667e-06
4.600000e-16
5.854840e-03

# clustering robust
pvalues <- c(1.57E-10,0.360298,1.83E-05,0.00852127)

p.adjust(pvalues, method = "fdr")
6.280000e-10
3.602980e-01
3.660000e-05
1.136169e-02




(((CPacific,Indian),(Caribbean,Atlantic)),Outgroup);

#!/usr/bin/env bash
#SBATCH --job-name=dsuite_mabo
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/dsuite_mabo%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/MABO

VCF=RFBO_MABO.autosomes.vcf.gz
~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt -t MABO.nwk --ABBAclustering -g -n MABOtree
~/programs/Dsuite/Build/Dsuite Fbranch MABO.nwk MABOtree_tree.txt > MABO_Fbranch.txt
