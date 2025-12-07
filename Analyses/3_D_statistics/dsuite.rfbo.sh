# # Dsuite RFBO
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


bcftools view \
  -r $(paste -sd, /xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt) -s BRBO201,BRBO202,BRBO203,BRBO205,RFBO101,RFBO102,RFBO103,RFBO104,RFBO105,RFBO106 \
  /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered.recode.vcf.gz \
  -Oz -o RFBO_BRBO.autosomes.vcf.gz

# Dsuite
# Save this as SETS.txt

BRBO201 Outgroup
BRBO202 Outgroup
BRBO203 Outgroup
BRBO205 Outgroup
RFBO101 RFBO_CPacific
RFBO102 RFBO_CPacific
RFBO103 RFBO_EPacific
RFBO104 RFBO_Caribbean
RFBO105 RFBO_Atlantic
RFBO106 RFBO_Indian

#!/usr/bin/env bash
#SBATCH --job-name=dsuite_rfbo
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --output=slurm_output/dsuite_rfbo%A_%a.out
#SBATCH --mail-type=ALL

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/RFBO

VCF=RFBO_BRBO.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dtrios $VCF SETS.txt --ABBAclustering -g





# likelihood f_branch 
pvalues <- c(1,1,1,0.611558542,1,0,0,1,1,0,1,1,1,0.000852933,1,0,5.08482E-14,1,1,0,0,1,1,0,0,1,1,0,0,1,1,4.55858E-13,0,1,1,0,0,1,1,0,1,0,0,1,1,0,1,0,0,1,1,5.60889E-09,0,1,0,1,1,0.005135488,0,1)
p.adjust(pvalues, method = "fdr")


# p value
pvalues <- c(2.30E-16,0.0123448,2.30E-16,2.30E-16,7.52E-06,2.30E-16,2.30E-16,0.000133288,2.30E-16,0.00113385)
3.833333e-16
1.234480e-02
3.833333e-16
3.833333e-16
1.074286e-05
3.833333e-16
3.833333e-16
1.666100e-04
3.833333e-16
1.259833e-03

# clustering sensitive
pvalues <- c(0.000111429,3.58E-08,3.26E-09,2.91E-07,1.23E-09,2.30E-16,3.89E-11,0.000516467,6.50E-08,2.48E-07)
1.238100e-04
7.160000e-08
8.150000e-09
3.637500e-07
4.100000e-09
2.300000e-15
1.945000e-10
5.164670e-04
1.083333e-07
3.542857e-07

# clustering_robust
pvalues <- c(0.7774,4.02E-05,0.525925,0.805883,1.05E-14,0.390454,0.904593,5.89E-08,0.881108,8.36E-11)

p.adjust(pvalues, method = "fdr")

