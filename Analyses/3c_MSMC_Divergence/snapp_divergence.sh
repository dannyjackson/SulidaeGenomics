# SNAPP
module load ruby/3.3.4
module load beast2/2.7.0
module load plink
module load bcftools htslib vcftools samtools
packagemanager -add SNAPP    # installs the SNAPP add-on
packagemanager -add CA
packagemanager -add SA
packagemanager -add ORC

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/chromosome1

# Generate input file for SNAPP
# Prepare files
awk '{print $2,$1}' /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.txt | tr ' ' '\t' | sed 's/sample/specimen/' > /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.snapp.txt
echo 'NOGA reference' >> /xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.snapp.txt
echo "(NOGA,(RFBO,(BRBO,((MABO,NABO),(PEBO,BFBO)))));" > sp.tree
# echo -e "lognormal(16,2.07944,0.6)\tstem\tRFBO,BRBO,MABO,NABO,BFBO,PEBO" > sulidae.con.txt
# echo -e "normal(0,6,3)\tstem\tRFBO,BRBO,MABO,NABO,BFBO,PEBO" > sulidae.con.txt
# echo -e "lognormal(16,2.07944,0.6)\tcrown\tNOGA,RFBO,BRBO,MABO,NABO,BFBO,PEBO" > sulidae.con.txt
# echo -e "cladeage(6,34,0.01,0.02,0.2,0.5,0.001,0.01)\tstem\tRFBO,BRBO,MABO,NABO,BFBO,PEBO"  > sulidae.con.txt
# lognormal(16, ln(20-16), 0.6) crown NOGA,RFBO,BRBO,MABO,NABO,BFBO,PEBO

echo -e "normal(16,9,3)\tstem\tRFBO,BRBO,MABO,NABO,BFBO,PEBO" > sulidae.con.txt

SPECIES_TABLE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.snapp.txt
# VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062567.1.qualitysort_filtered_mind2.vcf
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf
bgzip $VCF
bcftools index ${VCF}.gz
mkdir -p files
cd files 
REGIONSFILE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/AUTOSOMES.hiconf.txt
bcftools view -r $(cat $REGIONSFILE | tr '\n' ',') ${VCF}.gz -Oz -o autosomes.vcf.gz
bcftools index autosomes.vcf.gz


###
plink \
  --vcf autosomes.vcf.gz \
  --double-id \
  --allow-extra-chr \
  --snps-only just-acgt \
  --set-missing-var-ids @:#:\$1:\$2 \
  --make-bed \
  --out tmp


plink \
  --bfile tmp \
  --indep-pairwise 100kb 1 0.1 \
  --out pruned --allow-extra-chr --thin 0.01

plink \
  --bfile tmp \
  --extract pruned.prune.in \
  --recode vcf-iid \
  --out autosomes --allow-extra-chr



#!/usr/bin/env bash 
#SBATCH --job-name=makephy
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --ntasks=32
#SBATCH --time=10:00:00
#SBATCH --output=slurm_output/makephy
#SBATCH --mail-type=ALL
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes/files

module load ruby/3.3.4
module load beast2/2.7.0
module load plink
module load bcftools htslib vcftools samtools

python3 ~/programs/vcf2phylip.py3 \
    -i  autosomes.vcf.gz \
    -o autosomes.phy \
    -r 

cd ../ 

# ruby ~/programs/snapp_prep/snapp_prep.rb -p files/autosomes.phy -t $SPECIES_TABLE -c sulidae.con.txt -x autosomes.snapp.xml -s sp.tree -l 50000
(NOGA,(RFBO,((BRBO_AtlCar,BRBO_Pacific),(((MABO_AtlCar,MABO_IndoPacific),NABO),(PEBO,(BFBO_GofCA,BFBO_southern))))));
normal(16,9,3)  stem    RFBO,BRBO_AtlCar,BRBO_Pacific,MABO_AtlCar,MABO_IndoPacific,NABO,BFBO_GofCA,BFBO_southern,PEBO

#!/usr/bin/env bash 
#SBATCH --job-name=snapp_prep
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --ntasks=2
#SBATCH --time=20:00:00
#SBATCH --output=slurm_output/snapp_prep
#SBATCH --mail-type=ALL
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes_pops/
module load ruby/3.3.4
module load beast2/2.7.0
module load plink
module load bcftools htslib vcftools samtools

SPECIES_TABLE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.snapp.pops.txt

ruby ~/programs/snapp_prep/snapp_prep.rb -p ../autosomes/files/autosomes.phy -t $SPECIES_TABLE -c sulidae.con.txt -x autosomes.snapp.xml -s sp.tree -l 2000

echo "Validating XML..."
beast -validate autosomes.snapp.xml

#!/usr/bin/env bash 
#SBATCH --job-name=snapp_prep
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --ntasks=32
#SBATCH --time=100:00:00
#SBATCH --output=slurm_output/snapp_prep
#SBATCH --mail-type=ALL
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes_pops/
module load ruby/3.3.4
module load beast2/2.7.0
module load plink
module load bcftools htslib vcftools samtools

beast autosomes.snapp.xml

# Redo with
#!/usr/bin/env bash 
#SBATCH --job-name=snapp_prep
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --ntasks=2
#SBATCH --time=20:00:00
#SBATCH --output=slurm_output/snapp_prep
#SBATCH --mail-type=ALL
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes_pops_notree
module load ruby/3.3.4
module load beast2/2.7.0
module load plink
module load bcftools htslib vcftools samtools

SPECIES_TABLE=/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.all.snapp.pops.txt

ruby ~/programs/snapp_prep/snapp_prep.rb -p ../autosomes/files/autosomes.phy \
  -t $SPECIES_TABLE -c sulidae.con.txt \
  -x autosomes.snapp.xml -l 2000


#!/usr/bin/env bash 
#SBATCH --job-name=SNAPP_autosomes
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --ntasks=64
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/SNAPP_autosomes
#SBATCH --mail-type=ALL


cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes

module load ruby/3.3.4
module load beast2/2.7.0
packagemanager -add SNAPP    # installs the SNAPP add-on
packagemanager -add CA
packagemanager -add SA
packagemanager -add ORC

echo "Running SNAPP..."
beast -threads 64 autosomes.snapp.xml


#!/usr/bin/env bash 
#SBATCH --job-name=SNAPP_autosomes_pops
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --ntasks=64
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/SNAPP_autosomes_pops
#SBATCH --mail-type=ALL


cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes_pops

module load ruby/3.3.4
module load beast2/2.7.0
packagemanager -add SNAPP    # installs the SNAPP add-on
packagemanager -add CA
packagemanager -add SA
packagemanager -add ORC

echo "Running SNAPP..."
beast -threads 64 autosomes.snapp.xml


treeannotator -file snapp.trees -file sp.MCC2.tree -target sp.tree 
treeannotator -file snapp.trees -file sp.MCC.tree 
treeannotator -target sp.tree -height mean -file snapp.trees -topology Usertargettree sp.MCC.tree

#!/usr/bin/env bash 
#SBATCH --job-name=autosomes_pops_notree
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --ntasks=64
#SBATCH --time=60:00:00
#SBATCH --output=slurm_output/autosomes_pops_notree
#SBATCH --mail-type=ALL
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/divergence_snapp/autosomes_pops_notree
module load ruby/3.3.4
module load beast2/2.7.0
packagemanager -add SNAPP    # installs the SNAPP add-on
packagemanager -add CA
packagemanager -add SA
packagemanager -add ORC

echo "Running SNAPP..."
beast -threads 64 autosomes.snapp.xml

ruby ~/programs/snapp_prep/add_theta_to_log.rb -l snapp.log -t snapp.trees -g 8

treeannotator -file snapp.trees -file sp.MCC2.tree -target sp.tree 
treeannotator -file snapp.trees -file sp.MCC.tree 
treeannotator -target sp.tree -height mean -file snapp.trees -topology Usertargettree sp.MCC.tree


library(ggtree)
library(treeio)
library(ape)

tr <- read.beast("sp.MCC.tree")  # reads BEAST annotations like height, height_95%_HPD

phy <- as.phylo(tr)


p <- ggtree(tr) + 
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = signif(height, 3)),
             hjust = -0.2, size = 3)


H <- max(node.depth.edgelength(phy))               # tree height
node_age <- H - node.depth.edgelength(phy)         # ages for tips+internal nodes

# internal node numbers are (Ntip+1):(Ntip+Nnode)
internal_nodes <- (Ntip(phy) + 1):(Ntip(phy) + phy$Nnode)


ggsave("snapp_autosomes_pops.png", plot = p, width = 6, height = 8, units = "in")


library(ggtree)
library(treeio)

head(as.data.frame(read.beast("sp.MCC.tree")))

tr <- read.beast("sp.MCC.tree")

# Make sure annotations (incl. height) are attached to the ggtree data:
ann <- as.data.frame(tr)

mu <- 1.913e-9   # subs/site/generation
g  <- 8          # years/generation (example)

node_age_years <- node_age / mu * g


p <- ggtree(tr) %<+% ann +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = signif(height, 3)),
             hjust = -0.2, size = 3)

ggsave("snapp_autosomes_pops.png", plot = p, width = 6, height = 8, units = "in")

