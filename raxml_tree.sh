
sula_filtered_mind2.vcf


#!/bin/bash

cd /xdisk/mcnew/dannyjackson/sulidae/datafiles/genotype_calls

~/programs/vcf2phylip/vcf2phylip.py -i sula_filtered_mind2.vcf




/xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt

while read old new; do
  sed -i "s/$old/$new/g" sula_filtered_mind2.min4.phy
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt

sbatch --account=mcnew \
    --job-name=make_phy \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.make_phy.%j \
    --nodes=1 \
    --time=24:00:00 \
    make_phy.sh

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml/

mkdir pruned

#!/bin/bash
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml/pruned
python ~/programs/sula/filter_invariants_all.py /xdisk/mcnew/dannyjackson/sulidae/datafiles/genotype_calls/sula_filtered_mind2.min4.phy
mv variantsites.phy ../variantsites_mind2.phy
mv variantsites_kept.txt ../variantsites_mind2_kept.txt

# first, check that the alignment file can be read
raxml-ng --check --msa variantsites_mind2.phy --model GTR+G --prefix T1
raxml-ng --parse --msa variantsites_mind2.phy --model GTR+G --prefix T2


# infer tree with default settings
 raxml-ng --msa variantsites_mind2.phy --model GTR+G --prefix T3 --threads 8 --seed 2 
