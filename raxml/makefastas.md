mkdir -p /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas

cat /xdisk/mcnew/dannyjackson/sulidae/referencelists/*samplecodes.txt > /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

while read -r SCAFFOLD; do

sbatch --account=mcnew \
    --job-name=makefasta.test \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.makefasta.test.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=1:00:00 \
    makefastas.sh

done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/SCAFFOLDS.txt

# run raxml
raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n sula_flightless_b -q /data5/sulidae/final/to_flightless/raxml/partition_file.txt -s /data5/sulidae/my_datasets/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V