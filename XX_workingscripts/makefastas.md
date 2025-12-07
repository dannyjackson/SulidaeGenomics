# phased vcfs are all in /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml

ls /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_fastas/

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

mkdir /xdisk/mcnew/dannyjackson/sulidae/datafiles/phased_phylip

~/programs/vcf2phylip/vcf2phylip.py -i /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/vcf3/BFBO501.NC_087513.1..samtools.vcf.gz

# run raxml
raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n sula -q /data5/sulidae/final/to_flightless/raxml/partition_file.txt -s test.phy -T 6 -p 12345 -N 10 ­-b 12345 -V