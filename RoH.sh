# Runs of homozygosity
# Install and fix roh-viz script
wget https://raw.githubusercontent.com/samtools/bcftools/refs/heads/develop/misc/roh-viz
perl -pi -e 's/scalar %gt == 1/keys(%gt) == 1/' roh-viz

git clone https://github.com/samtools/bcftools.git
perl -pi -e 's/scalar %gt == 1/keys(%gt) == 1/' bcftools/misc/roh-viz
chmod +x bcftools/misc/roh-viz

for vcf in *.samtools.vcf; do
  # specimen name is the bit before the first dot in the filename
  base="$(basename "$vcf")"
  newsample="${base%%.*}"
  echo "Renaming sample in $vcf to $newsample"
  # sanity check: expect exactly 1 sample in the VCF
  ns=$(bcftools query -l "$vcf" | wc -l)
  if [ "$ns" -ne 1 ]; then
    echo "Skipping $vcf: has $ns samples (expected 1)" >&2
    continue
  fi

  # replace the single sample name with the new one
  out=renamedvcfs/"$vcf"
  printf "%s\n" "$newsample" > "${out%.vcf.gz}.names"
  bcftools reheader -s "${out%.vcf.gz}.names" -o "$out" "$vcf"
  bgzip "$out"
  tabix -f "$out.gz"

  rm -f "${out%.vcf.gz}.names"
  echo "Renamed sample in $vcf -> $newsample"
done

module load micromamba
micromamba activate bcftools_elgato

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/roh


while read -r SAMPLE;
do
VCF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/vcf/renamedvcfs/${SAMPLE}.CM062567.1.samtools.vcf.gz
VCF=PEBO/chr1/input/${SAMPLE}.CM062567.1.samtools.vcf.gz

bcftools roh -G30 --AF-dflt 0.4 $VCF -o ${SAMPLE}.CM062567.out.txt
./roh-viz -i ${SAMPLE}.CM062567.out.txt -v $VCF -o ${SAMPLE}.CM062567.html
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

mkdir PEBO PEBO MABO NABO RFBO BRBO
mkdir -p PEBO/chr1/input PEBO/chr1/output

cp /xdisk/mcnew/dannyjackson/sulidae/datafiles/snpable_masks/vcf/renamedvcfs/PEBO*CM062567* PEBO/chr1/input
./bcftools/misc/run-roh.pl -i PEBO/chr1/input -o PEBO/chr1/output
plot-roh.py -i output-dir






######## Try subsetting wgs and all indv file to chr 1 only ##########

module load micromamba
micromamba activate bcftools_elgato
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/roh

VCF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.autosomes.vcf.gz
# bcftools index $VCF
mkdir -p PEBO/chr4/input PEBO/chr4/output

while read -r SAMPLE;
do
bcftools view -r CM062570.1 -s ${SAMPLE}  -Oz -o PEBO/chr4/input/$SAMPLE.chr4.vcf.gz $VCF
bcftools roh -G30 --AF-dflt 0.4 PEBO/chr4/input/$SAMPLE.chr4.vcf.gz -o PEBO/chr4/output/$SAMPLE.chr4.out.txt
./roh-viz -i PEBO/chr4/output/$SAMPLE.chr4.out.txt -v PEBO/chr4/input/$SAMPLE.chr4.vcf.gz -o PEBO/chr4/output/$SAMPLE.chr4.html
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt

./bcftools/misc/run-roh.pl -i BFBO/chr4/input -o BFBO/chr4/output2

while read -r SAMPLE;
do
bcftools view -r CM062570.1 -s ${SAMPLE}  -Oz -o PEBO/chr4/input/$SAMPLE.chr4.vcf.gz $VCF
bcftools roh -G30 --AF-dflt 0.4 PEBO/chr4/input/$SAMPLE.chr4.vcf.gz -o PEBO/chr4/output/$SAMPLE.chr4.out.txt
./roh-viz -i PEBO/chr4/output/$SAMPLE.chr4.out.txt -v PEBO/chr4/input/$SAMPLE.chr4.vcf.gz -o PEBO/chr4/output/$SAMPLE.chr4.html
done < /xdisk/mcnew/dannyjackson/sulidae/referencelists/PEBO_samplecodes.txt

./bcftools/misc/run-roh.pl -i PEBO/chr4/input -o PEBO/chr4/output2

./bcftools/misc/run-roh.pl \
  -i PEBO/chr4/input \
  -o PEBO/chr4/output2 \
  --roh-args '-G30 --AF-dflt 0.4'

./bcftools/misc/plot-roh.py -o PEBO_chr4_roh.pdf PEBO/chr4/output2
./bcftools/misc/plot-roh.py -o PEBO_chr4_roh.png PEBO/chr4/output2