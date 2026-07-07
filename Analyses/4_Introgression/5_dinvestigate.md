# Analyze windowed patterns of introgression
Set up directories
```
mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate
mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo
mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/bfbo_pebo
mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/brbo_4taxa
```
## Analyze masked and Nazca boobies
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo

echo -e 'MABO_AtlCar\tMABO_IndoPacific\tNABO' > test_trios.txt

grep 'MABO' /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/SETS.txt > SETS.txt
grep 'NABO' /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/SETS.txt >> SETS.txt
grep 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/SETS.txt >> SETS.txt

#!/usr/bin/env bash
#SBATCH --job-name=MABO_NABO_dinvestigate
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/MABO_NABO_dinvestigate.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch MABO_NABO_dinvestigate.sh

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dinvestigate -w 5000,200 $VCF SETS.txt test_trios.txt

```
Plot it

```
Rscript ../localFstats_manhattan.R \
  MABO_AtlCar_MABO_IndoPacific_NABO_localFstats__5000_200.txt \
  /xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt \
  0.999
```
Create a bedfile of the outlier regions
```
FILE=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/mabo_nabo/analyses/f_dM/MABO_AtlCar_MABO_IndoPacific_NABO/50_25/MABO_AtlCar_MABO_IndoPacific_NABO.f_dM.50_25.0.9999.outliers.tsv

awk -v OFS='\t' 'NR>1 {
    start = $5 - 50000
    if (start < 0) start = 0        # BED cannot have negative starts
    end   = $6 + 50000
    print $2, start, end
}' "$FILE" > outliers.bed

module load bedtools2
GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

bedtools intersect -a ${GFF} -b outliers.bed -wa > outliers.gff

grep 'ID=gene' outliers.gff | sort -u | awk -F'\t' '{print $9}'
 > outliers.genes.txt
grep 'CM' outliers.gff \
  | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID=gene' \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($9, a, ";");
      id=""; name="";
      for (i=1; i<=n; i++) {
        if (a[i] ~ /^description=/) id=a[i];
      }
      print id;
    }' |  awk '{FS = "="} {print $2}' | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MABO_NABO_introgression.tsv


```
## Analyze Peruvian and blue-footed boobies
```

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/bfbo_pebo

echo -e 'BFBO_GofCA\tBFBO_southern\tPEBO' >> test_trios.txt

grep 'BFBO' /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/SETS.txt > SETS.txt
grep 'PEBO' /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/SETS.txt >> SETS.txt
grep 'RFBO' /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/SETS.txt >> SETS.txt
```

```
#!/usr/bin/env bash
#SBATCH --job-name=BFBO_PEBO_dinvestigate
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BFBO_PEBO_dinvestigate.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BFBO_PEBO_dinvestigate.sh

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/bfbo_pebo

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dinvestigate -w 5000,200 $VCF SETS.txt test_trios.txt

```
Plot it
```
Rscript ../localFstats_manhattan.R \
  BFBO_GofCA_BFBO_southern_PEBO_localFstats__5000_200.txt \
  /xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt \
  0.999
```
## Analyze BRBO 4 Taxa
```

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/brbo_4taxa

echo -e 'BRBO_Pacific\tBRBO_AtlCar\t4taxa' > test_trios.txt

cp ../../calls/BRBO_4taxa/SETS.txt SETS.txt

```

```
#!/usr/bin/env bash
#SBATCH --job-name=BRBO_4taxa_dinvestigate
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_output/BRBO_4taxa_dinvestigate.%A_%a.out
#SBATCH --mail-type=ALL
# sbatch BRBO_4taxa_dinvestigate.sh

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/dinvestigate/brbo_4taxa

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

~/programs/Dsuite/Build/Dsuite Dinvestigate -w 5000,200 $VCF SETS.txt test_trios.txt

```
Plot it
```
Rscript ../localFstats_manhattan.R \
  BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__5000_200.txt \
  /xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt \
  0.999
```