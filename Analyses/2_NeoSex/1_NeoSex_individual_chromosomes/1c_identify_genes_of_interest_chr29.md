# Chr 29 Gene Identification
This script documents an interactive exploration of the genes present in different regions of chromosome 29

Regions are defined as consecutive runs of windows in which adjacent windows demonstrate either sex differences or sex similarities in heterozygosity.

To infer region boundaries:
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/hetplot_chr29

module load micromamba
micromamba activate r_ocelote

Rscript 1c_identify_regions.r

```
Subset gff files to chromosome of interest:

```
GFF_RE=/xdisk/mcnew/dannyjackson/sulidae/analyses/repeats/vertebrata/ChrsOfInterest.fa.out.gff
grep 'CM062595.1' $GFF_RE > chr29.RE.gff

GFF_MB=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
grep 'CM062595.1' $GFF_MB | grep 'ID\=gene' > chr29.MB.genes.gff
```
## Search with the PC genome
```
GFF=chr29.PC.genes.gff


awk '$4 >= 1 && $5 <= 1400000' $GFF > region0.PC.diff.gff # PAR
awk '$4 >= 1400000 && $5 <= 2150000' $GFF > region1.PC.diff.gff # NR 
awk '$4 >= 2250000 && $5 <= 2600000' $GFF > region2.PC.diff.gff # PAR
awk '$4 >= 2750000 && $5 <= 3200000' $GFF > region3.PC.diff.gff # NR
awk '$4 >= 3250000 && $5 <= 4100000' $GFF > region4.PC.diff.gff # NR
awk '$4 >= 4100000 && $5 <= 4750000' $GFF > region5.PC.diff.gff # PAR



# pull out genes
for f in region*.PC.diff.gff; do
  awk -F'\t' '
    $3=="gene" {
      name="NA"; desc="NA";
      if (match($9,/;?Name=([^;]+)/,n))  name=n[1];
      if (match($9,/;?description=([^;]+)/,d)) desc=d[1];
      print name "\t" desc
    }' "$f" > "${f%.gff}.PC.names_desc.tsv"
done


awk -v FS='\t' -v OFS='\t' '
$3=="gene"{
  name="NA"; desc="NA";

  # split attributes on semicolons
  n = split($9, a, ";");
  for (i=1; i<=n; i++) {
    key = a[i]; sub(/=.*/, "", key);          # text before first =
    val = a[i]; sub(/^[^=]*=/, "", val);      # text after first =
    if (key=="Name")        name=val;
    else if (key=="description") desc=val;
  }

  print $1, $4, $5, name, desc;
}' "$GFF" > genecoords.PC.tsv

grep -Ei 'PMVK|FDPS|FLAD1' genecoords.PC.tsv > isoprenoidgenes.PC.tsv
grep 'keratin' genecoords.PC.tsv > keratingenes.PC.tsv
grep 'nicastrin' genecoords.PC.tsv > nicastrin.PC.tsv
grep 'cathepsin S' genecoords.PC.tsv > cathepsin.PC.tsv


```

## Search with the Morus bassanus genome.
```
GFF=chr29.MB.genes.gff


awk '$4 >= 1 && $5 <= 1400000' $GFF > region0.MB.diff.gff # PAR -- actually, dropout
awk '$4 >= 1400000 && $5 <= 2150000' $GFF > region1.MB.diff.gff # NR 
awk '$4 >= 2250000 && $5 <= 2600000' $GFF > region2.MB.diff.gff # PAR
awk '$4 >= 2750000 && $5 <= 3200000' $GFF > region3.MB.diff.gff # NR
awk '$4 >= 3250000 && $5 <= 4100000' $GFF > region4.MB.diff.gff # NR
awk '$4 >= 4100000 && $5 <= 4750000' $GFF > region5.MB.diff.gff # PAR



# pull out genes
for f in region*.MB.diff.gff; do
  awk -F'\t' '
    $3=="gene" {
      name="NA"; desc="NA";
      if (match($9,/;?Name=([^;]+)/,n))  name=n[1];
      if (match($9,/;?description=([^;]+)/,d)) desc=d[1];
      print name "\t" desc
    }' "$f" > "${f%.gff}.MB.names_desc.tsv"
done


awk -v FS='\t' -v OFS='\t' '
$3=="gene"{
  name="NA"; desc="NA";

  # split attributes on semicolons
  n = split($9, a, ";");
  for (i=1; i<=n; i++) {
    key = a[i]; sub(/=.*/, "", key);          # text before first =
    val = a[i]; sub(/^[^=]*=/, "", val);      # text after first =
    if (key=="Name")        name=val;
    else if (key=="description") desc=val;
  }

  print $1, $4, $5, name, desc;
}' "$GFF" > genecoords.MB.tsv

grep -Ei 'PMVK|FDPS|FLAD1' genecoords.MB.tsv > isoprenoidgenes.MB.tsv
grep 'keratin' genecoords.MB.tsv > keratingenes.MB.tsv
grep 'nicastrin' genecoords.MB.tsv > nicastrin.MB.tsv
grep 'cathepsin S' genecoords.MB.tsv > cathepsin.MB.tsv

```


# Plot with genes
# chr 29 genes

# GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff
# grep 'CM062595.1' $GFF | grep 'ID\=gene' > chr29.PC.genes.gff


GFF_RE=/xdisk/mcnew/dannyjackson/sulidae/analyses/repeats/vertebrata/ChrsOfInterest.fa.out.gff
grep 'CM062595.1' $GFF_RE > chr29.RE.gff

GFF_MB=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
grep 'CM062595.1' $GFF_MB | grep 'ID\=gene' > chr29.MB.genes.gff
```
```

# Species differences in heterozygosity across Region 3, which demonstrates interspecific differences in recombination

```
HETS=hets_CM062595.1.species_sex_hets_per_window.tsv

awk -v lo=2750000 -v hi=3200000 '
BEGIN{FS=OFS="\t"}
NR==1{next}  # skip header
$4>=lo && $4<=hi {
  sp=$1
  if ($7!="NA" && $8!="NA" && $7!="" && $8!="") {
    d = $8 - $7
    if (!(sp in mn) || d < mn[sp]) mn[sp]=d
    if (!(sp in mx) || d > mx[sp]) mx[sp]=d
  }
}
END{
  print "species","min_female_minus_male","max_female_minus_male"
  for (sp in mn) print sp, mn[sp], mx[sp]
}' "$HETS" | sort -k1,1

```


# Combine NR and PAR genes
# awk '$4 >= 1 && $5 <= 1400000' region0, dropout

# NRR
cat region1.MB.diff.gff region4.MB.diff.gff region3.MB.diff.gff > nonrecombining.MB.diff.gff
cat region1.PC.diff.gff region4.PC.diff.gff region3.PC.diff.gff > nonrecombining.PC.diff.gff

# PAR
cat region2.MB.diff.gff region5.MB.diff.gff > PAR.MB.diff.gff
cat region2.PC.diff.gff region5.PC.diff.gff > PAR.PC.diff.gff

# Get gene coords
## Nonrecombining
for f in nonrecombining*.diff.gff; do
  awk -F'\t' '
    $3=="gene"{
      name="NA";
      if (match($9,/(^|;)Name=([^;]+)/,m)) name=m[2];
      if (name!="NA") print name
    }' "$f" \
  | sort -u > "${f%.gff}.genes.txt"
done

## PAR
for f in PAR*.diff.gff; do
  awk -F'\t' '
    $3=="gene"{
      name="NA";
      if (match($9,/(^|;)Name=([^;]+)/,m)) name=m[2];
      if (name!="NA") print name
    }' "$f" \
  | sort -u > "${f%.gff}.genes.txt"
done
