# Make Background Gene List

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

grep 'CM' "$GFF" | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID\=gene' \
  | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' \
  | sed 's/Name\=//g' \
  | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.txt

grep 'CM' "$GFF" \
  | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID=gene' \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($9, a, ";");
      id=""; name="";`
      for (i=1; i<=n; i++) {
        if (a[i] ~ /^Name=/) name=a[i];
      }
      print id, name;
    }' |  awk '{FS = "="} {print $2}' | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.txt


grep 'CM' "$GFF" \
  | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID=gene' | grep 'Name=bMorBas2_egapxtmp' | grep 'description' \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($9, a, ";");
      id=""; name="";
      for (i=1; i<=n; i++) {
        if (a[i] ~ /^Name=/) name=a[i];
        if (a[i] ~ /^description=/)   id=a[i];
      }
      print id, name;
    }' | sort -u > MorBas.autosome_genes.unnamed.txt



CONVERSIONFILE=/xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder/bMorBas2_to_PhaCar_proteins_genes.simple.tsv
GENELIST=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.txt

head $CONVERSIONFILE
head $GENELIST

awk 'NR==FNR { conv[$1] = $2; next } 
     { if ($1 in conv) print conv[$1]; else print $1 }' \
     "$CONVERSIONFILE" "$GENELIST" > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.converted.txt


awk '
  NR==FNR { conv[$1] = $2; next }
  {
    key = $1
    if (key in conv) {
      print conv[key] > "/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.converted.txt"        # output replaced entry
      print key, conv[key] > "/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.matched_pairs.txt"    # record match from conversion file
    } else {
      print key > "/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.converted.txt"              # no replacement
    }
  }
' "$CONVERSIONFILE" "$GENELIST"

sed -i 's/,//g' MorBas.autosome_genes.converted.txt