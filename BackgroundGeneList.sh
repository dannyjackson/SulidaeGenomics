# Make Background Gene List

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

grep 'CM' "$GFF" | grep -Ev 'CM062595|CM062599|CM062600|CM062610' | grep 'ID\=gene' | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[3])}' | sed 's/Name\=//g' | sort -u > /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/MorBas.autosome_genes.txt

grep 'CM' "$GFF" \
  | grep -Ev 'CM062595|CM062599|CM062600|CM062610' \
  | grep 'ID=gene' \
  | awk -F'\t' '{
      OFS = "\t";
      n = split($9, a, ";");
      id=""; name="";
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
