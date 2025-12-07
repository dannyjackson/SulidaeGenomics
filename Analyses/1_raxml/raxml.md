# Scripts to model phylogenetic relationships using RAxML

These scripts use the output from the 0_variant_calling module.

First, convert every chromosome vcf to a phylip file using vcf2phylip.py3 and remove invarant sites with ascbias.py.

```
sbatch --array=1-35 chrom_to_phylip.sh 
```

Then, run RAxML on the output files.
```
sbatch --array=1-35 1_raxml_chroms.sh
```

Combine output files (run interactively)

```
# collect all trees
cat /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/trees/autosomes/*bestTree \
  > /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/treecomparison/allbestTree.autosomes.nwk

# collect number of sites used in each tree

printf "chromosome\tnumber\tsites\tpatterns\n" > summary.tsv

i=0

for f in $(ls /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/trees/autosomes/*.log 2>/dev/null | sort); do
  i=$((i+1))
  # chromosome = filename without trailing ".raxml.log" or ".log"
  chr=$(basename "$f")
  chr=${chr%.raxml.log}
  chr=${chr%.log}

  # extract "Alignment sites / patterns: X / Y"
  read -r sites patterns < <(
    awk -F: '
      /Alignment sites[[:space:]]*\/[[:space:]]*patterns/ {
        split($2,a,"/");
        s=a[1]; p=a[2];
        gsub(/^[ \t]+|[ \t]+$/,"",s);
        gsub(/^[ \t]+|[ \t]+$/,"",p);
        print s, p;
        exit
      }' "$f"
  )

  # fallback if not found
  : "${sites:=NA}" ; : "${patterns:=NA}"

  printf "%s\t%d\t%s\t%s\n" "$chr" "$i" "$sites" "$patterns" >> summary.tsv
done

echo "Wrote summary.tsv"
```

Use R to generate various plots of the RAxML output trees (run interactively).
```
module load micromamba

micromamba activate /xdisk/mcnew/dannyjackson/.local/share/mamba/envs/r_ocelote

Rscript 2_raxml_plotting.r
```