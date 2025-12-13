# # Dsuite MABO

Set working directory
```
mkdir -p cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/MABO
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/MABO
```
Save this as SETS.txt:
```
MABO302 CPacific
MABO304 Caribbean
MABO305 Atlantic
MABO306 Indian
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup
```
Clean up the SETS file:
```
tr ' ' '\t' < SETS.txt > SETS.tmp
mv SETS.tmp SETS.txt
```
# filter vcf to just autosomes
```
sbatch 4a_filtervcf_mabo.sh 
```

Save this as SETS.txt
```
MABO302 CPacific
MABO304 Caribbean
MABO305 Atlantic
MABO306 Indian
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup
```

Clean up the SETS file and make a tree file:
```
tr ' ' '\t' < SETS.txt > SETS.tmp
mv SETS.tmp SETS.txt

echo '(((CPacific,Indian),(Caribbean,Atlantic)),Outgroup);' > MABO.nwk
```

```
sbatch 4b_mabo_dsuite.sh
sbatch 4c_mabo_dsuite.sh
```


module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
~/programs/Dsuite/Build/Dsuite Fbranch MABO.nwk SETS_MABOtree_tree.txt > MABO_Fbranch.txt
# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py --color-cutoff 0.3 MABO_Fbranch.txt MABO.nwk --tree-label-size 3 --ladderize --use_distances

~/programs/Dsuite/Build/Dsuite Fbranch --Pb-matrix MABO.nwk SETS_MABOtree_tree.txt > MABO_Fbranch.txt
```

Collect the p-value matrix for FDR correction:

```
awk '
BEGIN { grab=0 }
{
    if ($0 ~ /^# p-values:/) {
        grab=1
        next
    }
    if (grab) print
}
' MABO_Fbranch.txt > MABO_Fbranch.pvalues.tsv
```

Correct for false discovery rate using the p-values from the f_branch p-value statistic using R
```
R
# read file
x <- read.table("MABO_Fbranch.pvalues.tsv",
                header = TRUE,
                stringsAsFactors = FALSE,
                check.names = FALSE,
                fill = TRUE,
                comment.char = "")

# keep only numeric columns (drop metadata columns)
mat <- x[, -(1:3)]

# convert "-nan" to NA
mat[mat == "-nan"] <- NA

# coerce to numeric
mat <- apply(mat, 2, as.numeric)

# flatten row-wise, dropping NA
p.values <- as.vector(t(mat))
p.values <- p.values[!is.na(p.values)]

p.values

adj <- p.adjust(p.values, method = "fdr")
cat(adj, sep = "\n")
```
Output:
0
0
1
1
0
0
1
1

Repeat the above process with the clustering sensitive p-values

```

R
# read file
df <- read.table("SETS_MABOtree_tree.txt",
                header = TRUE)

df$p_fdr <- p.adjust(df$p.value, method = "fdr")
df$cs_fdr <- p.adjust(df$clustering_sensitive, method = "fdr")
df$cr_fdr <- p.adjust(df$clustering_robust, method = "fdr")

write.csv(df,file = "SETS_MABOtree_tree_fdr.csv", row.names = FALSE)
```
