# Dsuite brbo 

Set working directory
```
mkdir -p cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO
```
Save this as SETS.txt:
```
BRBO201 CPacific
BRBO202 Cocos
BRBO203 Caribbean
BRBO205 Atlantic
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

echo '(((CPacific,Cocos),(Caribbean,Atlantic)),Outgroup);' > BRBO.nwk
```
Run dsuite trios
```
sbatch 3a_brbo_dsuite.sh
sbatch 3b_brbo_dsuite.sh
```
module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
~/programs/Dsuite/Build/Dsuite Fbranch BRBO.nwk SETS_BRBOtree_tree.txt > BRBO_Fbranch.txt

# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py BRBO_Fbranch.txt BRBO.nwk --tree-label-size 3 --color-cutoff 0.3 --ladderize --use_distances

~/programs/Dsuite/Build/Dsuite Fbranch --Pb-matrix BRBO.nwk SETS_BRBOtree_tree.txt > BRBO_Fbranch.txt

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
' BRBO_Fbranch.txt > BRBO_Fbranch.pvalues.tsv
```

Correct for false discovery rate using the p-values from the f_branch p-value statistic using R
```
R
# read file
x <- read.table("BRBO_Fbranch.pvalues.tsv",
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
1
1
0
0
0
0
1
1

Repeat the above process with the clustering sensitive p-values

```
awk '{print $8}' SETS_BRBOtree_tree.txt

R
# read file
df <- read.table("SETS_BRBOtree_tree.txt",
                header = TRUE)

df$p_fdr <- p.adjust(df$p.value, method = "fdr")
df$cs_fdr <- p.adjust(df$clustering_sensitive, method = "fdr")
df$cr_fdr <- p.adjust(df$clustering_robust, method = "fdr")

write.csv(df,file = "SETS_BRBOtree_tree_fdr.csv", row.names = FALSE)
```
