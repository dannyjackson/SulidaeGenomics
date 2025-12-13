# # Dsuite RFBO
Set working directory
```
mkdir -p cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/RFBO
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/RFBO
```
Save this as SETS.txt:
```
BRBO201 Outgroup
BRBO202 Outgroup
BRBO203 Outgroup
BRBO205 Outgroup
RFBO101 RFBO_CPacific
RFBO102 RFBO_CPacific
RFBO103 RFBO_EPacific
RFBO104 RFBO_Caribbean
RFBO105 RFBO_Atlantic
RFBO106 RFBO_Indian
```
Clean SETS file:
```
tr ' ' '\t' < SETS.txt > SETS.tmp
mv SETS.tmp SETS.txt

```
# filter vcf to just samples of interest and autosomes
```
sbatch 2a_filtervcf_rfbo.sh
```
# run dsuite
```
echo '((RFBO_Caribbean,(RFBO_Indian,(RFBO_Atlantic,(RFBO_EPacific,RFBO_CPacific)))),Outgroup);' > rfbo.nwk

sbatch 2b_rfbo_dsuite.sh
```

Compute f-branch statistics and a plot to visualize the output using scripts from Dsuite:

```

~/programs/Dsuite/Build/Dsuite Fbranch rfbo.nwk SETS_RFBOtree_tree.txt  > rfbo_Fbranch.txt


module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py rfbo_Fbranch.txt rfbo.nwk --tree-label-size 3 --ladderize --use_distances --color-cutoff 0.3

# Generate p-values
~/programs/Dsuite/Build/Dsuite Fbranch --Pb-matrix rfbo.nwk SETS_RFBOtree_tree.txt > rfbo_Fbranch.txt

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
' rfbo_Fbranch.txt > rfbo_Fbranch.pvalues.tsv
```

Correct for false discovery rate using the p-values from the f_branch p-value statistic using R
```
R
# read file
x <- read.table("rfbo_Fbranch.pvalues.tsv",
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

1
0
0
0
1
1
0
8.402496e-11
8.9002e-08
1
1
1

```
Repeat the above process with the clustering sensitive p-values

```
awk '{print $8}' SETS_tree.txt

R
# read file
df <- read.table("SETS_tree.txt",
                header = TRUE)

df$p_fdr <- p.adjust(df$p.value, method = "fdr")
df$cs_fdr <- p.adjust(df$clustering_sensitive, method = "fdr")
df$cr_fdr <- p.adjust(df$clustering_robust, method = "fdr")

write.csv(df,file = "SETS_tree_fdr.csv", row.names = FALSE)
```

