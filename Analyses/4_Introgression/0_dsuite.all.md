# Dsuite all
Set working directory
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods
```
Save this as SETS.txt:
```
BFBO501 BFBO_southern
BFBO502 BFBO_southern
BFBO503 BFBO_southern
BFBO504 BFBO_GofCA
BFBO505 BFBO_southern
BRBO201 BRBO_Pacific
BRBO202 BRBO_Pacific
BRBO203 BRBO_AtlCar
BRBO205 BRBO_AtlCar
MABO302 MABO_IndoPacific
MABO304 MABO_AtlCar
MABO305 MABO_AtlCar
MABO306 MABO_IndoPacific
NABO402 NABO
NABO403 NABO
NABO404 NABO
NABO405 NABO
NABO406 NABO
PEBO601 PEBO
PEBO603 PEBO
PEBO604 PEBO
PEBO605 PEBO
PEBO606 PEBO
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup
```
Filter vcf to just autosomes and then run dsuite on the entire genome
```
sbatch 0a_dsuite_autosomes.sh
```
Compute f-branch statistics and a plot to visualize the output using scripts from Dsuite:

```

~/programs/Dsuite/Build/Dsuite Fbranch --Pb-matrix all.nwk SETS_tree.txt  > all_Fbranch.txt


module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py all_Fbranch.txt all.nwk --tree-label-size 3 --ladderize --use_distances
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
' all_Fbranch.txt > all_Fbranch.pvalues.tsv
```

Correct for false discovery rate using the p-values from the f_branch p-value statistic using R
```
R
# read file
x <- read.table("all_Fbranch.pvalues.tsv",
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
Repeat the above process with the clustering sensitive p-values

```
awk '{print $8}' SETS_tree.txt

R
# read file
df <- read.table("SETS_Dmin.txt",
                header = TRUE)

df$p_fdr <- p.adjust(df$p.value, method = "fdr")
df$cs_fdr <- p.adjust(df$clustering_sensitive, method = "fdr")
df$cr_fdr <- p.adjust(df$clustering_robust, method = "fdr")

write.csv(df,file = "SETS_Dmin_fdr.csv", row.names = FALSE)
```
