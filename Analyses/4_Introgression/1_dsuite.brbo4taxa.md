# Dsuite brbo vs 4 taxa

Set working directory
```
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO_4taxa_likelihoods
```
Save this as SETS.txt:
```
BFBO501 4taxa
BFBO502 4taxa
BFBO503 4taxa
BFBO504 4taxa
BFBO505 4taxa
BRBO201 BRBO_Pacific
BRBO202 BRBO_Pacific
BRBO203 BRBO_AtlCar
BRBO205 BRBO_AtlCar
MABO302 4taxa
MABO304 4taxa
MABO305 4taxa
MABO306 4taxa
NABO402 4taxa
NABO403 4taxa
NABO404 4taxa
NABO405 4taxa
NABO406 4taxa
PEBO601 4taxa
PEBO603 4taxa
PEBO604 4taxa
PEBO605 4taxa
PEBO606 4taxa
RFBO101 Outgroup
RFBO102 Outgroup
RFBO103 Outgroup
RFBO104 Outgroup
RFBO105 Outgroup
RFBO106 Outgroup
```
Create tree file, and clean SETS file:
```
tr ' ' '\t' < SETS.txt > SETS.tmp
mv SETS.tmp SETS.txt
echo "(((BRBO_AtlCar,BRBO_Pacific),4taxa),Outgroup);" > BRBO_4taxa.nwk
```
Run dsuite on the filtered genome generated in 0_dsuite.all.md
```
sbatch 1a_BRBO4taxa_dsuite.sh
```
Create fbranch plot, and then redo fbranch script to generate p-values. Note: dtools.py does not run on a matrix generated with --Pb-matrix, so it must be run twice.
```
~/programs/Dsuite/Build/Dsuite Fbranch BRBO_4taxa.nwk SETS_tree.txt > BRBO4taxa_Fbranch.txt

module load micromamba
# micromamba create -n pyplots -c conda-forge python=3.10 matplotlib
micromamba activate pyplots
# micromamba install -c conda-forge pandas numpy matplotlib seaborn
python3 ~/programs/Dsuite/utils/dtools.py BRBO_4taxa_Fbranch.txt BRBO_4taxa.nwk --color-cutoff 0.3 --tree-label-size 3 --ladderize --use_distances

~/programs/Dsuite/Build/Dsuite Fbranch --Pb-matrix BRBO_4taxa.nwk SETS_tree.txt > BRBO4taxa_Fbranch.txt
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
' BRBO4taxa_Fbranch.txt > BRBO4taxa_Fbranch.pvalues.tsv
```

Correct for false discovery rate using the p-values from the f_branch p-value statistic using R
```
micromamba activate r_ocelote
R
# read file
x <- read.table("BRBO4taxa_Fbranch.pvalues.tsv",
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
1

Repeat the above process with the clustering sensitive p-values

```

R
# read file
df <- read.table("SETS_tree.txt",
                header = TRUE)

df$p_fdr <- p.adjust(df$p.value, method = "fdr")
df$cs_fdr <- p.adjust(df$clustering_sensitive, method = "fdr")
df$cr_fdr <- p.adjust(df$clustering_robust, method = "fdr")

write.csv(df,file = "SETS_tree_fdr.csv", row.names = FALSE)
```




Use Dinvestigate to infer regions of introgression:
```
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all_likelihoods/Sula_MorusBassanus.qualitysort_filtered.autosomes.vcf.gz

echo -e 'BRBO_Pacific\tBRBO_AtlCar\t4taxa' > test_trios.txt
~/programs/Dsuite/Build/Dsuite Dinvestigate -g -w 50,25 $VCF SETS.txt test_trios.txt

gunzip /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf.gz
```

Plot the data in sliding windows and extract sits in the top 0.9999th percentile of f_d, f_dM, and d_f.

```
Rscript ../localFstats_manhattan.R \
  BRBO_Pacific_BRBO_AtlCar_4taxa_localFstats__50_25.txt \
  /xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt \
  0.9999
```

Prior gene analyses (keep to guide scripts for new round)
```
chromo	chr	position	BPcum	windowStart	windowEnd	D	f_d	f_dM	d_f	cutoff_used
4	CM062570.1	76150455.5	499231318.5	76150228	76150683	0.14719	0.575061	0.55087	0.067784	0.5003591676001057
13	CM062579.1	26413271	921799107	26412558	26413984	0.249313	0.625807	0.607195	0.105941	0.5003591676001057
22	CM062588.1	930854.5	1065245245	920676	941033	0.524932	0.620804	0.519498	0.253805	0.5003591676001057
22	CM062588.1	940591.5	1065254982	939221	941962	0.475157	0.626907	0.522118	0.229384	0.5003591676001057
23	CM062589.1	292582	1073795662.5	292358	292806	0.81585	0.533434	0.516381	0.587031	0.5003591676001057


GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff
############################################
# Are there genes within these signals? 
############################################
# first signal of introgression
awk '$1 == "CM062570.1" && $4 >= 76150228 && $5 <= 76150683' \
  $GFF | grep 'ID=gene'

# second signal of introgression
awk '$1 == "CM062579.1" && $4 >= 26412558 && $5 <= 26413984' \
  $GFF | grep 'ID=gene'

# third signal of introgression
awk '$1 == "CM062588.1" && $4 >= 920676 && $5 <= 941033' \
  $GFF | grep 'ID=gene'
# CM062588.1      Liftoff gene    932580  932874  .       +       .       ID=gene-LOC135317588;Dbxref=GeneID:135317588;Name=LOC135317588;description=uncharacterized LOC135317588;gbkey=Gene;gene=LOC135317588;gene_biotype=lncRNA;coverage=0.540;sequence_ID=0.470;extra_copy_number=0;copy_num_ID=gene-LOC135317588_0;low_identity=True

# fourth signal of introgression
awk '$1 == "CM062588.1" && $4 >= 939221 && $5 <= 941962' \
  $GFF | grep 'ID=gene'

# fifth signal of introgression
awk '$1 == "CM062589.1" && $4 >= 292358 && $5 <= 292806' \
  $GFF | grep 'ID=gene'```

  Potential analyses: 
I was thinking it might make sense to look at which SNPs are fixed between brown booby populations and to then 
```
VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/Sula_MorusBassanus.qualitysort_filtered_mind2.vcf

ruby get_fixed_site_gts.rb $VCF BRBO_4taxa.fixed.txt BRBO201,BRBO202 BFBO501,BFBO502,BFBO503,BFBO504,BFBO505,MABO302,MABO304,MABO305,MABO306,NABO402,NABO403,NABO404,NABO405,NABO406,PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_4taxa.fixed.txt BRBO_4taxa.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_BFBO.fixed.txt BRBO201,BRBO202 BFBO501,BFBO502,BFBO503,BFBO504,BFBO505 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_BFBO.fixed.txt BRBO_BFBO.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_MABO.fixed.txt BRBO201,BRBO202 MABO302,MABO304,MABO305,MABO306 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_MABO.fixed.txt BRBO_MABO.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_NABO.fixed.txt BRBO201,BRBO202 NABO402,NABO403,NABO404,NABO405,NABO406 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_NABO.fixed.txt BRBO_NABO.fixed.svg 1.0 1000

ruby get_fixed_site_gts.rb $VCF BRBO_PEBO.fixed.txt BRBO201,BRBO202 PEBO601,PEBO603,PEBO604,PEBO605,PEBO606 BRBO203,BRBO205 1.0
ruby plot_fixed_site_gts.rb BRBO_PEBO.fixed.txt BRBO_PEBO.fixed.svg 1.0 1000
```