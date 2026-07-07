# Analyzing genes
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist

module load micromamba
micromamba activate r_ocelote

## Overlap across species
```
library(ggvenn)

df_MABO<-read.csv("MABO/MABO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BRBO<-read.csv("BRBO/BRBO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  MABO= df_MABO,
  BRBO = df_BRBO
)

intersect_all <- Reduce(intersect, list(df_MABO, df_BRBO))
write(intersect_all, file="BRBO_MABOpops.fixedsites.intersection.csv")

union_all <- Reduce(union, list(df_MABO, df_BRBO))
write(union_all, file="BRBO_MABOpops.fixedsites.union.csv")

pdf(file = "BRBO_MABOpops.fixedsites.intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("MABO", "BRBO"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()
# repeat with FST windows

df_MABO<-read.csv("MABOatlcar_MABOindopac.fst_50000.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BRBO<-read.csv("BRBOatlcar__BRBOpac_.fst_50000.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  MABO= df_MABO,
  BRBO = df_BRBO
)

intersect_all <- Reduce(intersect, list(df_MABO, df_BRBO))
write(intersect_all, file="MABO.BRBO.fst.50kb.intersection.csv")

union_all <- Reduce(union, list(df_MABO, df_BRBO))
write(union_all, file="MABO.BRBO.fst.50kb.union.csv")

pdf(file = "MABO.BRBO.fst.50kb.intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("MABO", "BRBO"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()



cd /Users/danjack/Documents/SulidaePaper/2025_revision/AnalysisFiles/PantherAnalysis
ls *.txt > filelist.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r BFBO_PEBO.fixedsites_in_genes.GO.Gallus
grep -o 'GO:[0-9]\{7\}' BFBO_PEBO.fixedsites_in_genes.GO.Gallus.fdr.05.txt > BFBO_PEBO.fixedsites_in_genes.GO.Gallus.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r BFBO_PEBO.fixedsites_in_genes.GO.Gallus.fdr.GO

while read -r file;
do
echo $file

sed "1s/.*/GO_Term\tBackground\tEmpirical\tExpected\tOver_Under\tfoldEnrichment\tP.value/" $file > ${file%.txt}.newhead.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r ${file%.txt}.newhead

grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.05.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO
mv plot/${file%.txt}* plot/fdr05/

grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.1.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO
mv plot/${file%.txt}* plot/fdr1/

grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.2.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO
mv plot/${file%.txt}* plot/fdr2/


file=BRBO.fst_50000.GO.Gallus.GO.txt
grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.2.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO


while read -r file;
do
file=BFPE_MANA.fixedsites.intersection.txt
awk -F'\t' -v OFS='\t' '
NR==1 {
    for (i=1; i<=NF; i++)
        if ($i == "Empirical") col = i
    print
    next
}
$col > 0
' ${file%.txt}.newhead.txt > ${file%.txt}.min1gene.txt
grep -o 'GO:[0-9]\{7\}' ${file%.txt}.min1gene.txt > ${file%.txt}.min1gene.GO.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.min1gene.GO

done < filelist.txt



BFBO_PEBO.fixedsites_in_genes.GO.Gallus.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO

sed "1s/.*/GO_Term\tBackground\tEmpirical\tExpected\tOver_Under\tfoldEnrichment\tP.value/" BFBO_PEBO.fixedsites_in_genes.GO.Gallus.txt


mkdir scatter heat tree
mv *treemap* tree 
mv *scatter.pdf scatter 
mv *heatmap.pdf heat 





cd /Users/danjack/Documents/SulidaePaper/2025_revision/AnalysisFiles/GeneLists

library(ggvenn)

df_MABO<-read.csv("MABO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BRBO<-read.csv("BRBO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_MABONABO<-read.csv("MABO_NABO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BFBOPEBO<-read.csv("BFBO_PEBO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]


gene_sets <- list(
  MABO= df_MABO,
  BRBO = df_BRBO,
  MABONABO = df_MABONABO,
  BFBOPEBO = df_BFBOPEBO
)

intersect_all <- Reduce(intersect, list(df_MABO, df_BRBO, df_MABONABO, df_BFBOPEBO))
write(intersect_all, file="All.fixedsites.intersection.csv")

pdf(file = "All.fixedsites.intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("MABO", "BRBO", "MABONABO", "BFBOPEBO"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

# repeat with just species comparisons

library(ggvenn)

df_MABONABO<-read.csv("MABO_NABO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BFBOPEBO<-read.csv("BFBO_PEBO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]


gene_sets <- list(
  MABONABO = df_MABONABO,
  BFBOPEBO = df_BFBOPEBO
)

intersect_all <- Reduce(intersect, list(df_MABONABO, df_BFBOPEBO))
write(intersect_all, file="BFPE.MANA.fixedsites.intersection.csv")

pdf(file = "BFPE.MANA.fixedsites.intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("MABONABO", "BFBOPEBO"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


file="All.fixedsites_in_genes.GO.Gallus.min1gene.GO.txt"
file="BFPE_MANA.fixedsites.intersection.txt"
file="MABO_NABO.fst_50000.GO.Gallus.txt"
file="BRBO.fixedsites_in_genes.GO.Gallus.txt"
file="MABO.fixedsites_in_genes.GO.Gallus.txt"

sed "1s/.*/GO_Term\tBackground\tEmpirical\tExpected\tOver_Under\tfoldEnrichment\tP.value/" $file > ${file%.txt}.newhead.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r ${file%.txt}.newhead

grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.05.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO
mv plot/${file%.txt}* plot/fdr05/

grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.1.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO
mv plot/${file%.txt}* plot/fdr1/

grep -o 'GO:[0-9]\{7\}' ${file%.txt}.newhead.fdr.2.txt > ${file%.txt}.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r ${file%.txt}.fdr.GO
mv plot/${file%.txt}* plot/fdr2/


# BFBO_PEBO vs MABO_NABO
# fixed sites
library(ggvenn)
df_MABONABO<-read.csv("MABO_NABO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BFBOPEBO<-read.csv("BFBO_PEBO.fixedsites_in_genes.genenames.tsv", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  MABONABO = df_MABONABO,
  BFBOPEBO = df_BFBOPEBO
)

intersect_all <- Reduce(intersect, list(df_MABONABO, df_BFBOPEBO))
write(intersect_all, file="BFPE_MANA.fixedsites.intersection.csv")

pdf(file = "BFPE_MANA.fixedsites.intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c( "MABONABO", "BFBOPEBO"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

# fst windows
df_MABONABO<-read.csv("MABO_NABO.fst_50000.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_BFBOPEBO<-read.csv("BFBO_PEBO.fst_50000.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  MABONABO = df_MABONABO,
  BFBOPEBO = df_BFBOPEBO
)

intersect_all <- Reduce(intersect, list(df_MABONABO, df_BFBOPEBO))
write(intersect_all, file="BFPE_MANA.fst_50000.intersection.csv")

pdf(file = "BFPE_MANA.fixedsites.intersection.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c( "MABONABO", "BFBOPEBO"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()



# analyze GO terms for region 2 of the neo sex
cd /Users/danjack/Documents/SulidaePaper/2025_revision/AnalysisFiles/PantherAnalysis/Gallus_GO_Lists/Chr29
grep -o 'GO:[0-9]\{7\}' Region2.txt > Region2.fdr.GO.txt

sed "1s/.*/GO_Term\tBackground\tEmpirical\tExpected\tOver_Under\tfoldEnrichment\tP.value\tQ.value/" Region2.GO.txt > Region2.newhead.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r Region2.newhead
grep -o 'GO:[0-9]\{7\}' Region2.newhead.fdr.05.txt > Region2.newhead.fdr.05.GO.txt

mkdir -p plot/fdr05 plot/fdr1 plot/fdr2
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r Region2.newhead.fdr.05.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r Region2.newhead.fdr.1.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r Region2.newhead.fdr.2.GO
file="MABO.fixedsites_in_genes.GO.Gallus.txt"
file="Region2.txt"

grep -o 'GO:[0-9]\{7\}' Region2.fdr05.txt > Region2.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r Region2.fdr.GO
