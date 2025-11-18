# chromosome 29 gene content analysis
# analyze genes from PhaCro in syntenic chromosome 28; NC_087540.1
GFF=/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/genomeRepo/PhaCar/PhaCar.gff

mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/chr29_geneanalysis
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/chr29_geneanalysis

mkdir -p /xdisk/mcnew/dannyjackson/sulidae/analyses/chr29_geneanalysis/PhaCar
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/chr29_geneanalysis/PhaCar

# make list of background genes
grep 'ID\=gene' $GFF | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > background.PhaCar.genes.txt

# make list of genes on chr Z2
grep 'NC_087540.1' $GFF > NC_087540.1.gff

grep 'ID\=gene' NC_087540.1.gff | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > NC_087540.1.PhaCar.genes.txt

GFF_RE=/xdisk/mcnew/dannyjackson/sulidae/analyses/repeats/vertebrata/ChrsOfInterest.fa.out.gff
GFF_MB=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
