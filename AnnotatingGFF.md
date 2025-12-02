# For each unknown gene in MorBas GFF:
1. Attach description / product to gene tree
2. Clean $2 to remove "description" and "-like"
3. Match remaining $2 value to lines in Chicken $GFF and 2F GFF, retain multiples but remove duplicates

1. Attach description / product to gene tree and 2. Clean $2 to remove "description" and "-like"
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes

DESCRIPT_FILE=/xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/MorBas.autosome_genes.unnamed.txt

awk -v FS='\t' -v OFS='\t' '{
    a = $2;                     # future col1
    b = $1;                     # future col2

    sub(/^Name=/, "", a);       # remove Name= only at start of col2
    gsub(/description=/, "", b); # remove description= anywhere
    sub(/-like$/, "", b)        # remove -like from end of col2

    print a, b
}' "$DESCRIPT_FILE" > "${DESCRIPT_FILE}.revised"


3. Match remaining $2 value to lines in various GFFs, retain multiples but remove duplicates


# Run this in /xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/human
# ~/programs/datasets download genome accession GCF_000001405.40 --include gff3,rna,cds,protein,genome,seq-report
# Run this in /xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/mouse
# ~/programs/datasets download genome accession GCF_000001635.27 --include gff3,rna,cds,protein,genome,seq-report
# Run this in /xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/zebrafish
# ~/programs/datasets download genome accession GCF_049306965.1 --include gff3,rna,cds,protein,genome,seq-report
# Run this in /xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/drosophila
# ~/programs/datasets download genome accession GCF_000001215.4 --include gff3,rna,cds,protein,genome,seq-report

GFF_Chicken=/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/GalGal/ncbi_dataset/data/GCF_000002315.5/genomic.gff
GFF_ZebraFinch=/xdisk/mcnew/dannyjackson/sulidae/analyses/synteny/deepspace/fivegenomes/referencegenomes/TaeGut/ncbi_dataset/data/GCA_048771995.1/genomic.gff
GFF_GreatCormorant=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff
GFF_Human=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/human/ncbi_dataset/data/GCF_000001405.40/genomic.gff
GFF_Mouse=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/mouse/ncbi_dataset/data/GCF_000001635.27/genomic.gff
GFF_Fly=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/drosophila/ncbi_dataset/data/GCF_000001215.4/genomic.gff

chmod +x match_descriptions_to_genes.py

# Chicken
./match_descriptions_to_genes.py \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/MorBas.autosome_genes.unnamed.txt.revised \
    $GFF_Chicken \
    -o MorBas.autosome_genes.unnamed.gene_matches.Chicken.txt

# Great Cormorant
./match_descriptions_to_genes.py \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/MorBas.autosome_genes.unnamed.txt.revised \
    $GFF_GreatCormorant \
    -o MorBas.autosome_genes.unnamed.gene_matches.GreatCormorant.txt

# Human 
./match_descriptions_to_genes.py \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/MorBas.autosome_genes.unnamed.txt.revised \
    $GFF_Human \
    -o MorBas.autosome_genes.unnamed.gene_matches.Human.txt

# Mouse
./match_descriptions_to_genes.py \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/MorBas.autosome_genes.unnamed.txt.revised \
    $GFF_Mouse \
    -o MorBas.autosome_genes.unnamed.gene_matches.Mouse.txt

# Drosophila
./match_descriptions_to_genes.py \
    /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/MorBas.autosome_genes.unnamed.txt.revised \
    $GFF_Fly \
    -o MorBas.autosome_genes.unnamed.gene_matches.Fly.txt

# Combine output
python3 combine_by_first_col.py \
  MorBas.autosome_genes.unnamed.gene_matches \
  -d /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes \
  -o MorBas.autosome_genes.unnamed.gene_matches.combined.tsv

python3 consensus_partition.py MorBas.autosome_genes.unnamed.gene_matches.combined.tsv



# First, subset this into two separate groups:
awk -F'\t' '
NR==1 {
    header = $0
    print header > "with_commas.tsv"
    print header > "no_commas.tsv"
    next
}
{
    if ($7 ~ /,/)          # column 7 contains a comma
        print >> "MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_TERM.internalconflict.tsv"
    else
        print >> "MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_TERM.singlematch.tsv"
}
' MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_TERM.tsv


FILE=MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_TERM.internalconflict.tsv
ORTH1=/xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder/files/for_genespace/orthofinder/Results_Nov18/Orthologues/bMorBas2.tsv  
ORTH2=/xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder/files/for_genespace/orthofinder/Results_Nov18/Orthologues/bPhaCar2.tsv 




# Resolve with Great Cormorant values
awk -F'\t' -v OFS='\t' '
NR==1 {
    # find column indices by name (avoids hard-coding)
    for (i=1; i<=NF; i++) {
        if ($i == "GreatCormorant.txt") c_value = i
        else if ($i == "consensus") c_cons = i
    }
    print; next
}
{
    if ($c_cons == "CONFLICT") {
        if ($c_value ~ /,/) {
            # leave as CONFLICT if there is a comma in GC column
        }
        else if ($c_value != "NA") {
            if ($c_value ~ /,/)
                $c_cons = "NA"
            else
                $c_cons = $c_value
        }
        # else leave as CONFLICT
    }
    print
}
' MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.tsv \
> MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.GreatCormorant_Resolved.tsv

# Resolve with Chicken values
awk -F'\t' -v OFS='\t' '
NR==1 {
    # find column indices by name (avoids hard-coding col numbers)
    for (i=1; i<=NF; i++) {
        if ($i == "Chicken.txt") c_value = i
        else if ($i == "consensus") c_cons = i
    }
    print; next
}
{
    if ($c_cons == "CONFLICT") {
        if ($c_value != "NA") {
            if ($c_value ~ /,/)
                $c_cons = "NA"
            else
                $c_cons = $c_value
        }
        # else leave as CONFLICT
    }
    print
}
' MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.GreatCormorant_Resolved.tsv \
> MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.GreatCormorant_Chicken_Resolved.tsv

# Resolve the following by hand:
# Rule: default to more general gene name
bMorBas2_egapxtmp_004691        NA      NA      NA      MT2A    MT2     CONFLICT # MT2
bMorBas2_egapxtmp_004782        NA      NA      NA      LTB4R   LTB4R1  CONFLICT # LTB4R
bMorBas2_egapxtmp_005229        NA      TBC1D7  NA      TBC1D7-LOC100130357     NA      CONFLICT # TBC1D7
bMorBas2_egapxtmp_006010        NA      NA      NA      C10ORF88        2310057M21RIK   CONFLICT # C10ORF88
bMorBas2_egapxtmp_012353        NA      NA      NA      RNF113A RNF113A1        CONFLICT # RNF113A

awk -F'\t' -v OFS='\t' '
BEGIN {
    # ID → replacement value
    repl["bMorBas2_egapxtmp_004691"] = "MT2"
    repl["bMorBas2_egapxtmp_004782"] = "LTB4R"
    repl["bMorBas2_egapxtmp_005229"] = "TBC1D7"
    repl["bMorBas2_egapxtmp_006010"] = "C10ORF88"
    repl["bMorBas2_egapxtmp_012353"] = "RNF113A"
}
NR==1 {
    # find consensus column index
    for (i=1; i<=NF; i++)
        if ($i == "consensus") c_cons = i
    print; next
}
{
    id = $1
    if ($c_cons == "CONFLICT" && id in repl)
        $c_cons = repl[id]
    print
}
' MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.GreatCormorant_Chicken_Resolved.tsv > MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.GreatCormorant_Chicken_Manual_Resolved.tsv


> Dictionary.tsv

awk -F'\t' -v OFS='\t' '{print $1, $7}' MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_TERM.tsv | grep -v ',' >> Dictionary.tsv


awk -F'\t' -v OFS='\t' '{print $1, $7}' MorBas.autosome_genes.unnamed.gene_matches.combined.tsv.consensus_CONFLICT.GreatCormorant_Chicken_Manual_Resolved.tsv  | tail -n 1 >> Dictionary.tsv

# Use dictionary to convert genes
MorBas.autosome_genes.txt
python3 replace_by_dictionary.py Dictionary.tsv MorBas.autosome_genes.txt

wc -l MorBas.autosome_genes.txt.replaced.tsv

sort -u MorBas.autosome_genes.txt.replaced.tsv > MorBas.autosome_genes.replaced.backgroundlistfinal.tsv

####################################################################################

awk -F'\t' '
{
    OFS = "\t";
    n = split($9, a, ";");
    id = ""; name = "";

    for (i=1; i<=n; i++) {
        if (a[i] ~ /^ID=/)   id = a[i];
        if (a[i] ~ /^Name=/) name = a[i];
    }

    sub(/^ID=/, "", id);
    sub(/^Name=/, "", name);

    # Keep only rows where ID begins with "gene"
    # Skip everything else (ID=rna-..., ID=transcript..., etc.)
    if (id !~ /^gene/)
        next;

    print name;
}
' BFBO_PEBO.fst_50000.genelist.txt | sort -u > BFBO_PEBO.fst_50000.genenames.2.txt

python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/replace_by_dictionary.py /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/Dictionary.tsv BFBO_PEBO.fst_50000.genenames.2.txt


####################################################################################

awk -F'\t' '
{
    OFS = "\t";
    n = split($9, a, ";");
    id = ""; name = "";

    for (i=1; i<=n; i++) {
        if (a[i] ~ /^ID=/)   id = a[i];
        if (a[i] ~ /^Name=/) name = a[i];
    }

    sub(/^ID=/, "", id);
    sub(/^Name=/, "", name);

    # Keep only rows where ID begins with "gene"
    # Skip everything else (ID=rna-..., ID=transcript..., etc.)
    if (id !~ /^gene/)
        next;

    print name;
}
' MABO_NABO.fst_50000.genelist.txt | sort -u > MABO_NABO.fst_50000.genenames.2.txt

python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/replace_by_dictionary.py /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/Dictionary.tsv MABO_NABO.fst_50000.genenames.2.txt

####################################################################################

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/BRBO_4taxa/analyses/f_dM/BRBO_Pacific_BRBO_AtlCar_4taxa/50_25
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/dsuite/all/analyses/f_dM/BFBO_GofCA_BFBO_southern_PEBO/50_25


awk '{print $4}' gene_window_overlaps.10kb.tsv > genenames.txt

python3 /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/replace_by_dictionary.py /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist/backgroundgenes/Dictionary.tsv genenames.txt

