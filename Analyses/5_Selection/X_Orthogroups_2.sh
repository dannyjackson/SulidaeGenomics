
###############################################################################
# STEP 1 -- Ensure that gene dictionarys are tab separated
###############################################################################


MOR_IN="bMorBas2_genes.tsv"
PHA_IN="bPhaCar2_genes.tsv"

MOR_OUT="bMorBas2_genes.norm.tsv"
PHA_OUT="bPhaCar2_genes.norm.tsv"

echo "Normalizing MOR -> $MOR_OUT"
awk 'NF>=2 {print $1 "\t" $NF}' "$MOR_IN" > "$MOR_OUT"

echo "Normalizing PHA -> $PHA_OUT"
awk 'NF>=2 {print $1 "\t" $NF}' "$PHA_IN" > "$PHA_OUT"

echo "Done step 1."

###############################################################################
# STEP 2a -- Clean Orthogroups (remove -R* and dedup Mor IDs)
###############################################################################

ORTH_DIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder/files/for_genespace/results

FILE1=${ORTH_DIR}/bMorBas2__v__bPhaCar2.tsv

OG_OUT="bMorBas2__v__bPhaCar2.clean.tsv"

echo "Cleaning $FILE1 -> $OG_OUT"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    # header
    print
    next
}
{
    og  = $1
    mb  = $2   # bMorBas2 list
    pc  = $3   # bPhaCar2 list (leave as-is)

    # split MorBas list on comma + optional spaces
    n = split(mb, a, /, */)

    delete seen
    out = ""
    for (i=1; i<=n; i++) {
        id = a[i]
        gsub(/^[ \t]+|[ \t]+$/, "", id)  # trim spaces
        sub(/-R[0-9]+$/, "", id)         # remove -R1, -R2, etc

        if (id == "") continue
        if (!(id in seen)) {
            seen[id] = 1
            out = (out=="" ? id : out ", " id)
        }
    }

    $2 = out
    print
}' "$FILE1" > "$OG_OUT"

echo "Done step 2a."


###############################################################################
# STEP 2a -- Clean Orthogroups (remove -R* and dedup Mor IDs)
###############################################################################

ORTH_DIR=/xdisk/mcnew/dannyjackson/sulidae/analyses/orthofinder/files/for_genespace/results

FILE2=${ORTH_DIR}/bPhaCar2__v__bMorBas2.tsv

OG_OUT="bPhaCar2__v__bMorBas2.clean.tsv"

echo "Cleaning $FILE2 -> $OG_OUT"

awk 'BEGIN{FS=OFS="\t"}
NR==1 { print; next }
{
    og  = $1
    pc  = $2
    mb  = $3

    n = split(mb, a, /, */)

    delete seen
    out = ""
    for (i=1; i<=n; i++) {
        id = a[i]
        gsub(/^[ \t]+|[ \t]+$/, "", id)
        gsub(/-R[0-9]+/, "", id)   # <--- changed from sub(...$) to gsub(...) 

        if (id == "") continue
        if (!(id in seen)) {
            seen[id] = 1
            out = (out=="" ? id : out ", " id)
        }
    }

    $3 = out
    print
}' "$FILE2" > "$OG_OUT"


echo "Done step 2b."


###############################################################################
# STEP 3a -- Explode OGs into Mor–protein pairs
###############################################################################

OG_CLEAN="bMorBas2__v__bPhaCar2.clean.tsv"
PAIRS_OUT="bMorBas2__v__bPhaCar2.protein_pairs.tsv"

echo "Making Mor–protein pairs -> $PAIRS_OUT"

awk 'BEGIN{FS=OFS="\t"}
NR==1 { next }  # skip header
{
    og  = $1
    mb  = $2
    pcs = $3

    # split mor list
    mcount = split(mb, ms, /, */)

    # split protein list
    pcount = split(pcs, ps, /, */)

    for (mi = 1; mi <= mcount; mi++) {
        m = ms[mi]
        gsub(/^[ \t]+|[ \t]+$/, "", m)
        if (m=="") continue

        for (pi = 1; pi <= pcount; pi++) {
            p = ps[pi]
            gsub(/^[ \t]+|[ \t]+$/, "", p)
            if (p=="") continue

            print m, p
        }
    }
}' "$OG_CLEAN" > "$PAIRS_OUT"

echo "Done step 3a."


###############################################################################
# STEP 3b -- Explode OGs into Mor–protein pairs
###############################################################################

OG_CLEAN="bPhaCar2__v__bMorBas2.clean.tsv"
PAIRS_OUT="bPhaCar2__v__bMorBas2.protein_pairs.tsv"

echo "Making Mor–protein pairs -> $PAIRS_OUT"

awk 'BEGIN{FS=OFS="\t"}
NR==1 { next }  # skip header
{
    og  = $1
    mb  = $2
    pcs = $3

    # split mor list
    mcount = split(mb, ms, /, */)

    # split protein list
    pcount = split(pcs, ps, /, */)

    for (mi = 1; mi <= mcount; mi++) {
        m = ms[mi]
        gsub(/^[ \t]+|[ \t]+$/, "", m)
        if (m=="") continue

        for (pi = 1; pi <= pcount; pi++) {
            p = ps[pi]
            gsub(/^[ \t]+|[ \t]+$/, "", p)
            if (p=="") continue

            print m, p
        }
    }
}' "$OG_CLEAN" > "$PAIRS_OUT"

echo "Done step 3b."


###############################################################################
# STEP 4a -- Aggregate per Mor + map proteins → gene names
###############################################################################

PAIRS_IN="bMorBas2__v__bPhaCar2.protein_pairs.tsv"
PHA_NORM="bPhaCar2_genes.norm.tsv"
OUT="bMorBas2_to_PhaCar_proteins_genes.tsv"

echo "Aggregating and mapping -> $OUT"

awk -v PHA="$PHA_NORM" '
BEGIN { FS=OFS="\t" }

# First, read PHA file: protein ID -> gene symbol
FNR==NR {
    prot = $1
    gene = $2
    prot2gene[prot] = gene
    next
}

# Then read Mor_protein_pairs.tsv
{
    mor = $1
    p   = $2

    # accumulate proteins (deduplicate per Mor)
    key_p = mor "|" p
    if (!(key_p in seen_p)) {
        seen_p[key_p] = 1
        if (protlist[mor] == "") protlist[mor] = p
        else                      protlist[mor] = protlist[mor] ", " p
    }

    # map protein -> gene symbol and accumulate (deduplicate per Mor)
    g = prot2gene[p]
    if (g != "") {
        key_g = mor "|" g
        if (!(key_g in seen_g)) {
            seen_g[key_g] = 1
            if (genelist[mor] == "") genelist[mor] = g
            else                      genelist[mor] = genelist[mor] ", " g
        }
    }
}

END {
    for (mor in protlist) {
        print mor, protlist[mor], genelist[mor]
    }
}
' "$PHA_NORM" "$PAIRS_IN" > "$OUT"

echo "Done step 4a."


###############################################################################
# STEP 4b -- Aggregate per Mor + map proteins → gene names
###############################################################################

PAIRS_IN="bPhaCar2__v__bMorBas2.protein_pairs.tsv"
PHA_NORM="bMorBas2_genes.norm.tsv"
OUT="bPhaCar2__v__bMorBas2_proteins_genes.tsv"

echo "Aggregating and mapping -> $OUT"

awk -v PHA="$PHA_NORM" '
BEGIN { FS=OFS="\t" }

# First, read PHA file: protein ID -> gene symbol
FNR==NR {
    prot = $1
    gene = $2
    prot2gene[prot] = gene
    next
}

# Then read Mor_protein_pairs.tsv
{
    mor = $1
    p   = $2

    # accumulate proteins (deduplicate per Mor)
    key_p = mor "|" p
    if (!(key_p in seen_p)) {
        seen_p[key_p] = 1
        if (protlist[mor] == "") protlist[mor] = p
        else                      protlist[mor] = protlist[mor] ", " p
    }

    # map protein -> gene symbol and accumulate (deduplicate per Mor)
    g = prot2gene[p]
    if (g != "") {
        key_g = mor "|" g
        if (!(key_g in seen_g)) {
            seen_g[key_g] = 1
            if (genelist[mor] == "") genelist[mor] = g
            else                      genelist[mor] = genelist[mor] ", " g
        }
    }
}

END {
    for (mor in protlist) {
        print mor, protlist[mor], genelist[mor]
    }
}
' "$PHA_NORM" "$PAIRS_IN" > "$OUT"

echo "Done step 4b."


###############################################################################
# STEP 5 -- Aggregate per Mor + map proteins → gene names, simple output version
###############################################################################


OUT="bMorBas2_to_PhaCar_proteins_genes.simple.tsv"

echo "Aggregating and mapping -> $OUT"

awk -v PHA="$PHA_NORM" '
BEGIN { FS=OFS="\t" }

# First, read PHA file: protein ID -> gene symbol
FNR==NR {
    prot = $1
    gene = $2
    prot2gene[prot] = gene
    next
}

# Then read Mor_protein_pairs.tsv
{
    mor = $1
    p   = $2

    # accumulate proteins (deduplicate per Mor)
    key_p = mor "|" p
    if (!(key_p in seen_p)) {
        seen_p[key_p] = 1
        if (protlist[mor] == "") protlist[mor] = p
        else                      protlist[mor] = protlist[mor] ", " p
    }

    # map protein -> gene symbol and accumulate (deduplicate per Mor)
    g = prot2gene[p]
    if (g != "") {
        key_g = mor "|" g
        if (!(key_g in seen_g)) {
            seen_g[key_g] = 1
            if (genelist[mor] == "") genelist[mor] = g
            else                      genelist[mor] = genelist[mor] ", " g
        }
    }
}

END {
    for (mor in protlist) {
        print mor, genelist[mor]
    }
}
' "$PHA_NORM" "$PAIRS_IN" > "$OUT"


sed -i 's/rna-//g' "$OUT"

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/genelist

# sed -i 's/ //g' bMorBas2_to_PhaCar_proteins_genes.tsv
sed -i 's/,//g' MorBas.autosome_genes.converted.txt

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff
FASTA=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
GFF2=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

grep '^LOC' MorBas.autosome_genes.converted.txt > MorBas.autosome_genes.LOC.txt

awk 'NR==FNR { keep[$1]; next }
BEGIN { FS = OFS = "\t" }

    $3=="gene" {
        split($9, a, ";");
        loc=""; desc="";
        for(i in a){
            if(a[i] ~ /^ID=gene-LOC/) loc = substr(a[i], 9);      # remove "ID=gene-"
            if(a[i] ~ /^description=/) desc = substr(a[i], 13);    # remove "description="
        }
        if(loc in keep) print loc "\t" desc
}' MorBas.autosome_genes.LOC.txt "$GFF2" > MorBas.autosome_genes.LOC_description.tsv

