
#!/usr/bin/env bash

MOR="bMorBas2_genes.tsv"
PHA="bPhaCar2_genes.tsv"
FILE1_CLEAN="Orthogroups.clean.tsv"   # cleaned OG file (no -R*, deduped)

OUT="bMorBas2_to_PhaCar_proteins_genes.tsv"

# 1) Normalize MOR & PHA to tab-separated: <ID>\t<symbol>
awk 'NF>=2 {print $1 "\t" $NF}' "$MOR" > tmp.mor
awk 'NF>=2 {print $1 "\t" $NF}' "$PHA" > tmp.pha

# 2) For each Mor transcript, gather all XP IDs from FILE1_CLEAN and map to gene names
awk -v PHA="tmp.pha" -v MOR="tmp.mor" '
BEGIN { FS = OFS = "\t" }

# ------- First file: PHA (protein ID -> gene symbol) -------
FILENAME == PHA {
    prot = $1
    gene = $2
    pha_gene[prot] = gene
    next
}

# ------- Second file: MOR (list of Mor transcripts; preserve order) -------
FILENAME == MOR {
    mor_id = $1
    mor_order[++n_mor] = mor_id
    next
}

# ------- Third file: FILE1_CLEAN (Orthogroup, Mor, Pha proteins) -------
{
    # Expect: OG, MorID, XP list
    og  = $1
    mor = $2
    pcs = $3

    # there might be multiple Mor IDs in col2, but often 1
    mcount = split(mor, ms, /, */)

    for (mi = 1; mi <= mcount; mi++) {
        m = ms[mi]
        gsub(/^[ \t]+|[ \t]+$/, "", m)
        if (m == "") continue

        # split XP list
        pcount = split(pcs, aa, /, */)
        for (i = 1; i <= pcount; i++) {
            p = aa[i]
            gsub(/^[ \t]+|[ \t]+$/, "", p)
            if (p == "") continue

            # de-duplicate proteins per Mor
            k = m "|" p
            if (!(k in seen_prot)) {
                seen_prot[k] = 1
                if (prot_list[m] == "") prot_list[m] = p
                else                    prot_list[m] = prot_list[m] ", " p
            }

            # map protein -> gene symbol (from PHA) and de-duplicate
            gene = pha_gene[p]
            if (gene != "") {
                kg = m "|" gene
                if (!(kg in seen_gene)) {
                    seen_gene[kg] = 1
                    if (gene_list[m] == "") gene_list[m] = gene
                    else                    gene_list[m] = gene_list[m] ", " gene
                }
            }
        }
    }
}

END {
    # For each Mor transcript in the original MOR file order
    for (i = 1; i <= n_mor; i++) {
        m = mor_order[i]
        # Only print if it appears in at least one orthogroup
        if (prot_list[m] != "") {
            print m, prot_list[m], gene_list[m]
        }
    }
}
' tmp.pha tmp.mor "$FILE1" > "$OUT"

echo "✓ Wrote: $OUT"
