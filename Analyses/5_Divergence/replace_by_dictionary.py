#!/usr/bin/env python3
import sys

dict_file = sys.argv[1]          # e.g., Dictionary.txt
gene_file = sys.argv[2]          # e.g., MorBas.autosome_genes.txt
out_file = sys.argv[3] if len(sys.argv) > 3 else gene_file + ".replaced.tsv"

# Load dictionary
mapping = {}
with open(dict_file, encoding="utf-8") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            key, val = parts[0], parts[1]
            mapping[key] = val

# Replace in gene file
with open(gene_file, encoding="utf-8") as inp, open(out_file, "w", encoding="utf-8") as out:
    for line in inp:
        line = line.rstrip("\n")
        if not line:
            out.write("\n")
            continue
        parts = line.split("\t")
        # replace every field that matches a dictionary key
        parts = [mapping.get(x, x) for x in parts]
        out.write("\t".join(parts) + "\n")

print(f"Done. Output written to: {out_file}")
