#!/usr/bin/env python3
"""
Lookup genes by exact product description in a GFF.

Given:
  1) A descript_file with 2 tab-separated columns:
       ID <TAB> DESCRIPTION
  2) A GFF file with CDS features and attributes including product= and gene=

This script:
  - Parses the GFF once, collecting for each product string
    the set of gene names (gene=) for CDS features.
  - For each DESCRIPTION in the descript_file, looks up that exact
    product string in the GFF-derived map.
  - Writes: ID, DESCRIPTION, comma-separated gene names (or NA)
"""

import argparse
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Match product descriptions to genes from a GFF (CDS only, exact match)."
    )
    p.add_argument(
        "descript_file",
        help="Tab-delimited file: ID <TAB> DESCRIPTION",
    )
    p.add_argument(
        "gff",
        help="GFF file (e.g. genomic.gff) containing CDS features with product= and gene= attributes.",
    )
    p.add_argument(
        "-o", "--output",
        help="Output TSV file (default: DESCRIPT_FILE + '.geneHits.tsv')",
    )
    return p.parse_args()


def parse_gff_products_to_genes(gff_path):
    """
    Parse GFF and return: dict[product_string] -> set of gene names

    Only uses lines where:
      - column 3 == 'CDS'
      - attributes contain product= and gene=
      - gene does not start with 'LOC'
    """
    product_to_genes = defaultdict(set)

    with open(gff_path, "r", encoding="utf-8") as gff:
        for line in gff:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != "CDS":
                continue

            attributes = parts[8]
            product = None
            gene = None

            # parse attributes like key=value;key2=value2;...
            for attr in attributes.split(";"):
                if not attr:
                    continue
                if "=" not in attr:
                    continue
                key, value = attr.split("=", 1)
                key = key.strip()
                value = value.strip()

                if key == "product":
                    product = value
                elif key == "gene":
                    gene = value

            if product is None or gene is None:
                continue

            if gene.startswith("LOC"):
                continue

            product_to_genes[product.lower()].add(gene)

    return product_to_genes


def main():
    args = parse_args()
    descript_file = args.descript_file
    gff_path = args.gff
    out_path = args.output or f"{descript_file}.geneHits.tsv"

    # Build map: product -> set(genes)
    sys.stderr.write(f"Parsing GFF: {gff_path}\n")
    product_to_genes = parse_gff_products_to_genes(gff_path)
    sys.stderr.write(f"Found {len(product_to_genes)} unique product descriptions with genes.\n")

    # Process descript file
    sys.stderr.write(f"Processing descriptions: {descript_file}\n")
    with open(descript_file, "r", encoding="utf-8") as infile, \
         open(out_path, "w", encoding="utf-8") as out:

        for line in infile:
            line = line.rstrip("\n")
            if not line:
                continue

            # Expect: ID<TAB>DESC
            parts = line.split("\t", 1)
            if len(parts) < 2:
                # no description; treat as NA
                id_val = parts[0]
                desc = ""
            else:
                id_val, desc = parts

            genes = product_to_genes.get(desc.lower(), None)

            if not genes:
                gene_str = "NA"
            else:
                gene_str = ",".join(sorted(genes))

            out.write(f"{id_val}\t{desc}\t{gene_str}\n")

    sys.stderr.write(f"Done. Output written to: {out_path}\n")


if __name__ == "__main__":
    main()
