#!/usr/bin/env python3
"""
Combine files with a common base name, joining by column 1.

Example:
  Base name: MorBas.autosome_genes.unnamed.gene_matches

  Files:
    MorBas.autosome_genes.unnamed.gene_matches.Chicken
    MorBas.autosome_genes.unnamed.gene_matches.Fly

Output (TSV):
  ID    Chicken     Fly
  id1   value_from_Chicken   value_from_Fly
  ...

Notes:
  - All input files are assumed to be tab-separated.
  - Column 1 is the key.
  - The value taken from each file is the LAST column.
  - Output header is "ID" + one column per file, where the column
    header is the suffix after the base name (e.g. Chicken, Fly).
  - Missing values are written as "NA".
"""

import argparse
import glob
import os
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Combine tab-separated files with a common base name, joining by column 1."
    )
    p.add_argument(
        "base",
        help="Base filename, e.g. 'MorBas.autosome_genes.unnamed.gene_matches'",
    )
    p.add_argument(
        "-o", "--output",
        help="Output TSV file (default: BASE + '.combined.tsv')",
    )
    p.add_argument(
        "-d", "--directory",
        default=".",
        help="Directory to search for files (default: current directory).",
    )
    return p.parse_args()


def find_files(base, directory):
    pattern = os.path.join(directory, base + ".*")
    files = sorted(
        f for f in glob.glob(pattern)
        if os.path.isfile(f)
    )
    return files


def main():
    args = parse_args()
    base = args.base
    out_path = args.output or f"{base}.combined.tsv"
    directory = args.directory

    files = find_files(base, directory)
    # Avoid accidentally including the output file if it matches the pattern
    files = [f for f in files if os.path.abspath(f) != os.path.abspath(out_path)]

    if not files:
        sys.stderr.write(f"No files found matching '{base}.*' in {directory}\n")
        sys.exit(1)

    sys.stderr.write("Found files:\n")
    for f in files:
        sys.stderr.write(f"  {f}\n")

    # key -> {label: value}
    data = defaultdict(dict)
    labels = []

    for path in files:
        basename = os.path.basename(path)
        # Label: part after "base." in the filename
        if not basename.startswith(base + "."):
            # Shouldn't happen if glob matched correctly, but be safe
            continue
        label = basename[len(base) + 1 :]  # everything after base + "."
        labels.append(label)

        sys.stderr.write(f"Parsing {path} as label '{label}'\n")

        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if not parts:
                    continue

                key = parts[0]
                value = parts[-1]  # last column as the value
                value = value.upper()  # uppercase everything except column 1
                data[key][label] = value

    # Sort keys for stable output
    all_keys = sorted(data.keys())

    # Write output
    with open(out_path, "w", encoding="utf-8") as out:
        # Header
        out.write("ID\t" + "\t".join(labels) + "\n")

        for key in all_keys:
            row_vals = [key]
            for label in labels:
                val = data[key].get(label, "NA")
                row_vals.append(val)
            out.write("\t".join(row_vals) + "\n")

    sys.stderr.write(f"Done. Wrote combined table to: {out_path}\n")


if __name__ == "__main__":
    main()
