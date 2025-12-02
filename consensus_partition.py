#!/usr/bin/env python3
"""
Given a combined TSV (ID + 1+ data columns), compute a consensus column and
write four files:

1) <input>.withConsensus.tsv                -> full table
2) <input>.consensus_NA.tsv                 -> only rows where consensus == NA
3) <input>.consensus_CONFLICT.tsv           -> only rows where consensus == CONFLICT
4) <input>.consensus_TERM.tsv               -> rows where consensus is a real term
                                               (not NA and not CONFLICT)
"""

import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Add consensus column and partition rows.")
    p.add_argument("input", help="Combined TSV file (first column = ID)")
    p.add_argument("-o", "--output", help="Prefix for output files (default: input filename)")
    return p.parse_args()


def compute_consensus(values, na_str="NA"):
    non_na = [v for v in values if v != na_str]
    if not non_na:
        return na_str
    uniq = set(non_na)
    if len(uniq) == 1:
        return non_na[0]
    return "CONFLICT"


def main():
    args = parse_args()
    in_path = args.input
    prefix = args.output or in_path

    # Output paths
    out_all      = f"{prefix}.withConsensus.tsv"
    out_na       = f"{prefix}.consensus_NA.tsv"
    out_conflict = f"{prefix}.consensus_CONFLICT.tsv"
    out_term     = f"{prefix}.consensus_TERM.tsv"

    with open(in_path, "r", encoding="utf-8") as inp, \
         open(out_all, "w", encoding="utf-8") as f_all, \
         open(out_na, "w", encoding="utf-8") as f_na, \
         open(out_conflict, "w", encoding="utf-8") as f_conf, \
         open(out_term, "w", encoding="utf-8") as f_term:

        header = inp.readline().rstrip("\n").split("\t")
        if len(header) < 2:
            sys.stderr.write("Error: expected ID + data columns.\n")
            sys.exit(1)

        # Write headers
        full_header = header + ["consensus"]
        f_all.write("\t".join(full_header) + "\n")
        f_na.write("\t".join(full_header) + "\n")
        f_conf.write("\t".join(full_header) + "\n")
        f_term.write("\t".join(full_header) + "\n")

        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            parts += ["NA"] * max(0, len(header) - len(parts))   # pad if needed

            id_val = parts[0]
            vals = parts[1:]

            consensus = compute_consensus(vals, na_str="NA")
            row = parts + [consensus]
            out_line = "\t".join(row) + "\n"

            # write to all files
            f_all.write(out_line)

            if consensus == "NA":
                f_na.write(out_line)
            elif consensus == "CONFLICT":
                f_conf.write(out_line)
            else:
                f_term.write(out_line)

    sys.stderr.write(
        "Done.\n"
        f"  Full file:         {out_all}\n"
        f"  NA-only:           {out_na}\n"
        f"  CONFLICT-only:     {out_conflict}\n"
        f"  Consensus matches: {out_term}\n"
    )


if __name__ == "__main__":
    main()
