#!/usr/bin/env python3
import os, re, glob, csv, sys, argparse
import numpy as np

parser = argparse.ArgumentParser(description="Extract split time quantiles from MSMC-IM fittingdetails files")
parser.add_argument("--base", required=True, help="Base directory containing analyses (e.g., /xdisk/.../msmc/divergence)")
parser.add_argument("--comparison", required=True, help="Comparison name (e.g., BFBO_PEBO)")
args = parser.parse_args()

BASE = args.base
COMPARISON = args.comparison
GEN = 8

pattern = os.path.join(BASE, COMPARISON, "MSMC-IM", "*_*.b1_1e-08.b2_1e-06.MSMC_IM.fittingdetails.txt")
files = sorted(glob.glob(pattern))
out = os.path.join(BASE, COMPARISON, "MSMC-IM", "split_times_quantiles.csv")

rx = re.compile(r"\[([0-9eE\.\+\-]+)\s*,\s*([0-9eE\.\+\-]+)\s*,\s*([0-9eE\.\+\-]+)\]")

rows = []
for f in files:
    with open(f) as fin:
        first = fin.readline()
    m = rx.search(first)
    if not m:
        print(f"[warn] Could not parse quantiles in: {f}", file=sys.stderr)
        continue
    q1, q2, q3 = map(float, m.groups())
    comp = os.path.basename(f).split(".b1_")[0]
    row = [comp, q1, q2, q3, q1*GEN, q2*GEN, q3*GEN]
    rows.append(row)

# Compute summary stats
Q1_vals = [r[1] for r in rows]
Q2_vals = [r[2] for r in rows]
Q3_vals = [r[3] for r in rows]
Q1y_vals = [r[4] for r in rows]
Q2y_vals = [r[5] for r in rows]
Q3y_vals = [r[6] for r in rows]

summary = [
    "MEAN±SD",
    f"{np.mean(Q1_vals):.2f} ± {np.std(Q1_vals, ddof=1):.2f}",
    f"{np.mean(Q2_vals):.2f} ± {np.std(Q2_vals, ddof=1):.2f}",
    f"{np.mean(Q3_vals):.2f} ± {np.std(Q3_vals, ddof=1):.2f}",
    f"{np.mean(Q1y_vals):.2f} ± {np.std(Q1y_vals, ddof=1):.2f}",
    f"{np.mean(Q2y_vals):.2f} ± {np.std(Q2y_vals, ddof=1):.2f}",
    f"{np.mean(Q3y_vals):.2f} ± {np.std(Q3y_vals, ddof=1):.2f}"
]

# Write output table
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["COMPARISON","Q1","Q2","Q3","Q1y","Q2y","Q3y"])
    w.writerows(rows)
    w.writerow([])
    w.writerow(summary)

print(f"Wrote: {out}")
print(f"Processed {len(rows)} files.")
