#!/usr/bin/env python3
# Usage:
#   python3 plot_rccr.py BRBO201_BRBO202
# Optional args:
#   --mu 1.913e-9           # mutation rate per site per generation
#   --g 14.11               # years per generation
#   --results-dir results   # directory for MSMC2 outputs
#   --outdir .              # where to write TSV/CSV/PNG
#   --msmc-tools ~/programs/msmc-tools   # directory containing plot_utils.py
#   --suffix .msmc2.final.txt            # input filename suffix; fallback to .final.txt if not found

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # safe for headless nodes
import matplotlib.pyplot as plt

def find_crossing(time, y, target):
    time = np.asarray(time); y = np.asarray(y)
    # find first index where we cross the target (either direction)
    idx = np.where(((y[:-1] >= target) & (y[1:] < target)) |
                   ((y[:-1] <= target) & (y[1:] > target)))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    x0, x1, y0, y1 = time[i], time[i+1], y[i], y[i+1]
    if y1 == y0:
        return x0
    return x0 + (target - y0) * (x1 - x0) / (y1 - y0)

def thresholds_from_years(df_years, g):
    # df_years has columns: time_years, rccr
    df_years = df_years.sort_values("time_years")
    t50y = find_crossing(df_years["time_years"], df_years["rccr"], 0.5)
    t25y = find_crossing(df_years["time_years"], df_years["rccr"], 0.25)
    t10y = find_crossing(df_years["time_years"], df_years["rccr"], 0.10)
    out = {
        "T50_years": t50y,
        "T25_years": t25y,
        "T10_years": t10y,
        "T50_generations": (t50y / g) if np.isfinite(t50y) else np.nan,
        "T25_generations": (t25y / g) if np.isfinite(t25y) else np.nan,
        "T10_generations": (t10y / g) if np.isfinite(t10y) else np.nan,
    }
    return out

def main():
    ap = argparse.ArgumentParser(description="Compute and plot MSMC2 rCCR + thresholds for a pair.")
    ap.add_argument("pair", help="Pair ID, e.g. BRBO201_BRBO202")
    ap.add_argument("--mu", type=float, default=1.913e-9, help="Mutation rate per site per generation")
    ap.add_argument("--g", type=float, default=14.11, help="Generation time (years)")
    ap.add_argument("--results-dir", default="results", help="Directory with MSMC2 outputs")
    ap.add_argument("--outdir", default=".", help="Output directory for TSV/CSV/PNG")
    ap.add_argument("--msmc-tools", default=os.path.expanduser("~/programs/msmc-tools"),
                    help="Path to directory containing plot_utils.py")
    ap.add_argument("--suffix", default=".msmc2.final.txt",
                    help="Filename suffix for joint MSMC2 output (fallback to .final.txt if not found)")
    args = ap.parse_args()

    # import plot_utils
    sys.path.append(args.msmc_tools)
    try:
        import plot_utils as pu
    except Exception as e:
        raise SystemExit(f"Could not import plot_utils from {args.msmc_tools}:\n{e}")

    # locate input file
    in1 = os.path.join(args.results_dir, f"{args.pair}{args.suffix}")
    if not os.path.exists(in1):
        fallback = os.path.join(args.results_dir, f"{args.pair}.final.txt")
        if os.path.exists(fallback):
            in1 = fallback
        else:
            raise SystemExit(f"Input not found: {in1} or {fallback}")

    os.makedirs(args.outdir, exist_ok=True)

    # compute rCCR using plot_utils (time is in YEARS)
    x_years, y_rccr = pu.crossCoalPlot(in1, mu=args.mu, gen=args.g)

    # build DataFrame and write TSVs
    df_years = pd.DataFrame({"time_years": x_years, "rccr": y_rccr})
    tsv_years = os.path.join(args.outdir, f"{args.pair}_rccr_years.tsv")
    df_years.to_csv(tsv_years, sep="\t", index=False)

    df_gens = df_years.assign(time_generations=df_years["time_years"] / args.g)[
        ["time_generations", "rccr"]
    ]
    tsv_gens = os.path.join(args.outdir, f"{args.pair}_rccr_generations.tsv")
    df_gens.to_csv(tsv_gens, sep="\t", index=False)

    # thresholds
    th = thresholds_from_years(df_years, args.g)
    th_row = {"PAIR": args.pair, **th, "mu": args.mu, "g": args.g}
    csv_th = os.path.join(args.outdir, f"{args.pair}_divergence_thresholds.csv")
    pd.DataFrame([th_row]).to_csv(csv_th, index=False)

    # plot PNG
    plt.figure(figsize=(7, 5))
    plt.plot(x_years, y_rccr, lw=2, color="black")
    plt.xscale("log")
    plt.xlabel("Time (years, log scale)")
    plt.ylabel("Relative cross-coalescence rate")
    plt.title(f"{args.pair} cross-coalescence")
    plt.ylim(0, 1.05)
    # optional reference lines
    for h in (0.5, 0.25, 0.10):
        plt.axhline(h, ls=":", lw=1, color="grey")
    plt.grid(True, which="both", ls="--", lw=0.5)
    plt.tight_layout()
    png = os.path.join(args.outdir, f"{args.pair}_rccr_plot.png")
    plt.savefig(png, dpi=300)

    # stdout summary
    def fmt(x):
        return f"{x:,.0f}" if np.isfinite(x) else "NA"
    print(
        f"[{args.pair}] T50={fmt(th['T50_years'])} y "
        f"({fmt(th['T50_generations'])} gen), "
        f"T25={fmt(th['T25_years'])} y "
        f"({fmt(th['T25_generations'])} gen), "
        f"T10={fmt(th['T10_years'])} y "
        f"({fmt(th['T10_generations'])} gen) | "
        f"mu={args.mu}, g={args.g}\n"
        f"Outputs:\n  {tsv_years}\n  {tsv_gens}\n  {csv_th}\n  {png}"
    )

if __name__ == "__main__":
    main()
