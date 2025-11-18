#!/usr/bin/env python3
import argparse, os, re, glob, sys
import numpy as np
import pandas as pd

def parse_fittingdetails(path):
    """
    First line example:
    The split time is estimated to be around [53678.3587, 68377.4874, 83076.6161] gens (i.e. 0.25,0.5,0.75 quantile)
    """
    with open(path, 'r') as f:
        first = f.readline().strip()
    m = re.search(r'\[([^\]]+)\]\s*gens', first)
    if not m:
        raise ValueError(f"Could not find split-time quantiles in first line of {path}")
    nums = [float(x.strip()) for x in m.group(1).split(',')]
    if len(nums) != 3:
        raise ValueError(f"Expected 3 quantiles, got: {nums}")
    q1, q2, q3 = nums
    return q1, q2, q3

def find_fitting_file(base, pair):
    """
    Looks under: {base}/{COMPARISON}/MSMC-IM/ for *MSMC_IM.fittingdetails.txt
    Accepts various prefixes (e.g., b1_1e-08.b2_1e-06) and orderings.
    """
    comp = pair.split('/')[-1]
    search_dir = os.path.join(base, "MSMC-IM")
    candidates = glob.glob(os.path.join(search_dir, f"{pair}*.MSMC_IM.fittingdetails.txt"))
    if not candidates:
        # fallback: accept any file for this comparison that ends correctly
        candidates = glob.glob(os.path.join(search_dir, f"*{comp}*.MSMC_IM.fittingdetails.txt"))
    if not candidates:
        raise FileNotFoundError(f"No MSMC_IM.fittingdetails.txt found under {search_dir}")
    # pick the shortest name (usually the most standard)
    return sorted(candidates, key=len)[0]


def maybe_find_migration_table(fitfile):
    """
    Try common siblings of the fittingdetails file for migration history.
    We look in the same directory for files that likely hold m(t):
      * *.MSMC_IM.migration.txt
      * *.MSMC_IM.m_of_t.tsv / .txt
      * *migration*.tsv / .txt
    Return path or None.
    """
    d = os.path.dirname(fitfile)
    stem = re.sub(r'\.MSMC_IM\.fittingdetails\.txt$', '', os.path.basename(fitfile))
    patterns = [
        os.path.join(d, f"{stem}.MSMC_IM.migration.txt"),
        os.path.join(d, f"{stem}.MSMC_IM.m_of_t.tsv"),
        os.path.join(d, f"{stem}.MSMC_IM.m_of_t.txt"),
        os.path.join(d, f"{stem}*.migration*.tsv"),
        os.path.join(d, f"{stem}*.migration*.txt"),
    ]
    for p in patterns:
        hits = glob.glob(p)
        if hits:
            return hits[0]
    return None

def migration_cessation_time(mig_path, g_years, threshold=1e-5, streak=3):
    """
    Heuristic: earliest time (in years) after which migration rate m(t) < threshold
    for 'streak' consecutive bins.
    Accept columns with time boundaries and a migration column; we try to auto-detect.
    """
    # Try TSV first, then whitespace
    try:
        df = pd.read_csv(mig_path, sep='\t')
    except Exception:
        df = pd.read_csv(mig_path, sep=r'\s+', engine='python')

    # Find time columns
    time_cols = [c for c in df.columns if 'time' in c.lower() or 'left' in c.lower() or 'right' in c.lower()]
    m_cols = [c for c in df.columns if re.search(r'\bm[_-]?(12|21|mig|migration)\b', c.lower()) or c.lower()=='m']
    if not time_cols or not m_cols:
        return np.nan

    # Build a time-mid vector in "generations" if clearly in coalescent units; we’ll assume generations already.
    # Many MSMC-IM exports already have time in generations; if it's bins with left/right, we take mid.
    if {'left_time_boundary','right_time_boundary'}.issubset(set(df.columns)):
        tgen = 0.5*(df['left_time_boundary'].to_numpy() + df['right_time_boundary'].to_numpy())
    elif {'t_left','t_right'}.issubset(set(df.columns)):
        tgen = 0.5*(df['t_left'].to_numpy() + df['t_right'].to_numpy())
    elif 'time' in df.columns:
        tgen = df['time'].to_numpy()
    else:
        # fallback: first numeric col
        tgen = df[time_cols[0]].to_numpy()

    mvals = df[m_cols[0]].astype(float).to_numpy()

    # find earliest index from which we have 'streak' consecutive bins below threshold
    below = (mvals < threshold).astype(int)
    # rolling sum over window 'streak'
    if len(below) < streak:
        return np.nan
    win = np.convolve(below, np.ones(streak, dtype=int), 'valid')
    idx = np.where(win == streak)[0]
    if len(idx) == 0:
        return np.nan
    i0 = idx[0]
    tgen_cess = tgen[i0]
    return tgen_cess * g_years  # to years

def main():
    p = argparse.ArgumentParser(description="Extract MSMC-IM split time and optional migration-cessation time.")
    p.add_argument("pair", help="e.g., BFBO501_PEBO601")
    p.add_argument("--base", required=True, help="Base directory containing <COMPARISON>/MSMC-IM/")
    p.add_argument("--g", type=float, default=14.11, help="Generation time in years")
    p.add_argument("--mig-threshold", type=float, default=1e-5, help="Threshold for m(t) cessation")
    p.add_argument("--mig-streak", type=int, default=3, help="Consecutive bins below threshold")
    p.add_argument("--out", default="MSMC_IM_split_summary.csv", help="Output CSV (appends if exists)")
    args = p.parse_args()

    fitfile = find_fitting_file(args.base, args.pair)
    q1, q2, q3 = parse_fittingdetails(fitfile)

    row = {
        "PAIR": args.pair,
        "Q25_gen": q1,
        "Q50_gen": q2,
        "Q75_gen": q3,
        "Q25_years": q1 * args.g,
        "Q50_years": q2 * args.g,
        "Q75_years": q3 * args.g,
        "fittingdetails_path": fitfile
    }

    migfile = maybe_find_migration_table(fitfile)
    if migfile:
        t_stop_y = migration_cessation_time(migfile, g_years=args.g,
                                            threshold=args.mig_threshold, streak=args.mig_streak)
        row["m_stop_years"] = t_stop_y
        row["m_stop_generations"] = (t_stop_y / args.g) if np.isfinite(t_stop_y) else np.nan
        row["migration_table_path"] = migfile
    else:
        row["m_stop_years"] = np.nan
        row["m_stop_generations"] = np.nan
        row["migration_table_path"] = ""

    # append or create
    df = pd.DataFrame([row])
    if os.path.exists(args.out):
        old = pd.read_csv(args.out)
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(args.out, index=False)

    # print a compact summary
    def fmt(x): 
        return "NA" if (x is None or (isinstance(x, float) and not np.isfinite(x))) else f"{x:,.0f}"
    print(
        f"[{args.pair}] MSMC-IM split (gens): Q25={fmt(q1)}, Q50={fmt(q2)}, Q75={fmt(q3)} | "
        f"(years @ g={args.g}): Q25={fmt(q1*args.g)}, Q50={fmt(q2*args.g)}, Q75={fmt(q3*args.g)}"
    )
    if migfile:
        print(f"    migration cessation ≈ {fmt(row['m_stop_years'])} y "
              f"({fmt(row['m_stop_generations'])} gen)  from {os.path.basename(migfile)}")
    print(f"→ Wrote/updated: {args.out}")

if __name__ == "__main__":
    main()
