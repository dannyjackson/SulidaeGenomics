#!/usr/bin/env python3
python3 - <<'PY'
import pandas as pd
import numpy as np
import glob

# Collect all divergence_thresholds.csv files
files = glob.glob("*_divergence_thresholds.csv")
dfs = [pd.read_csv(f) for f in files]
all_df = pd.concat(dfs, ignore_index=True)

# Compute summary statistics (mean, SD, 95% CI)
def summary_ci(x):
    mean = np.mean(x)
    sd = np.std(x, ddof=1)
    n = len(x)
    ci95 = 1.96 * sd / np.sqrt(n) if n > 1 else np.nan
    return pd.Series({"mean": mean, "sd": sd, "ci95": ci95, "n": n})

# apply to numeric columns only
numeric_cols = ["T50_years","T25_years","T10_years",
                "T50_generations","T25_generations","T10_generations"]
summary = all_df[numeric_cols].apply(summary_ci)

# save summary and combined data
all_df.to_csv("all_divergence_thresholds.csv", index=False)
summary.to_csv("divergence_thresholds_summary.csv")

print("\nCombined data saved to: all_divergence_thresholds.csv")
print("Summary (mean ± 95% CI):")
print(summary)
PY