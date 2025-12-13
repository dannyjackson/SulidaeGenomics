#!/usr/bin/env bash
# Parallel mean/sd over column 3 of *_depthstats.txt files
# Output: depthstats.tsv with columns: bird<TAB>depth<TAB>stdev

set -euo pipefail

module load parallel

# ---- Config (edit if needed) ----
LIST=/xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt
STATDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats/sexing_by_depth
OUT=${STATDIR}/depthstats.parallel.tsv

# Use all allocated CPUs if under SLURM, otherwise all local cores
JOBS="${SLURM_CPUS_ON_NODE:-$(nproc)}"

cd "$STATDIR"

# Header
printf "bird\tdepth\tstdev\n" > "$OUT"

# Speed tips: C locale for faster awk; guard against small negative var from FP errors
export LC_ALL=C

# Run in parallel; keep output in the same order as LIST (-k)
parallel -k --jobs "$JOBS" --no-run-if-empty '
  bird={}
  f="'"$STATDIR"'/${bird}_depthstats.txt"

  if [[ -s "$f" ]]; then
    awk -v bird="$bird" '"'"'
      { sum += $3; sumsq += $3*$3 }
      END {
        n = NR
        if (n > 0) {
          mean = sum / n
          var  = (sumsq / n) - (mean * mean)
          if (var < 0) var = 0  # FP guard
          sd = sqrt(var)
          printf "%s\t%.6f\t%.6f\n", bird, mean, sd
        } else {
          printf "%s\tNA\tNA\n", bird
        }
      }
    '"'"' "$f"
  else
    # Missing or empty file → NA row (still emit a line to keep order)
    printf "%s\tNA\tNA\n" "$bird"
  fi
' :::: "$LIST" >> "$OUT"

echo "Wrote $OUT"