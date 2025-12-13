#!/usr/bin/env bash
#SBATCH --job-name=contigdepth
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/contigdepth.%A_%a.out
#SBATCH --mail-type=ALL
set -euo pipefail

module load parallel

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt
BAMSTATS=/xdisk/mcnew/dannyjackson/sulidae/datafiles/bamstats
CONTIGS=/xdisk/mcnew/dannyjackson/sulidae/referencelists/CONTIGS.txt
OUTDIR="$BAMSTATS/mean_by_scaffold"
JOBS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUTDIR"

# Write awk logic to a temp file to avoid quoting issues
AWK_PROG=$(mktemp)
cat > "$AWK_PROG" <<'AWK'
BEGIN {
  # read allowed contigs
  while ((getline c < list) > 0) keep[c]=1
  close(list)
}
($1 in keep) { sum[$1]+=$3; n[$1]++ }
END {
  tsv = outdir "/" sample "_means.tsv"
  for (sc in keep) {
    mean = (n[sc] ? sum[sc]/n[sc] : 0)
    # per-sample TSV (sample \t contig \t mean)
    print sample, sc, mean > tsv
    # optional single-value files (comment out to avoid many small files)
    fn = outdir "/" sample "_" sc "_mean_coverage.txt"
    printf("%.6f\n", mean) > fn
  }
}
AWK

export BAMSTATS CONTIGS OUTDIR AWK_PROG
LC_ALL=C parallel --bar --jobs "$JOBS" \
  'awk -v list="$CONTIGS" -v OFS="\t" -v sample={1} -v outdir="$OUTDIR" -v FS="[ \t]+" -f "$AWK_PROG" "$BAMSTATS"/{1}_depthstats.txt' \
  :::: "$LIST"

# Clean up (optional)
rm -f "$AWK_PROG"
