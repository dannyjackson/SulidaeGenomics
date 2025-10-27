#!/usr/bin/env bash
set -euo pipefail

SCAF="$1"
START="$2"
END="$3"
VCF_DIR="$4"
VCFDIR="$5"
PHYDIR="$6"
VCF2PHYLIP="$7"
ASC_BIAS="$8"
CSV="$9"
LOCK="${10}"


VCFGZ="${VCF_DIR}/${SCAF}.filtered.vcf.gz"
OUTPFX="${SCAF}.${START}.${END}"

# 1) Extract window to uncompressed VCF (stdout from vcftools)


vcftools \
    --gzvcf "$VCFGZ" \
    --chr "$SCAF" \
    --from-bp "$START" \
    --to-bp "$END" \
    --recode --stdout | bgzip > "$VCFDIR/${OUTPFX}.vcf.gz"

tabix -p vcf "$VCFDIR/${OUTPFX}.vcf.gz" 

# Convert to PHYLIP
python "$VCF2PHYLIP" -i "$VCFDIR/${OUTPFX}.vcf.gz" -o "$PHYDIR/${OUTPFX}.phy" -r

# Count sites
if [ -s "$PHYDIR/${OUTPFX}.phy" ]; then
  SITES="$(head -n 1 "$PHYDIR/${OUTPFX}.phy" | awk '{print $2}')" || SITES=0
else
  SITES=0
fi

NSNPS=0
if [ "$SITES" -eq 0 ]; then
  # ensure no stale outputs
  rm -f "$PHYDIR/${OUTPFX}.noINV.phy" \
        "$PHYDIR/${OUTPFX}.noINV.phy.felsenstein" \
        "$PHYDIR/${OUTPFX}.noINV.phy.stamatakis" || true
else
  # Ascertainment bias correction
  python "$ASC_BIAS" -p "$PHYDIR/${OUTPFX}.phy" -o "$PHYDIR/${OUTPFX}.noINV.phy"
  rm -f "$PHYDIR/${OUTPFX}.noINV.phy.felsenstein" "$PHYDIR/${OUTPFX}.noINV.phy.stamatakis" || true

  if [ -s "$PHYDIR/${OUTPFX}.noINV.phy" ]; then
    NSNPS="$(head -n 1 "$PHYDIR/${OUTPFX}.noINV.phy" | awk '{print $2}')" || NSNPS=0
  else
    NSNPS=0
  fi
fi

# Append to CSV atomically
{
  flock -x 200
  printf "%s,%s,%s,%s\n" "$SCAF" "$START" "$END" "$NSNPS" >> "$CSV"
} 200>"$LOCK"
