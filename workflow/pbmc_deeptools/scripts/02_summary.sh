#!/bin/bash
# Build a signal matrix over genome-wide 10kb bins across all tracks:
#   ENCODE bigwigs (just made) + scPBMC pseudobulk bigwigs (pipeline bam_to_bw).
# H3K27ac + H3K27me3 only (ENCODE has no ATAC). All hg38.
set -u
BIN=/date/gcb/gcb_MZ/nanoCTAR_YD/MS/.snakemake/conda/67977e66f74713228faaa1470a4bdd30_/bin
BASE=/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS
OUT=$BASE/encode_correlation
BW=$OUT/bigwigs
THREADS=16

bws=(); labels=()
# --- ENCODE tracks ---
for f in "$BW"/ENCODE_H3K27ac_*.bw "$BW"/ENCODE_H3K27me3_*.bw; do
  [ -s "$f" ] || continue
  bws+=("$f"); labels+=("$(basename "$f" .bw)")
done
# --- scPBMC pseudobulk tracks (existing pipeline bigwigs) ---
for s in PBMC_CT PBMC_DynaTag; do
  for m in H3K27ac_CCTATCCT H3K27me3_ATAGAGGC; do
    f="$BASE/$s/$m/bigwig/all_reads.bw"
    [ -s "$f" ] || { echo "WARN: missing $f"; continue; }
    bws+=("$f"); labels+=("scPBMC_${s}_${m%%_*}")
  done
done

echo "=== ${#bws[@]} tracks ==="
printf '  %s\n' "${labels[@]}"

"$BIN/multiBigwigSummary" bins \
  -b "${bws[@]}" --labels "${labels[@]}" \
  --binSize 10000 -p "$THREADS" \
  -o "$OUT/summary_10kb.npz" \
  --outRawCounts "$OUT/counts_10kb.tab"

echo "=== wrote $OUT/summary_10kb.npz + counts_10kb.tab ==="
wc -l "$OUT/counts_10kb.tab"
