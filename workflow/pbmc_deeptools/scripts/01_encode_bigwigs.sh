#!/bin/bash
# Index the ENCODE PBMC ChIP BAMs and make bigwigs with the SAME bamCoverage
# command as the pipeline's bam_to_bw rule. Runs the 5 BAMs in parallel.
# Single-end BAMs get --extendReads 200 (paired-end use the pipeline's bare
# --extendReads, which extends to the mate).
set -u
BIN=/date/gcb/gcb_MZ/nanoCTAR_YD/MS/.snakemake/conda/67977e66f74713228faaa1470a4bdd30_/bin
REFS=/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/Refs
OUT=/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/encode_correlation/bigwigs
THREADS=8
mkdir -p "$OUT"

make_bw () {
  local bam="$1" name="$2"
  # work on a BAM in the output dir (never write into Refs); symlink if sorted.
  local work="$OUT/$name.input.bam"
  if "$BIN/samtools" view -H "$bam" | grep -q "SO:coordinate"; then
    ln -sf "$bam" "$work"
  else
    echo "[$(date +%H:%M:%S)] $name: sorting"
    "$BIN/samtools" sort -@ "$THREADS" -o "$work" "$bam"
  fi
  bam="$work"
  "$BIN/samtools" index -@ "$THREADS" "$bam"
  # paired-end? count paired reads on chr1 via the index (no head/SIGPIPE)
  local paired ext
  paired=$("$BIN/samtools" view -c -f 1 "$bam" chr1 2>/dev/null || echo 0)
  if [ "$paired" -eq 0 ]; then ext="--extendReads 200"; else ext="--extendReads"; fi
  echo "[$(date +%H:%M:%S)] $name: bamCoverage (paired/1000=$paired -> '$ext')"
  # >>> identical to pipeline bam_to_bw, plus SE fragment length when needed <<<
  "$BIN/bamCoverage" -b "$bam" -o "$OUT/$name.bw" -p "$THREADS" --minMappingQuality 5 \
    --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates $ext \
    && echo "[$(date +%H:%M:%S)] $name: DONE" || echo "[$(date +%H:%M:%S)] $name: FAILED"
}

pids=()
for bam in "$REFS"/PBMC_*.bam; do
  base=$(basename "$bam" .bam)            # PBMC_H3K27ac_ENCSR105EMQ
  make_bw "$bam" "ENCODE_${base#PBMC_}" & # -> ENCODE_H3K27ac_ENCSR105EMQ.bw
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

echo "=== ENCODE bigwigs produced ==="
ls -la "$OUT"/ENCODE_*.bw 2>/dev/null || echo "NONE"
