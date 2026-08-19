#!/usr/bin/env bash
# Export per-cell-type bigwigs for MS H3K27ac (grouped by predicted.id_Human_MS)
# and MS H3K27me3 (same label transferred to matched cells). Run inside tmux so
# it survives logout:
#   tmux new -s bigwig
#   bash /home/mattia/Multi_nanoCTRNA/workflow/pbmc_deeptools/scripts/launch_MS_HumanMS_bigwig.sh
#   # detach: Ctrl-b then d   |   reattach: tmux attach -t bigwig
set -euo pipefail

source /home/mattia/miniconda3_n/etc/profile.d/conda.sh
conda activate nanoctarna-analysis
cd /home/mattia/Multi_nanoCTRNA

RDS=/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/first_eval/objects/merged_MS.rds
OREF=/date/gcb/gcb_MZ/multiNanoCT/Analysis/MS_H3K27ac_HumanMS_bigwig_rc
OQ=/date/gcb/gcb_MZ/multiNanoCT/Analysis/MS_H3K27me3_HumanMS_bigwig_rc
LOG=/date/gcb/gcb_MZ/multiNanoCT/Analysis/MS_HumanMS_bigwig_rc.log

mkdir -p "$(dirname "$OREF")"
echo "[$(date +%H:%M)] starting export -> log: $LOG"

Rscript workflow/pbmc_deeptools/scripts/run_ExportGroupBW_matched.R \
  --rds "$RDS" --refmark H3K27ac --qmark H3K27me3 \
  --group predicted.id_Human_MS \
  --outref "$OREF" --outq "$OQ" \
  --tilesize 100 --cutoff 4 --mincells 5 2>&1 | tee "$LOG"

echo "[$(date +%H:%M)] cleaning intermediate beds (/date is tight on space)..."
rm -f "$OREF"/*.bed "$OQ"/*.bed

echo "H3K27ac  bigwigs: $(ls "$OREF"/*.bw 2>/dev/null | wc -l)/9  -> $OREF"
echo "H3K27me3 bigwigs: $(ls "$OQ"/*.bw   2>/dev/null | wc -l)/9  -> $OQ"
echo "[$(date +%H:%M)] ALL DONE"
