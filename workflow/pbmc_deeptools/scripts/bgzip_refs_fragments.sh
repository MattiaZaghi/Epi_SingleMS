#!/usr/bin/env bash
# Make the GEO scATAC fragment files under Refs/ usable by Signac
# (CreateFragmentObject needs BGZF compression + a .tbi tabix index).
# The GEO downloads are STANDARD gzip and grouped by barcode, not
# coordinate-sorted, so tabix -p bed rejects them ("not BGZF form").
#
# For each *.bed.gz: decompress -> coordinate-sort -> re-compress with bgzip
# -> tabix index. Filenames are kept IDENTICAL (the markdown's crosswalk maps
# on the GSMxxxxxxx_..._fragments.bed.gz name). Original is only replaced after
# its index builds, so a failure leaves the source intact.
#
# All 6 files are processed CONCURRENTLY (16 cores / 50G RAM available); each
# gets its own sort scratch dir on the local disk, buffer, and per-file log.
set -uo pipefail

export PATH=/home/mattia/miniconda3_n/envs/nanoscope_general/bin:$PATH
REFS=/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/Refs
cd "$REFS"

convert_one() {
    local f=$1
    local tmp="/tmp/frag_sort_${f}_$$"
    mkdir -p "$tmp"
    if [ -f "$f.tbi" ] && bgzip -t "$f" 2>/dev/null; then
        echo "[$(date +%H:%M)] $f already BGZF + indexed, skipping"; rm -rf "$tmp"; return 0
    fi
    echo "[$(date +%H:%M)] $f : sort -> bgzip -> tabix ..."
    if zcat "$f" | grep -v '^#' \
         | sort -k1,1 -k2,2n -T "$tmp" -S 6G --parallel=2 \
         | bgzip -@ 2 > "$f.bgz.tmp" \
       && tabix -p bed "$f.bgz.tmp"; then
        mv "$f.bgz.tmp"     "$f"
        mv "$f.bgz.tmp.tbi" "$f.tbi"
        echo "[$(date +%H:%M)] $f DONE ($(ls -lh "$f" | awk '{print $5}') + .tbi)"
    else
        echo "[$(date +%H:%M)] $f FAILED -- original left intact" >&2
        rm -f "$f.bgz.tmp" "$f.bgz.tmp.tbi"
    fi
    rm -rf "$tmp"
}
export -f convert_one
export PATH

# fan out: all files at once
pids=()
for f in GSM*_fragments.bed.gz; do
    convert_one "$f" &
    pids+=($!)
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done

echo "[$(date +%H:%M)] ============ ALL FILES PROCESSED (rc=$rc) ============"
ls -lh "$REFS"/GSM*_fragments.bed.gz "$REFS"/GSM*_fragments.bed.gz.tbi
exit $rc