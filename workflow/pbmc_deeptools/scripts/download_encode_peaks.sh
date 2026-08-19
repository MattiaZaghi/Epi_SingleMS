#!/usr/bin/env bash
# Download ENCODE peak files - using direct file IDs
set -euo pipefail

OUT_DIR="/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/peaks_analysis/encode_peaks"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "Downloading ENCODE peak files to $OUT_DIR"
echo "Note: If downloads fail, download manually from:"
echo "  https://www.encodeproject.org/experiments/ENCSR156XNC/ (H3K27ac)"
echo "  https://www.encodeproject.org/experiments/ENCSR615HXA/ (H3K27ac)"
echo "  https://www.encodeproject.org/experiments/ENCSR442BHP/ (H3K27me3)"
echo "  https://www.encodeproject.org/experiments/ENCSR553XBX/ (H3K27me3)"
echo ""

# Try with increased timeout and verbose error output
download_file() {
    local url=$1
    local output=$2
    echo "Downloading: $output"
    wget --timeout=30 --tries=3 -O "$output" "$url" 2>&1 || {
        echo "  Download failed for: $url"
        return 1
    }
}

# Known file IDs from ENCODE (this would need to be updated if files change)
# Format: ENCODE_<mark>_<accession>.peaktype

# Try using ENCODE portal API to resolve file URLs
echo "Attempting to fetch file URLs from ENCODE..."

# H3K27ac
for acc in ENCSR156XNC ENCSR615HXA; do
    output="ENCODE_H3K27ac_${acc}.narrowPeak.gz"
    if [ -f "$output" ] && [ -s "$output" ]; then
        echo "✓ Already exists: $output"
        continue
    fi

    # Query ENCODE API for this experiment
    echo "Querying ENCODE for $acc..."
    python3 << EOFPYTHON
import urllib.request, json, sys

try:
    with urllib.request.urlopen(f"https://www.encodeproject.org/experiments/${sys.argv[1]}/?format=json", timeout=10) as r:
        data = json.loads(r.read())
        for f in data.get("files", []):
            if f.get("file_type") == "narrowPeak" and f.get("status") == "released":
                url = "https://www.encodeproject.org" + f["href"]
                print(url)
                break
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
EOFPYTHON acc 2>/dev/null | while read url; do
        if [ -n "$url" ]; then
            download_file "$url" "$output" || rm -f "$output"
        fi
    done
done

# H3K27me3
for acc in ENCSR442BHP ENCSR553XBX; do
    output="ENCODE_H3K27me3_${acc}.broadPeak.gz"
    if [ -f "$output" ] && [ -s "$output" ]; then
        echo "✓ Already exists: $output"
        continue
    fi

    echo "Querying ENCODE for $acc..."
    python3 << EOFPYTHON
import urllib.request, json, sys

try:
    with urllib.request.urlopen(f"https://www.encodeproject.org/experiments/${sys.argv[1]}/?format=json", timeout=10) as r:
        data = json.loads(r.read())
        for f in data.get("files", []):
            if f.get("file_type") == "broadPeak" and f.get("status") == "released":
                url = "https://www.encodeproject.org" + f["href"]
                print(url)
                break
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
EOFPYTHON acc 2>/dev/null | while read url; do
        if [ -n "$url" ]; then
            download_file "$url" "$output" || rm -f "$output"
        fi
    done
done

echo ""
echo "Decompressing files..."
for f in *.gz; do
    [ -f "$f" ] && [ -s "$f" ] && gunzip -v "$f" || true
done

echo ""
echo "=== Final files ==="
ls -lh ENCODE_* 2>/dev/null || echo "✗ No ENCODE peak files downloaded. Download manually or check network connection."
