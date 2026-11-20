#!/bin/bash
#SBATCH --job-name=dl_phased_vcfs
#SBATCH -c 2
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=dl_phased_vcfs_%j.out
#SBATCH --error=dl_phased_vcfs_%j.err

# Download every *.vcf.gz (plus any .tbi/.csi index and md5/checksum files, if present) from the
# Sanger APPG phased-VCF directory.
#
# Usage:   sbatch download_phased_vcfs.sh [outdir]
#          bash   download_phased_vcfs.sh [outdir]
#
#   outdir   default: $SLURM_SUBMIT_DIR/PHASED_VCFs (or ./PHASED_VCFs)
#
# The file list is scraped from the directory index at run time, so nothing is hard-coded.
# Downloads are resumable (wget -c) and skip files that are already complete, so you can simply
# resubmit the job if it is interrupted or hits the time limit.
# Note: compute nodes need outbound internet access; if yours don't, run this on a login/transfer node
# with `bash download_phased_vcfs.sh`.

set -euo pipefail

BASE_URL="https://ngs.sanger.ac.uk/production/appg/PHASED_VCFs/"
OUTDIR="${1:-${SLURM_SUBMIT_DIR:-$PWD}/PHASED_VCFs}"
INCLUDE_INDEXES=1   # 1 = also fetch .vcf.gz.tbi / .csi and checksum files; 0 = only .vcf.gz

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "[$(date)] Fetching directory index: $BASE_URL"
wget -q --tries=5 --waitretry=10 --timeout=60 -O index.html "$BASE_URL"

# Pull every href out of the index page, drop query strings/anchors and parent/sort links.
if (( INCLUDE_INDEXES )); then
    PATTERN='\.vcf\.gz(\.tbi|\.csi|\.md5)?$|(md5|MD5|checksums?)[^/]*$'
else
    PATTERN='\.vcf\.gz$'
fi

grep -oE 'href="[^"?#]+"' index.html \
    | sed -E 's/^href="//; s/"$//' \
    | grep -vE '^(/|\.\./|\?)' \
    | grep -E "$PATTERN" \
    | sort -u > file_list.txt

N=$(wc -l < file_list.txt)
if [[ "$N" -eq 0 ]]; then
    echo "ERROR: no matching files found in $BASE_URL (see $OUTDIR/index.html)" >&2
    exit 1
fi
echo "[$(date)] Found $N files; downloading to $OUTDIR"

FAILED=0
while read -r f; do
    # Hrefs may be relative (usual) or absolute; handle both.
    case "$f" in
        http*) url="$f"; f=$(basename "$f") ;;
        *)     url="${BASE_URL}${f}" ;;
    esac
    # Percent-decoded name for the local file
    local_name=$(printf '%b' "${f//%/\\x}")

    echo "[$(date)] $local_name"
    # -c resume partial files; --timestamping is not combined with -c, so skip completed files explicitly
    # by comparing the remote Content-Length with the local size.
    remote_size=$(wget --spider --server-response --tries=3 --timeout=60 "$url" 2>&1 | awk 'tolower($1)=="content-length:"{s=$2} END{gsub(/\r/,"",s); print s}')
    if [[ -n "$remote_size" && -f "$local_name" && "$(stat -c %s "$local_name")" == "$remote_size" ]]; then
        echo "    already complete, skipping"
        continue
    fi
    if ! wget -c -q --show-progress --tries=10 --waitretry=30 --timeout=60 "$url"; then
        echo "    FAILED: $local_name" >&2
        FAILED=$((FAILED + 1))
    fi
done < file_list.txt

# Verify gzip integrity of everything we got.
echo "[$(date)] Verifying gzip integrity"
for v in *.vcf.gz; do
    [[ -e "$v" ]] || continue
    gzip -t "$v" 2>/dev/null || { echo "    CORRUPT/INCOMPLETE: $v" >&2; FAILED=$((FAILED + 1)); }
done

if [[ "$FAILED" -gt 0 ]]; then
    echo "[$(date)] Finished with $FAILED problem(s); resubmit to resume." >&2
    exit 1
fi
echo "[$(date)] All done."
