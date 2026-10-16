#!/bin/bash

#SBATCH -c 4
#SBATCH -t 02:00:00

# `source`d scripts inherit our arguments and conda's activate rejects them (ArgumentError): hide, then restore
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro
set -- "${ARGS[@]}"

#Loading File as myenviro contains PLINK.

# Trims everything from the first ".chr" onward in the IIDs (column 2 of the .fam),
# e.g.  I0698.chr_trimmed.sorted.nochr  ->  I0698
# and writes a new .bed/.bim/.fam (via PLINK --update-ids) so the .bed set carries the
# new IDs. FIDs and all other columns are unchanged; the input is never modified.
#
# Usage: sbatch trim_iids.sh <input_prefix> <output_prefix> [trim_from=.chr]
#   prefixes may be paths in any directory (no .bed/.bim/.fam)
#   trim_from: IIDs are cut at the first occurrence of this text (default ".chr")
#
# Also writes <output_prefix>.iid_changes.txt (old IID -> new IID).

set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <output_prefix> [trim_from=.chr]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

IN=$(strip_ext "$1")
mkdir -p "$(dirname "$(strip_ext "$2")")"
IN=$(abspath "$IN")
OUT=$(abspath "$(strip_ext "$2")")
CUT=${3:-.chr}

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }

TMP="${OUT}.update_ids.tmp"
trap 'rm -f "$TMP"' EXIT

# PLINK --update-ids format: old_FID old_IID new_FID new_IID (FID kept as is)
awk -v cut="$CUT" '{ new=$2; i=index($2,cut); if (i>1) new=substr($2,1,i-1); print $1, $2, $1, new }' "${IN}.fam" > "$TMP"

# Safety: trimmed IIDs must stay unique and non-empty
dups=$(awk '{print $3" "$4}' "$TMP" | sort | uniq -d | wc -l)
if [ "$dups" -gt 0 ]; then
	echo "ERROR: trimming would create ${dups} duplicated FID/IID pair(s); nothing written." >&2
	awk '{print $3" "$4}' "$TMP" | sort | uniq -d | head >&2
	exit 1
fi

awk '$2!=$4 {print $2, $4}' "$TMP" > "${OUT}.iid_changes.txt"
nchg=$(wc -l < "${OUT}.iid_changes.txt")

plink --bfile "$IN" --update-ids "$TMP" --allow-no-sex --make-bed --out "$OUT"

# Verify: same rows in same order, IIDs equal the trimmed names
[ "$(wc -l < "${OUT}.fam")" -eq "$(wc -l < "${IN}.fam")" ] || { echo "ERROR: fam row count changed" >&2; exit 1; }
paste -d' ' "$TMP" <(awk '{print $1, $2}' "${OUT}.fam") |
	awk '$3!=$5 || $4!=$6 {bad++} END{ if (bad) {print "ERROR: " bad " rows not renamed as expected" > "/dev/stderr"; exit 1} }'

echo "done: trimmed ${nchg} IID(s) at '${CUT}'; list: ${OUT}.iid_changes.txt"
