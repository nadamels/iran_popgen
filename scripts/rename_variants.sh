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

#Loading File as myenviro contains PLINK as well as most dependanceis for popgen.

# Renames the variant ID (.bim column 2) to <chromosome>_<position>
# (column 1 _ column 4), e.g.  "1  .  0  10542  T  C"  ->  "1  1_10542  0  10542  T  C"
#
# Variant order is unchanged, so the .bed is byte-identical and is simply
# copied alongside the new .bim and the .fam. The input is never modified.
#
# Usage: sbatch rename_variants.sh <input_prefix> [output_prefix]
#   prefix = path to the fileset without .bed/.bim/.fam (any directory)
#   output_prefix defaults to <input_prefix>.renamed

set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
	echo "Usage: $0 <input_prefix> [output_prefix]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }

IN=$(strip_ext "$1")
OUT=$(strip_ext "${2:-${IN}.renamed}")

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"

# bim is tab-delimited; keep tabs. Fail on malformed rows.
awk 'BEGIN{FS=OFS="\t"} NF!=6 {print "ERROR: line " NR " has " NF " columns" > "/dev/stderr"; exit 1} {$2=$1"_"$4; print}' \
	"${IN}.bim" > "${OUT}.bim"

nin=$(wc -l < "${IN}.bim")
nout=$(wc -l < "${OUT}.bim")
[ "$nin" -eq "$nout" ] || { echo "ERROR: row count mismatch ($nin vs $nout)" >&2; rm -f "${OUT}.bim"; exit 1; }

cp "${IN}.bed" "${OUT}.bed"
cp "${IN}.fam" "${OUT}.fam"

# Duplicate chr_pos IDs happen at multi-allelic sites (same position, different alleles)
dups=$(cut -f 2 "${OUT}.bim" | sort | uniq -d | wc -l)
if [ "$dups" -gt 0 ]; then
	echo "WARNING: ${dups} chr_pos IDs occur more than once (multi-allelic positions)." >&2
	echo "         PLINK merging/extract may complain about duplicate IDs." >&2
fi

echo "done: wrote ${OUT}.{bed,bim,fam} (${nout} variants)"
