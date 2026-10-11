#!/bin/bash

#SBATCH -c 10
#SBATCH -t 02:00:00

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

#Loading File as myenviro contains PLINK.

# Applies a MAF filter to an existing (already filtered / merged) fileset.
# MAF is recomputed on the samples in the input, and the input is not modified.
#
# Usage: sbatch maf_filter.sh <input_prefix> <output_prefix> [maf=0.05]
#   prefixes may be paths in any directory (no .bed/.bim/.fam)

set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <output_prefix> [maf=0.05]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
IN=$(strip_ext "$1")
OUT=$(strip_ext "$2")
MAF=${3:-0.05}

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

plink --bfile "$IN" --maf "$MAF" --allow-no-sex --make-bed --out "$OUT"

echo "variants before: $(wc -l < "${IN}.bim")  after: $(wc -l < "${OUT}.bim")"
