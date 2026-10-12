#!/bin/bash

#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH --mail-user=nme3923@nyu.edu
#SBATCH --mail-type=END,FAIL

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

#Loading File as myenviro contains PLINK and python; the KING executable ("king") is in the scripts directory.

# Removes the MINIMUM number of samples needed so that no pair is related up to
# the 2nd degree (parent-offspring, full sibs, 2nd-degree relatives), and lists
# who was removed. The input fileset is never modified.
#
# KING (--related) gives the pairwise kinship table; prune_related.py then finds
# the smallest set of samples whose removal breaks every related pair (e.g. in a
# trio/cluster only the "hub" sample is removed, not everyone). Ties are broken
# by removing the sample with more missing data.
#
# Usage: sbatch remove_related.sh <input_prefix> <output_prefix> [degree=2]
#   prefixes may be paths in any directory (no .bed/.bim/.fam)
#
# Outputs:
#   <output_prefix>.{bed,bim,fam}                 unrelated samples only
#   <output_prefix>.removed_related.txt           FID IID of every removed sample
#   <output_prefix>.removed_related_details.tsv   why each was removed (kinship, relatives kept)
#   <output_prefix>.king.kin / .kin0              KING's pairwise tables
#
# Note: KING groups people by FID. If every FID is 0 all samples are treated as
# one "family" and relatedness is inferred from the genotypes alone.
export PATH="/scratch/nme3923/IranGen/scripts:$PATH"

set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <output_prefix> [degree=2]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

IN=$(strip_ext "$1")
mkdir -p "$(dirname "$(strip_ext "$2")")"
IN=$(abspath "$IN")
OUT=$(abspath "$(strip_ext "$2")")
DEGREE=${3:-2}
DEGREE=${DEGREE#degree=}   # accept "2" or "degree=2"
[[ "$DEGREE" =~ ^[1-3]$ ]] || { echo "ERROR: degree must be 1, 2 or 3 (got '$3')" >&2; exit 1; }
CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-10}}
# the helper python script lives next to this one (sbatch runs a spooled copy)
PRUNE_PY=/scratch/nme3923/IranGen/scripts/prune_related.py
[ -f "$PRUNE_PY" ] || PRUNE_PY="$(dirname "${BASH_SOURCE[0]}")/prune_related.py"
[ -f "$PRUNE_PY" ] || { echo "ERROR: prune_related.py not found" >&2; exit 1; }

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }

# ---- KING: pairwise kinship up to the chosen degree ------------------------
king -b "${IN}.bed" --related --degree "$DEGREE" --cpus "$CPUS" --prefix "${OUT}.king"

# ---- missingness (used only to break ties) ---------------------------------
plink --bfile "$IN" --missing --allow-no-sex --out "${OUT}.tmpmiss"

# ---- minimum set of samples to remove --------------------------------------
python3 "$PRUNE_PY" "$DEGREE" "${OUT}.tmpmiss.imiss" "${OUT}.king.kin" "${OUT}.king.kin0" -- "$OUT"
rm -f "${OUT}.tmpmiss".*

nrem=$(wc -l < "${OUT}.removed_related.txt")

# ---- write the pruned fileset ----------------------------------------------
if [ "$nrem" -gt 0 ]; then
	plink --bfile "$IN" --remove "${OUT}.removed_related.txt" --allow-no-sex --make-bed --out "$OUT"
else
	plink --bfile "$IN" --allow-no-sex --make-bed --out "$OUT"
fi

# ---- verify ----------------------------------------------------------------
nin=$(wc -l < "${IN}.fam"); nout=$(wc -l < "${OUT}.fam")
[ "$nout" -eq $(( nin - nrem )) ] || { echo "ERROR: expected $((nin - nrem)) samples, got $nout" >&2; exit 1; }

echo "done: removed ${nrem} sample(s) to clear relatedness up to degree ${DEGREE}; ${nin} -> ${nout}"
echo "      removed list: ${OUT}.removed_related.txt"
echo "      details     : ${OUT}.removed_related_details.tsv"
