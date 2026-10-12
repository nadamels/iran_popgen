#!/bin/bash

#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH --mail-user=nme3923@nyu.edu
#SBATCH --mail-type=END,FAIL

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

#Loading File as myenviro contains PLINK; the KING executable ("king") sits in the directory you run/submit from.

# Removes individuals related up to the 2nd degree (parent-offspring, full sibs,
# 2nd-degree relatives) using KING, and lists who was removed.
# KING's --unrelated picks a maximal set of mutually unrelated individuals, so as
# few samples as possible are dropped. The input fileset is never modified.
#
# Usage: sbatch remove_related.sh <input_prefix> <output_prefix> [degree=2]
#   prefixes may be paths in any directory (no .bed/.bim/.fam)
#
# Outputs:
#   <output_prefix>.{bed,bim,fam}          unrelated samples only
#   <output_prefix>.removed_related.txt    FID IID of every removed sample
#   <output_prefix>.king_unrelated.txt     KING's own kept-sample list
#
# Note: KING groups people by FID. If every FID is 0 (as in your merged data) all
# samples are treated as one "family" and KING infers relatedness from the
# genotypes alone, which is what you want for degree-based pruning.

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
CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-10}}

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }

# KING is run as ./king, so submit/run this script from the directory containing it
[ -x ./king ] || { echo "ERROR: ./king not found or not executable in $(pwd); run from the directory containing king" >&2; exit 1; }

# ---- KING: list of unrelated individuals -----------------------------------
./king -b "${IN}.bed" --unrelated --degree "$DEGREE" --cpus "$CPUS" --prefix "${OUT}.king"

KEEP="${OUT}.king_unrelated.txt"
[ -s "$KEEP" ] || { echo "ERROR: KING did not produce $KEEP" >&2; exit 1; }

# ---- removed = everyone in the fam that KING did not keep (FID IID match) --
awk 'NR==FNR { keep[$1" "$2]=1; next } !(($1" "$2) in keep) { print $1, $2 }' "$KEEP" "${IN}.fam" > "${OUT}.removed_related.txt"

nrem=$(wc -l < "${OUT}.removed_related.txt")

# ---- write the pruned fileset ----------------------------------------------
plink --bfile "$IN" --keep "$KEEP" --allow-no-sex --make-bed --out "$OUT"

# ---- verify ----------------------------------------------------------------
nin=$(wc -l < "${IN}.fam"); nout=$(wc -l < "${OUT}.fam")
[ "$nout" -eq $(( nin - nrem )) ] || { echo "ERROR: expected $((nin - nrem)) samples, got $nout" >&2; exit 1; }

echo "done: removed ${nrem} related sample(s) (up to degree ${DEGREE}); ${nin} -> ${nout}"
echo "      removed list: ${OUT}.removed_related.txt"
