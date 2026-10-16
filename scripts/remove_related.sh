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
# After removing samples the script runs KING again on the rest and repeats until KING finds no
# related pair up to the degree (MAX_ROUNDS, default 5), so the output really is clear of relatives.
#
# Outputs:
#   <output_prefix>.{bed,bim,fam}                 unrelated samples only
#   <output_prefix>.removed_related.txt           FID IID of every removed sample (all rounds)
#   <output_prefix>.removed_related_details.tsv   why each was removed (round, kinship, relatives kept)
#   <output_prefix>.king.kin / .kin0              KING's pairwise tables for round 1
#   <output_prefix>.king_round2.kin ...           later rounds (only if round 1 left relatives)
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
[[ "$DEGREE" =~ ^[1-3]$ ]] || { echo "ERROR: degree must be 1, 2 or 3 (got '${3:-}')" >&2; exit 1; }
CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-10}}
# the helper python script lives next to this one (sbatch runs a spooled copy)
PRUNE_PY=/scratch/nme3923/IranGen/scripts/prune_related.py
[ -f "$PRUNE_PY" ] || PRUNE_PY="$(dirname "${BASH_SOURCE[0]}")/prune_related.py"
[ -f "$PRUNE_PY" ] || { echo "ERROR: prune_related.py not found" >&2; exit 1; }

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }

# ---- rounds of: KING -> minimum removal set -> remove -> re-check -------------------------------
# Kinship estimates shift a little once samples are removed (and KING's relationship call comes from
# shared IBD segments, not only the Kinship column), so after each removal KING is run again on what
# is left; this repeats until KING finds NO related pair up to the degree (or MAX_ROUNDS is hit).
MAX_ROUNDS=${MAX_ROUNDS:-5}
CUR="$IN"
: > "${OUT}.removed_related.txt"
printf 'round\tFID\tIID\tn_related\tmax_kinship\trelated_kept_samples\n' > "${OUT}.removed_related_details.tsv"
for ((r = 1; r <= MAX_ROUNDS; r++)); do
	KP="${OUT}.king"; [ "$r" -gt 1 ] && KP="${OUT}.king_round${r}"
	RP="${OUT}.round${r}"
	echo "=== round ${r}: KING on $(wc -l < "${CUR}.fam") samples ==="
	king -b "${CUR}.bed" --related --degree "$DEGREE" --cpus "$CPUS" --prefix "$KP"

	# missingness (used only to break ties)
	plink --bfile "$CUR" --missing --allow-no-sex --out "${OUT}.tmpmiss"
	# minimum set of samples to remove (pairs KING lists as related up to DEGREE, or kinship over the cutoff)
	python3 "$PRUNE_PY" "$DEGREE" "${OUT}.tmpmiss.imiss" "${KP}.kin" "${KP}.kin0" -- "$RP"
	rm -f "${OUT}.tmpmiss".*

	n=$(wc -l < "${RP}.removed_related.txt")
	cat "${RP}.removed_related.txt" >> "${OUT}.removed_related.txt"
	tail -n +2 "${RP}.removed_related_details.tsv" | awk -v r="$r" 'BEGIN{FS=OFS="\t"} {print r, $0}' >> "${OUT}.removed_related_details.tsv"
	rm -f "${RP}.removed_related.txt" "${RP}.removed_related_details.tsv"

	if [ "$n" -eq 0 ]; then
		echo "round ${r}: KING finds no related pair up to degree ${DEGREE}; done"
		break
	fi
	tail -n "$n" "${OUT}.removed_related.txt" > "${OUT}.tmp_remove${r}.txt"      # this round's removals
	plink --bfile "$CUR" --remove "${OUT}.tmp_remove${r}.txt" --allow-no-sex --make-bed --out "${OUT}.tmp_round${r}"
	CUR="${OUT}.tmp_round${r}"
	if [ "$r" -eq "$MAX_ROUNDS" ]; then
		echo "WARNING: still removing samples after ${MAX_ROUNDS} rounds; re-run KING on ${OUT} to check what is left" >&2
	fi
done
nrem=$(wc -l < "${OUT}.removed_related.txt")

# ---- write the pruned fileset ----------------------------------------------
plink --bfile "$CUR" --allow-no-sex --make-bed --out "$OUT"
rm -f "${OUT}".tmp_round* "${OUT}".tmp_remove*

# ---- verify ----------------------------------------------------------------
nin=$(wc -l < "${IN}.fam"); nout=$(wc -l < "${OUT}.fam")
[ "$nout" -eq $(( nin - nrem )) ] || { echo "ERROR: expected $((nin - nrem)) samples, got $nout" >&2; exit 1; }

echo "done: removed ${nrem} sample(s) to clear relatedness up to degree ${DEGREE}; ${nin} -> ${nout}"
echo "      removed list: ${OUT}.removed_related.txt"
echo "      details     : ${OUT}.removed_related_details.tsv"
