#!/bin/bash

#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH --mail-user=nme3923@nyu.edu
#SBATCH --mail-type=END,FAIL

# `source`d scripts inherit our arguments and conda's activate rejects them (ArgumentError): hide, then restore
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro
set -- "${ARGS[@]}"

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
# If MAX_ROUNDS runs out right after a removal, one more KING run checks the result; if related pairs
# remain the script exits with status 2 and lists them (never reports success on an unverified set).
# Every KING run is also cross-checked against KING's own "Relationship summary": if KING reports related
# pairs that are missing from its .kin table, the script stops (status 3) instead of calling the set clean.
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

# KING prints "Inference  MZ  PO  FS  2nd  3rd  OTHER"; return how many related pairs it reports up to $DEGREE
# (empty if the summary line is not found, e.g. a different KING version)
king_expected_pairs() {
	awk -v d="$DEGREE" '/^[[:space:]]*Inference[[:space:]]/ { n=$2+$3+$4; if (d>=2) n+=$5; if (d>=3) n+=$6; print n; exit }' "$1"
}

# king_and_prune <fileset> <king prefix> <prune prefix>
#   1. KING --related; 2. minimum removal set from its tables; 3. cross-check against KING's own summary.
#   KING's .kin table can come out INCOMPLETE (seen: cut off at exactly 2048 rows, i.e. only the first
#   2048 of all pairs), so pairs KING counts can be missing from it. If that happens, KING --ibdseg (which
#   writes only the related pairs) is run and its .seg table is added, then the removal set is recomputed.
#   Sets: N_REMOVE (samples to remove), EXPECTED (KING's pair count), UNPARSED (pairs still missing from the tables).
king_and_prune() {
	local fs="$1" kp="$2" rp="$3" npairs
	king -b "${fs}.bed" --related --degree "$DEGREE" --cpus "$CPUS" --prefix "$kp" 2>&1 | tee "${kp}.log"
	plink --bfile "$fs" --missing --allow-no-sex --out "${OUT}.tmpmiss"
	python3 "$PRUNE_PY" "$DEGREE" "${OUT}.tmpmiss.imiss" "${kp}.kin" "${kp}.kin0" -- "$rp"
	npairs=$(cat "${rp}.npairs" 2>/dev/null || echo 0)
	EXPECTED=$(king_expected_pairs "${kp}.log")
	UNPARSED=0
	if [ -n "$EXPECTED" ] && [ "$EXPECTED" -gt "$npairs" ]; then
		echo "WARNING: KING's own summary reports ${EXPECTED} related pair(s) up to degree ${DEGREE}, but only ${npairs} were read from ${kp}.kin / .kin0" >&2
		echo "         (the table is incomplete: $(($(wc -l < "${kp}.kin" 2>/dev/null || echo 1) - 1)) rows for $(wc -l < "${fs}.fam") samples). Trying KING --ibdseg to recover the missing pairs." >&2
		king -b "${fs}.bed" --ibdseg --degree "$DEGREE" --cpus "$CPUS" --prefix "${kp}_seg" 2>&1 | tee "${kp}_seg.log"
		python3 "$PRUNE_PY" "$DEGREE" "${OUT}.tmpmiss.imiss" "${kp}.kin" "${kp}.kin0" "${kp}_seg.seg" -- "$rp"
		npairs=$(cat "${rp}.npairs" 2>/dev/null || echo 0)
		if [ "$EXPECTED" -gt "$npairs" ]; then
			UNPARSED=$((EXPECTED - npairs))
			echo "WARNING: still ${UNPARSED} pair(s) short after --ibdseg ($npairs of ${EXPECTED} found)." >&2
		else
			echo "recovered: ${npairs} related pair(s) found (KING reports ${EXPECTED})"
		fi
	fi
	rm -f "${rp}.npairs" "${OUT}.tmpmiss".*
	N_REMOVE=$(wc -l < "${rp}.removed_related.txt")
}

# ---- rounds of: KING -> minimum removal set -> remove -> re-check -------------------------------
# Kinship estimates shift a little once samples are removed (and KING's relationship call comes from
# shared IBD segments, not only the Kinship column), so after each removal KING is run again on what
# is left; this repeats until KING finds NO related pair up to the degree (or MAX_ROUNDS is hit).
MAX_ROUNDS=${MAX_ROUNDS:-5}
CUR="$IN"
CLEAN=0                     # 1 once a KING run on the CURRENT fileset has found no related pair
: > "${OUT}.removed_related.txt"
printf 'round\tFID\tIID\tn_related\tmax_kinship\trelated_kept_samples\n' > "${OUT}.removed_related_details.tsv"
for ((r = 1; r <= MAX_ROUNDS; r++)); do
	KP="${OUT}.king"; [ "$r" -gt 1 ] && KP="${OUT}.king_round${r}"
	RP="${OUT}.round${r}"
	echo "=== round ${r}: KING on $(wc -l < "${CUR}.fam") samples ==="
	king_and_prune "$CUR" "$KP" "$RP"
	n=$N_REMOVE
	expected=$EXPECTED
	cat "${RP}.removed_related.txt" >> "${OUT}.removed_related.txt"
	tail -n +2 "${RP}.removed_related_details.tsv" | awk -v r="$r" 'BEGIN{FS=OFS="\t"} {print r, $0}' >> "${OUT}.removed_related_details.tsv"
	rm -f "${RP}.removed_related.txt" "${RP}.removed_related_details.tsv"

	if [ "$n" -eq 0 ]; then
		if [ "$UNPARSED" -gt 0 ]; then
			echo "ERROR: round ${r}: nothing to remove according to the table, but KING itself reports ${expected} related pair(s) in ${CUR}." >&2
			echo "       Not claiming the result is clean. Look at ${KP}.log and ${KP}.kin (are those pairs in the table? what do their InfType / ID columns look like?)." >&2
			exit 3
		fi
		echo "round ${r}: KING finds no related pair up to degree ${DEGREE}; done"
		CLEAN=1
		break
	fi
	tail -n "$n" "${OUT}.removed_related.txt" > "${OUT}.tmp_remove${r}.txt"      # this round's removals
	plink --bfile "$CUR" --remove "${OUT}.tmp_remove${r}.txt" --allow-no-sex --make-bed --out "${OUT}.tmp_round${r}"
	CUR="${OUT}.tmp_round${r}"
done
nrem=$(wc -l < "${OUT}.removed_related.txt")

# ---- write the pruned fileset ----------------------------------------------
plink --bfile "$CUR" --allow-no-sex --make-bed --out "$OUT"
rm -f "${OUT}".tmp_round* "${OUT}".tmp_remove*

# ---- final verification ------------------------------------------------------
# If the loop ended because MAX_ROUNDS ran out, the last removal has NOT been checked by KING yet: check now.
# A run that ended "clean" was already verified (its last KING run was on exactly this fileset).
if [ "$CLEAN" -ne 1 ]; then
	echo "=== final check: MAX_ROUNDS (${MAX_ROUNDS}) reached after a removal; running KING on the result ==="
	king_and_prune "$OUT" "${OUT}.king_final" "${OUT}.final_check"
	nleft=$N_REMOVE
	if [ "$nleft" -eq 0 ] && [ "$UNPARSED" -gt 0 ]; then
		echo "ERROR: final check: KING reports ${EXPECTED} related pair(s) in ${OUT} that could not be read from its tables; result NOT verified clean (see ${OUT}.king_final.log)." >&2
		exit 3
	fi
	if [ "$nleft" -gt 0 ]; then
		echo "ERROR: ${OUT} STILL contains related pairs up to degree ${DEGREE} (KING final check)." >&2
		echo "       ${nleft} more sample(s) would have to go; they are listed in ${OUT}.final_check.removed_related.txt" >&2
		echo "       KING's pairs: ${OUT}.king_final.kin. Re-run with a higher MAX_ROUNDS, e.g.  MAX_ROUNDS=10 sbatch ..." >&2
		exit 2
	fi
	rm -f "${OUT}.final_check.removed_related.txt" "${OUT}.final_check.removed_related_details.tsv"
	echo "final check: KING finds no related pair up to degree ${DEGREE}"
fi

# ---- verify ----------------------------------------------------------------
nin=$(wc -l < "${IN}.fam"); nout=$(wc -l < "${OUT}.fam")
[ "$nout" -eq $(( nin - nrem )) ] || { echo "ERROR: expected $((nin - nrem)) samples, got $nout" >&2; exit 1; }

echo "done: removed ${nrem} sample(s) to clear relatedness up to degree ${DEGREE}; ${nin} -> ${nout}"
echo "      removed list: ${OUT}.removed_related.txt"
echo "      details     : ${OUT}.removed_related_details.tsv"
