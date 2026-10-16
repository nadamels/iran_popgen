#!/bin/bash

#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=nme3923@nyu.edu
#SBATCH --mail-type=END,FAIL

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

#Loading File as myenviro contains PLINK (see the package list below).

# Classical multidimensional scaling (MDS / cmdscale) of a PLINK binary fileset.
#   1. (optional) LD-prune the variants (autosomes only)            plink --indep-pairwise
#   2. pairwise IBS distance matrix, D = 1 - IBS                     plink --distance square bin 1-ibs
#   3. classical MDS: double-centre -D^2/2, take the top eigenvectors (numpy / scipy)
# The input fileset is never modified.
#
# Usage: sbatch mds.sh <input_prefix> <output_prefix> [n_dims=10]
#   prefixes may be paths in any directory (no .bed/.bim/.fam)
#
# Outputs (all next to <output_prefix>):
#   <out>.mds.eigenvec   FID IID MDS1..MDSk   (PLINK 1.9 .eigenvec layout, no header)
#   <out>.mds.eigenval   the k eigenvalues, one per line
#   <out>.mds.fam        copy of the input .fam (same sample order)
#   -> so the plotting script works directly:  plot_pca.sh <out>.mds <mega_table> <pca_colors.txt>
#
# Settings (environment variables):
#   LD_PRUNE=1        1 = LD-prune first (default), 0 = use every variant as is
#   LD_WINDOW="200 25 0.4"   plink --indep-pairwise window / step / r^2
#   KEEP_TMP=0        1 = keep the intermediate pruned fileset and distance matrix
#
# Needs in the conda env: plink (1.9), python3, numpy; scipy is optional (faster eigensolver).
# Memory: the distance matrix is n x n float64 (about 8*n^2 bytes, e.g. 0.25 GB for 5,500 samples);
# the SBATCH --mem above is generous.
#
# Install (once):
#   conda install -n myenviro -c conda-forge -c bioconda plink numpy scipy pandas matplotlib openpyxl
#   (pandas / matplotlib / openpyxl are only for plot_pca.sh; the MDS itself needs plink + numpy [+ scipy])

set -euo pipefail

if [ "${#ARGS[@]}" -lt 2 ] || [ "${#ARGS[@]}" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <output_prefix> [n_dims=10]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

IN=$(strip_ext "${ARGS[0]}")
mkdir -p "$(dirname "$(strip_ext "${ARGS[1]}")")"
IN=$(abspath "$IN")
OUT=$(abspath "$(strip_ext "${ARGS[1]}")")
NDIMS=${ARGS[2]:-10}
THREADS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-4}}
LD_PRUNE=${LD_PRUNE:-1}
LD_WINDOW=${LD_WINDOW:-"200 25 0.4"}
KEEP_TMP=${KEEP_TMP:-0}

[[ "$NDIMS" =~ ^[0-9]+$ && "$NDIMS" -ge 2 ]] || { echo "ERROR: n_dims must be an integer >= 2 (got '$NDIMS')" >&2; exit 1; }
for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ "$IN" != "$OUT" ] || { echo "ERROR: output prefix must differ from input" >&2; exit 1; }
command -v plink >/dev/null || { echo "ERROR: plink not found in the active environment" >&2; exit 1; }
python3 -c "import numpy" 2>/dev/null || { echo "ERROR: python3 needs numpy (conda install numpy)" >&2; exit 1; }

dups=$(awk '{print $1" "$2}' "${IN}.fam" | sort | uniq -d | wc -l)
[ "$dups" -eq 0 ] || { echo "ERROR: ${dups} duplicated FID/IID pair(s) in ${IN}.fam; run remove_duplicates.sh first" >&2; exit 1; }
N=$(wc -l < "${IN}.fam")
echo "input  : $IN  ($N samples, $(wc -l < "${IN}.bim") variants)"
echo "output : ${OUT}.mds.*   dims: $NDIMS   threads: $THREADS   LD prune: $LD_PRUNE"

TMP=$(mktemp -d "${OUT}.mdstmp.XXXXXX")
if [ "$KEEP_TMP" = "1" ]; then echo "tmp kept: $TMP"; else trap 'rm -rf "$TMP"' EXIT; fi

# ---- 1. LD pruning (autosomes) ---------------------------------------------------------------
SRC="$IN"
if [ "$LD_PRUNE" = "1" ]; then
	# shellcheck disable=SC2086
	plink --bfile "$IN" --autosome --allow-no-sex --threads "$THREADS" \
		--indep-pairwise $LD_WINDOW --out "${TMP}/ld"
	plink --bfile "$IN" --extract "${TMP}/ld.prune.in" --make-bed --allow-no-sex --threads "$THREADS" --out "${TMP}/pruned"
	SRC="${TMP}/pruned"
	echo "variants after LD pruning: $(wc -l < "${SRC}.bim")"
fi

# ---- 2. IBS distance matrix (1 - IBS), square, binary float -----------------------------------
plink --bfile "$SRC" --autosome --allow-no-sex --threads "$THREADS" \
	--distance square bin 1-ibs --out "${TMP}/dist"
[ -s "${TMP}/dist.mdist.bin" ] || { echo "ERROR: plink did not write the distance matrix" >&2; exit 1; }

# ---- 3. classical MDS ------------------------------------------------------------------------
python3 - "${TMP}/dist.mdist.bin" "${TMP}/dist.mdist.id" "$NDIMS" "${OUT}.mds" <<'PY'
import sys
import numpy as np

dist_bin, id_file, k, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
ids = [l.split()[:2] for l in open(id_file) if l.strip()]
n = len(ids)
D = np.fromfile(dist_bin, dtype=np.float32)
if D.size != n * n:
    sys.exit(f"ERROR: distance matrix has {D.size} values, expected {n}x{n}={n*n}")
D = D.reshape(n, n).astype(np.float64)
if not np.isfinite(D).all():
    sys.exit("ERROR: non-finite distances (a pair of samples with no overlapping genotypes?)")
D = (D + D.T) / 2                                    # remove float asymmetry
k = min(k, n - 1)

# double centring of -D^2/2  (B = -1/2 J D^2 J)
D2 = D * D
rm = D2.mean(axis=1, keepdims=True)
B = -0.5 * (D2 - rm - rm.T + D2.mean())
trace = float(np.trace(B))

try:
    from scipy.sparse.linalg import eigsh
    vals, vecs = eigsh(B, k=k, which="LA")
except Exception:                                     # scipy missing or no convergence
    vals, vecs = np.linalg.eigh(B)
    vals, vecs = vals[-k:], vecs[:, -k:]
order = np.argsort(vals)[::-1]
vals, vecs = vals[order], vecs[:, order]
pos = vals > 0
if not pos.all():
    print(f"NOTE: only {int(pos.sum())} positive eigenvalues among the top {k}; keeping those")
vals, vecs = vals[pos], vecs[:, pos]
coords = vecs * np.sqrt(vals)
# make the sign deterministic: the largest-magnitude coordinate of each axis is positive
sign = np.sign(coords[np.abs(coords).argmax(axis=0), np.arange(coords.shape[1])])
coords *= np.where(sign == 0, 1, sign)

with open(out + ".eigenvec", "w") as fh:
    for (fid, iid), row in zip(ids, coords):
        fh.write(f"{fid} {iid} " + " ".join(f"{x:.6g}" for x in row) + "\n")
np.savetxt(out + ".eigenval", vals, fmt="%.8g")
print(f"MDS done: {n} samples, {coords.shape[1]} dimensions")
for i, v in enumerate(vals[:10], 1):
    print(f"  MDS{i}: eigenvalue {v:.5g}  ({100 * v / trace:.2f}% of total variance)")
PY

# fam in the same sample order as the eigenvec (needed by plot_pca.sh)
cp "${IN}.fam" "${OUT}.mds.fam"
cmp -s <(awk '{print $1" "$2}' "${IN}.fam") <(awk '{print $1" "$2}' "${OUT}.mds.eigenvec") \
	|| echo "NOTE: eigenvec order differs from the fam; plot_pca.sh joins on IID so this is fine"

echo "done: ${OUT}.mds.eigenvec  ${OUT}.mds.eigenval  ${OUT}.mds.fam"
echo "plot: bash plot_pca.sh ${OUT}.mds <mega_table> <pca_colors.txt>"
