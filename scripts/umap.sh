#!/bin/bash

#SBATCH -c 10
#SBATCH -t 12:00:00
#SBATCH --mem=64G

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/myenviro

# UMAP of a PLINK binary fileset, restricted to the samples that are in BOTH the fileset and the
# mega sample table. The input fileset can therefore be any subset of the samples: only the IIDs that
# are in <prefix>.fam AND in the table's IID column are used (the rest are reported and ignored).
#   1. intersect fam IIDs with the table IIDs -> keep list
#   2. subset + (optional) LD-prune the variants (autosomes only)       plink --keep, --indep-pairwise
#   3. PCA of the subset                                               plink --pca
#   4. UMAP on the top PCs                                             umap-learn
#   5. plot, coloured by Super-population (colours from Plotting/pca_colors.txt if available)
# The input fileset is never modified.
#
# Usage:   sbatch umap.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]
#          bash   umap.sh <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]
#   <input_prefix>   PLINK fileset without .bed/.bim/.fam
#   output_prefix    default: <input_prefix>.umap
#
# Outputs:  <out>.keep.txt        FID IID of the samples used (fam ∩ table)
#           <out>.coords.tsv      IID, UMAP1, UMAP2, metadata columns from the table
#           <out>.png / .pdf      UMAP plot
#
# Settings (environment variables):
#   NPCS=20           number of PCs fed into UMAP
#   N_NEIGHBORS=15    UMAP n_neighbors
#   MIN_DIST=0.1      UMAP min_dist
#   SEED=42           random seed
#   KEEP_PCA=0        1 = also keep the PLINK PCA (<out>.pca.eigenvec / .eigenval); by default it is
#                     written to a temporary directory and deleted, so it cannot clash with other eigenvecs
#   LD_PRUNE=1        1 = LD-prune first (default), 0 = use every variant
#   LD_WINDOW="200 25 0.4"   plink --indep-pairwise window / step / r^2
#   COLORS_FILE       "name: #hex" lines (default Plotting/pca_colors.txt under $SLURM_SUBMIT_DIR or cwd)
#   COLOR_BY="Super-population"   table column used for the colours
#
# Needs in the conda env: plink (1.9), python3, numpy, pandas, matplotlib, umap-learn
#   (openpyxl too if the table is .xlsx):
#   conda install -n myenviro -c conda-forge -c bioconda plink numpy pandas matplotlib openpyxl umap-learn

set -euo pipefail

if [ "${#ARGS[@]}" -lt 2 ] || [ "${#ARGS[@]}" -gt 3 ]; then
	echo "Usage: $0 <input_prefix> <mega_sample_table.xlsx|.csv> [output_prefix]" >&2
	exit 1
fi

strip_ext() { local f=${1%.bed}; f=${f%.bim}; f=${f%.fam}; echo "$f"; }
abspath() { local d; d=$(cd "$(dirname "$1")" && pwd) && echo "${d}/$(basename "$1")"; }

IN=$(abspath "$(strip_ext "${ARGS[0]}")")
TABLE="${ARGS[1]}"
OUT="${ARGS[2]:-${IN}.umap}"
mkdir -p "$(dirname "$OUT")"
OUT=$(abspath "$OUT")
BASE="${SLURM_SUBMIT_DIR:-$PWD}"
COLORS="${COLORS_FILE:-$BASE/Plotting/pca_colors.txt}"
COLOR_BY="${COLOR_BY:-Super-population}"
NPCS=${NPCS:-20}
N_NEIGHBORS=${N_NEIGHBORS:-15}
MIN_DIST=${MIN_DIST:-0.1}
SEED=${SEED:-42}
KEEP_PCA=${KEEP_PCA:-0}
THREADS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-4}}
LD_PRUNE=${LD_PRUNE:-1}
LD_WINDOW=${LD_WINDOW:-"200 25 0.4"}

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
	alt="${TABLE%.*}.csv"
	if [[ -s "$alt" ]]; then
		echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
		TABLE="$alt"
	fi
fi

for ext in bed bim fam; do
	[ -s "${IN}.${ext}" ] || { echo "ERROR: ${IN}.${ext} not found" >&2; exit 1; }
done
[ -s "$TABLE" ] || { echo "ERROR: table $TABLE not found" >&2; exit 1; }
[[ "$NPCS" =~ ^[0-9]+$ && "$NPCS" -ge 2 ]] || { echo "ERROR: NPCS must be an integer >= 2 (got '$NPCS')" >&2; exit 1; }
command -v plink >/dev/null || { echo "ERROR: plink not found in the active environment" >&2; exit 1; }
python3 -c "import numpy, pandas, matplotlib, umap" 2>/dev/null || {
	echo "ERROR: python3 needs numpy, pandas, matplotlib and umap-learn (conda install -c conda-forge umap-learn)" >&2; exit 1; }

echo "input  : $IN  ($(wc -l < "${IN}.fam") samples, $(wc -l < "${IN}.bim") variants)"
echo "table  : $TABLE"
echo "output : ${OUT}.*   PCs: $NPCS   n_neighbors: $N_NEIGHBORS   min_dist: $MIN_DIST   LD prune: $LD_PRUNE"

TMP=$(mktemp -d "${OUT}.tmp.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# ---- 1. fam IIDs  ∩  table IIDs -> keep list ---------------------------------------------------
python3 -I - "${IN}.fam" "$TABLE" "${OUT}.keep.txt" <<'PY'
import sys
import pandas as pd

fam_path, table_path, keep_path = sys.argv[1:4]
fam = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str, usecols=[0, 1], names=["FID", "IID"])
# Un-repaired fam where every IID is one constant string (ID stored in FID): fall back to FID.
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; matching the table on FID")
fam["sid"] = (fam["FID"] if use_fid else fam["IID"]).str.strip()

if table_path.lower().endswith((".xlsx", ".xls")):
    tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
else:
    tab = pd.read_csv(table_path, dtype=str)
if "IID" not in tab.columns:
    sys.exit(f"ERROR: {table_path} has no 'IID' column")
table_ids = set(tab["IID"].dropna().str.strip())

in_both = fam["sid"].isin(table_ids)
print(f"fam samples: {len(fam)}   table samples: {len(table_ids)}   intersection: {int(in_both.sum())}")
if (~in_both).any():
    miss = fam.loc[~in_both, "sid"]
    print(f"WARNING: {len(miss)} fam sample(s) not in the table are dropped, e.g. {', '.join(miss.head(10))}")
only_table = len(table_ids - set(fam["sid"]))
print(f"{only_table} table sample(s) are not in this fileset (ignored)")
if in_both.sum() < 10:
    sys.exit("ERROR: fewer than 10 samples in common between the fam and the table")
if fam.loc[in_both, "sid"].duplicated().any():
    sys.exit("ERROR: duplicated IIDs among the shared samples; run remove_duplicates.sh first")
fam.loc[in_both, ["FID", "IID"]].to_csv(keep_path, sep="\t", header=False, index=False)
PY

# ---- 2. subset + LD pruning (autosomes) --------------------------------------------------------
SRC="${TMP}/subset"
plink --bfile "$IN" --keep "${OUT}.keep.txt" --autosome --allow-no-sex --threads "$THREADS" \
	--make-bed --out "$SRC"
if [ "$LD_PRUNE" = "1" ]; then
	# shellcheck disable=SC2086
	plink --bfile "$SRC" --allow-no-sex --threads "$THREADS" --indep-pairwise $LD_WINDOW --out "${TMP}/ld"
	plink --bfile "$SRC" --extract "${TMP}/ld.prune.in" --make-bed --allow-no-sex --threads "$THREADS" --out "${TMP}/pruned"
	SRC="${TMP}/pruned"
fi
echo "variants used: $(wc -l < "${SRC}.bim")   samples: $(wc -l < "${SRC}.fam")"

# ---- 3. PCA of the subset ----------------------------------------------------------------------
plink --bfile "$SRC" --allow-no-sex --threads "$THREADS" --pca "$NPCS" --out "${TMP}/pca"

if [ "$KEEP_PCA" = "1" ]; then
	cp "${TMP}/pca.eigenvec" "${OUT}.pca.eigenvec"
	cp "${TMP}/pca.eigenval" "${OUT}.pca.eigenval"
fi

# ---- 4./5. UMAP + plot -------------------------------------------------------------------------
python3 -I - "${TMP}/pca.eigenvec" "$TABLE" "${IN}.fam" "$COLORS" "$COLOR_BY" "$NPCS" "$N_NEIGHBORS" "$MIN_DIST" "$SEED" "$OUT" <<'PY'
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import umap

eig_path, table_path, fam_path, colors_path, color_by = sys.argv[1:6]
npcs, n_neighbors, min_dist, seed = int(sys.argv[6]), int(sys.argv[7]), float(sys.argv[8]), int(sys.argv[9])
out = sys.argv[10]

ev = pd.read_csv(eig_path, sep=r"\s+", header=None, dtype=str)
pcs = [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
ev.columns = ["FID", "IID"] + pcs
ev[pcs] = ev[pcs].astype(float)
X = ev[pcs].to_numpy()

# Same matching rule as step 1 (IID, or FID when the fam's IIDs are not unique).
fam = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str, usecols=[0, 1], names=["FID", "IID"])
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
ev["sid"] = (ev["FID"] if use_fid else ev["IID"]).str.strip()

reducer = umap.UMAP(n_neighbors=min(n_neighbors, len(ev) - 1), min_dist=min_dist, n_components=2,
                    metric="euclidean", random_state=seed)
emb = reducer.fit_transform(X)
ev["UMAP1"], ev["UMAP2"] = emb[:, 0], emb[:, 1]

if table_path.lower().endswith((".xlsx", ".xls")):
    tab = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
else:
    tab = pd.read_csv(table_path, dtype=str)
tab["sid"] = tab["IID"].str.strip()
df = ev[["sid", "UMAP1", "UMAP2"]].merge(tab.drop(columns=["IID"]), on="sid", how="left")
df = df.rename(columns={"sid": "IID"})
df.to_csv(f"{out}.coords.tsv", sep="\t", index=False)

# ---- plot --------------------------------------------------------------------------------------
col_of = {}
if os.path.isfile(colors_path):
    for line in open(colors_path):
        if ":" in line:
            k, v = line.rsplit(":", 1)
            col_of[k.strip()] = v.strip()
else:
    print(f"NOTE: colours file {colors_path} not found; using default matplotlib colours")
if color_by not in df.columns:
    print(f"WARNING: column {color_by!r} not in the table; plotting without colours")
    df[color_by] = "all"
df[color_by] = df[color_by].fillna("unlabelled")
groups = sorted(df[color_by].unique(), key=lambda g: -(df[color_by] == g).sum())
cycle = plt.get_cmap("tab20")
fig, ax = plt.subplots(figsize=(10, 8))
for i, g in enumerate(groups):
    sub = df[df[color_by] == g]
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=6, alpha=0.8, linewidths=0,
               color=col_of.get(g, cycle(i % 20)), label=f"{g} (n={len(sub)})")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_title(f"{os.path.basename(out)}  |  n={len(df)}, {npcs} PCs, n_neighbors={n_neighbors}, min_dist={min_dist}")
ax.legend(markerscale=3, fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
fig.tight_layout()
fig.savefig(f"{out}.png", dpi=200)
fig.savefig(f"{out}.pdf")
print(f"wrote {out}.coords.tsv, {out}.png, {out}.pdf")
PY
