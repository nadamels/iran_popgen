#!/bin/bash

#SBATCH -c 2

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# Hierarchical clustering of samples on their top principal components, drawn as a dendrogram with
# "<sample ID> | <Population>" under every leaf (a colour strip and the text are coloured by Population).
#
# Usage:   sbatch cluster_pcs.sh <prefix> [mega_sample_table.xlsx|.csv] [n_pcs=10]
#          bash   cluster_pcs.sh <prefix> [mega_sample_table.xlsx|.csv] [n_pcs=10]
#
#   <prefix>   expects  <prefix>.fam  <prefix>.eigenvec  <prefix>.eigenval
#              (PLINK 1.9 or PLINK 2 eigenvec; with or without a header line)
#   table      default: $MASTER_TABLE, else metadata/mega_sample_table.xlsx under $SLURM_SUBMIT_DIR
#              (or the current directory); the .csv next to it is used if openpyxl is missing
#   n_pcs      how many leading PCs to cluster on (default 10)
#
# Which samples (environment variables):
#   SUPERPOP=Iran      Super-population to cluster (default Iran); SUPERPOP=all = every sample
#   ANCIENT=0          1 = include ancient samples too (default: modern only)
#   IRAN_SOURCE="Iranian cohort"   with SUPERPOP=Iran keep only the study's own Iranians, i.e. leave
#                      out Iran-labelled reference samples (e.g. SGDP); IRAN_SOURCE="" keeps them
# Clustering:
#   LINKAGE=ward       ward | average | complete | single | centroid | median | weighted
#                      (ward and the other tree methods use Euclidean distance on the PC coordinates)
#   N_CLUSTERS=        e.g. 5: cut the tree into 5 clusters (coloured branches + dashed cut line)
#   SCALE_PCS=0        1 = standardise every PC first (default 0: keep the PCA scale so PC1 counts most)
#
# Outputs (next to the prefix):  <prefix>.hclust_<superpop>.{png,pdf}   the dendrogram
#                                <prefix>.hclust_<superpop>.tsv         leaf order, Population and cluster
# A dendrogram with hundreds of leaves is unreadable: cluster a subset (SUPERPOP / ANCIENT).

set -euo pipefail

PREFIX="${ARGS[0]:-}"
if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <prefix> [mega_sample_table.xlsx|.csv] [n_pcs=10]" >&2
    exit 1
fi
BASE="${SLURM_SUBMIT_DIR:-$PWD}"
TABLE="${ARGS[1]:-${MASTER_TABLE:-$BASE/metadata/mega_sample_table.xlsx}}"
NPCS="${ARGS[2]:-10}"
[[ "$NPCS" =~ ^[0-9]+$ && "$NPCS" -ge 2 ]] || { echo "ERROR: n_pcs must be an integer >= 2 (got '$NPCS')" >&2; exit 1; }

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
    alt="${TABLE%.*}.csv"
    if [[ -s "$alt" ]]; then
        echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
        TABLE="$alt"
    fi
fi

for f in "${PREFIX}.fam" "${PREFIX}.eigenvec" "${PREFIX}.eigenval" "$TABLE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done
python3 -c "import numpy, pandas, scipy, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas, scipy and matplotlib" >&2; exit 1; }

SUPERPOP="${SUPERPOP:-Iran}"; ANCIENT="${ANCIENT:-0}"; LINKAGE="${LINKAGE:-ward}"
N_CLUSTERS="${N_CLUSTERS:-}"; SCALE_PCS="${SCALE_PCS:-0}"
echo "prefix : $PREFIX"
echo "table  : $TABLE"
echo "samples: Super-population=$SUPERPOP  ancient=$ANCIENT   PCs: $NPCS   linkage: $LINKAGE   clusters: ${N_CLUSTERS:-none}"

python3 - "$PREFIX" "$TABLE" "$NPCS" "$SUPERPOP" "$ANCIENT" "$LINKAGE" "$N_CLUSTERS" "$SCALE_PCS" <<'PY'
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, cophenet
from scipy.spatial.distance import pdist

prefix, table_path, npcs, superpop, ancient, method, n_clusters, scale = sys.argv[1:9]
npcs, ancient, scale = int(npcs), ancient == "1", scale == "1"
n_clusters = int(n_clusters) if n_clusters.strip() else None
if method not in ("ward", "average", "complete", "single", "centroid", "median", "weighted"):
    sys.exit(f"ERROR: unknown LINKAGE {method!r}")

# ---- fam / eigenvec ----------------------------------------------------------------------------
fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; using FID as the sample ID")
fam["sid"] = fam["FID"] if use_fid else fam["IID"]
with open(f"{prefix}.eigenvec") as fh:
    first = fh.readline().split()
if first and first[0].startswith("#"):
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", dtype=str)
    ev.columns = [c.lstrip("#") for c in ev.columns]
    pc_cols = [c for c in ev.columns if c.upper().startswith("PC")]
    if "FID" not in ev:
        ev["FID"] = "0"
else:
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", header=None, dtype=str)
    pc_cols = [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
    ev.columns = ["FID", "IID"] + pc_cols
ev[pc_cols] = ev[pc_cols].astype(float)
ev = ev.rename(columns={"FID" if use_fid else "IID": "sid"})
if ev["sid"].duplicated().any():
    sys.exit("ERROR: duplicate sample IDs in eigenvec")
ev = ev[ev["sid"].isin(set(fam["sid"]))]
if len(pc_cols) < npcs:
    print(f"NOTE: eigenvec has only {len(pc_cols)} PCs; using all of them")
    npcs = len(pc_cols)
pcs_used = pc_cols[:npcs]

# ---- mega table --------------------------------------------------------------------------------
if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass the .csv instead")
else:
    t = pd.read_csv(table_path, dtype=str)
for c in ("IID", "Population", "Super-population", "Undergone Imputation"):
    if c not in t.columns:
        sys.exit(f"ERROR: table lacks column {c!r}; columns are: {', '.join(t.columns)}")
t = t.drop_duplicates("IID")
df = ev.drop(columns=["IID"], errors="ignore").merge(t, left_on="sid", right_on="IID", how="inner")
print(f"{len(ev)} samples in the PCA, {len(df)} of them in the mega table")

# ---- choose the samples --------------------------------------------------------------------------
is_anc = df["Undergone Imputation"].str.strip().str.lower() == "yes"
if superpop.lower() != "all":
    df = df[df["Super-population"].str.lower() == superpop.lower()]
    is_anc = df["Undergone Imputation"].str.strip().str.lower() == "yes"
if not ancient:
    df = df[~is_anc]
src = os.environ.get("IRAN_SOURCE", "Iranian cohort")
if superpop.lower() == "iran" and src and "Source Publication(s)" in df.columns:
    keep = df["Source Publication(s)"].fillna("").str.contains(src, case=False, regex=False) | \
           (df["Undergone Imputation"].str.strip().str.lower() == "yes")
    if keep.any() and (~keep).any():
        print(f"leaving out {int((~keep).sum())} other Iran-labelled reference samples (Source lacks {src!r}): "
              + ", ".join(df.loc[~keep, "sid"]))
    if keep.any():
        df = df[keep]
df = df.reset_index(drop=True)
n = len(df)
if n < 3:
    sys.exit(f"ERROR: only {n} samples selected; need at least 3 (check SUPERPOP / ANCIENT)")
if n > 600:
    print(f"WARNING: {n} leaves - labels will be tiny; consider a subset (SUPERPOP / ANCIENT)")
df["pop"] = df["Population"].fillna("Population not recorded")

# ---- clustering --------------------------------------------------------------------------------
X = df[pcs_used].to_numpy(float)
if scale:
    X = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))
dists = pdist(X, metric="euclidean")
Z = linkage(dists, method=method)
coph = cophenet(Z, dists)[0]
print(f"{n} samples, {npcs} PCs, {method} linkage, cophenetic correlation {coph:.3f}")

# ---- dendrogram --------------------------------------------------------------------------------
pops = sorted(df["pop"].unique())
cmap = np.vstack([plt.cm.tab20(np.arange(20)), plt.cm.tab20b(np.arange(20))])
col_of = {p: matplotlib.colors.to_hex(cmap[i % len(cmap)]) for i, p in enumerate(pops)}

cut_h = None
cluster = np.ones(n, int)
if n_clusters and n_clusters > 1:
    k = min(n_clusters, n)
    cluster = fcluster(Z, k, criterion="maxclust")
    cut_h = (Z[n - k - 1, 2] + Z[n - k, 2]) / 2 if k < n else None

fs = float(np.clip(560 / n, 4.5, 9))                      # label font size shrinks with the number of leaves
width = float(np.clip(0.22 * n, 12, 80))
fig = plt.figure(figsize=(width, 10))
ax = fig.add_axes([0.04, 0.30, 0.80 if len(pops) <= 40 else 0.94, 0.66])
res = dendrogram(Z, ax=ax, no_labels=True, color_threshold=cut_h if cut_h else 0,
                 above_threshold_color="#444444", link_color_func=None)
order = res["leaves"]
xs = 5 + 10 * np.arange(n)
ax.set_xlim(0, 10 * n)
ax.set_ylabel(f"{method} linkage distance ({npcs} PCs{', standardised' if scale else ''})")
ax.set_title(f"{prefix.split('/')[-1]}  |  hierarchical clustering, n={n}  ({superpop}{'' if ancient else ', modern only'})")
for sp in ("top", "right", "bottom"):
    ax.spines[sp].set_visible(False)
ax.tick_params(axis="x", bottom=False, labelbottom=False)
if cut_h:
    ax.axhline(cut_h, color="red", ls="--", lw=0.8)
    ax.text(10 * n, cut_h, f" {int(cluster.max())} clusters", color="red", va="bottom", ha="right", fontsize=9)

ids = df["sid"].to_numpy()[order]
labs = df["pop"].to_numpy()[order]
colors = [col_of[p] for p in labs]
# text: a darkened version of the population colour (pastels are unreadable on white); strip keeps the full colour
dark = [tuple(0.6 * v for v in matplotlib.colors.to_rgb(c)) for c in colors]
# colour strip under the leaves + "<sample ID> | <Population>" text
ax.scatter(xs, [-0.012] * n, transform=ax.get_xaxis_transform(), marker="s", s=max(fs * 6, 14), c=colors,
           edgecolors="none", clip_on=False, zorder=5)
for x, i, p, c in zip(xs, ids, labs, dark):
    ax.text(x, -0.03, f"{i} | {p}", transform=ax.get_xaxis_transform(), rotation=90, ha="center", va="top",
            fontsize=fs, color=c, clip_on=False)

if len(pops) <= 40:
    h = [Line2D([], [], marker="s", ls="", color=col_of[p], markersize=8,
                label=f"{p} ({int((df['pop'] == p).sum())})") for p in pops]
    fig.legend(handles=h, loc="upper right", bbox_to_anchor=(0.995, 0.96), fontsize=8, frameon=False,
               title="Population", title_fontsize=9)
tag = superpop.lower().replace(" ", "_")
out = f"{prefix}.hclust_{tag}"
for ext in ("png", "pdf"):
    fig.savefig(f"{out}.{ext}", dpi=200)
plt.close(fig)
print("wrote", f"{out}.png", f"{out}.pdf")

res_df = pd.DataFrame({"leaf_order": np.arange(1, n + 1), "IID": ids, "Population": labs,
                       "cluster": cluster[order]})
res_df.to_csv(f"{out}.tsv", sep="\t", index=False)
print("wrote", f"{out}.tsv")
if n_clusters and n_clusters > 1:
    print("cluster sizes:", dict(zip(*np.unique(cluster, return_counts=True))))
    print(pd.crosstab(df["pop"], cluster).to_string())
PY
