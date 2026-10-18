#!/bin/bash

#SBATCH -c 4

# Save our arguments, then clear the positional parameters: `source`d scripts inherit them and
# conda's activate rejects extra arguments.
ARGS=("$@")
set --

#Activate your environments here
source ~/.bashrc
source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate /scratch/nme3923/conda-envs/glimpse2

# HDBSCAN clustering of samples on their top principal components. Density based: the number of
# clusters is not chosen up front, and samples that belong to no dense region are labelled noise (-1).
#
# Usage:   sbatch hdbscan.sh <prefix> [mega_sample_table.xlsx|.csv] [n_pcs=10]
#          bash   hdbscan.sh <prefix> [mega_sample_table.xlsx|.csv] [n_pcs=10]
#
#   <prefix>   expects  <prefix>.fam  <prefix>.eigenvec   (PLINK 1.9 or PLINK 2 eigenvec; with or
#              without a header line)
#   table      default: $MASTER_TABLE, else metadata/mega_sample_table.xlsx under $SLURM_SUBMIT_DIR
#              (or the current directory); the .csv next to it is used if openpyxl is missing
#   n_pcs      how many leading PCs to cluster on (default 10)
#
# Which samples (environment variables), same meaning as in cluster_pcs.sh:
#   SUPERPOP=all       Super-population to cluster (default all samples); e.g. SUPERPOP=Iran
#   ANCIENT=1          0 = leave out ancient samples (default: include them)
#   IRAN_SOURCE="Iranian cohort"   with SUPERPOP=Iran keep only the study's own Iranians
# Clustering:
#   MIN_CLUSTER_SIZE=10   smallest group of samples that may be called a cluster
#   MIN_SAMPLES=          core-point neighbourhood size (default = MIN_CLUSTER_SIZE); larger = more noise
#   SELECTION=eom         eom (default, fewer larger clusters) | leaf (many small, homogeneous clusters)
#   EPSILON=0             clusters closer than this distance are merged (0 = off)
#   SCALE_PCS=0           1 = standardise every PC first (default 0: keep the PCA scale so PC1 counts most)
#
# Outputs (next to the prefix):
#   <prefix>.hdbscan.tsv           IID, cluster (-1 = noise), probability, PCs used, Population, Super-population, ...
#   <prefix>.hdbscan_summary.tsv   per cluster: size, populations it contains, purity of the biggest one
#   <prefix>.hdbscan_PC1_PC2.{png,pdf}   (and PC3_PC4 when 4+ PCs are used)  samples coloured by cluster, noise grey
#
# Needs in the conda env: python3, numpy, pandas, matplotlib and scikit-learn >= 1.3 (or the `hdbscan` package).

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

for f in "${PREFIX}.fam" "${PREFIX}.eigenvec" "$TABLE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done
python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas and matplotlib" >&2; exit 1; }
python3 -c "from sklearn.cluster import HDBSCAN" 2>/dev/null || python3 -c "import hdbscan" 2>/dev/null || {
    echo "ERROR: needs scikit-learn >= 1.3 (conda install -c conda-forge scikit-learn) or the hdbscan package" >&2; exit 1; }

SUPERPOP="${SUPERPOP:-all}"; ANCIENT="${ANCIENT:-1}"
MIN_CLUSTER_SIZE="${MIN_CLUSTER_SIZE:-10}"; MIN_SAMPLES="${MIN_SAMPLES:-}"
SELECTION="${SELECTION:-eom}"; EPSILON="${EPSILON:-0}"; SCALE_PCS="${SCALE_PCS:-0}"
echo "prefix : $PREFIX"
echo "table  : $TABLE"
echo "samples: Super-population=$SUPERPOP  ancient=$ANCIENT   PCs: $NPCS"
echo "hdbscan: min_cluster_size=$MIN_CLUSTER_SIZE  min_samples=${MIN_SAMPLES:-$MIN_CLUSTER_SIZE}  selection=$SELECTION  epsilon=$EPSILON  scale=$SCALE_PCS"

python3 - "$PREFIX" "$TABLE" "$NPCS" "$SUPERPOP" "$ANCIENT" "$MIN_CLUSTER_SIZE" "$MIN_SAMPLES" "$SELECTION" "$EPSILON" "$SCALE_PCS" <<'PY'
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

prefix, table_path, npcs, superpop, ancient = sys.argv[1:6]
min_cluster_size, min_samples, selection, epsilon, scale = sys.argv[6:11]
npcs, ancient, scale = int(npcs), ancient == "1", scale == "1"
min_cluster_size = int(min_cluster_size)
min_samples = int(min_samples) if min_samples.strip() else None
epsilon = float(epsilon)
if selection not in ("eom", "leaf"):
    sys.exit(f"ERROR: SELECTION must be eom or leaf (got {selection!r})")

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

# ---- choose the samples (same rules as cluster_pcs.sh) -----------------------------------------
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
    if keep.any():
        df = df[keep]
df = df.reset_index(drop=True)
n = len(df)
if n < max(3, min_cluster_size):
    sys.exit(f"ERROR: only {n} samples selected; need at least MIN_CLUSTER_SIZE={min_cluster_size} "
             "(check SUPERPOP / ANCIENT)")
df["pop"] = df["Population"].fillna("Population not recorded")
df["superpop"] = df["Super-population"].fillna("Super-population not recorded")

# ---- HDBSCAN -----------------------------------------------------------------------------------
X = df[pcs_used].to_numpy(float)
if scale:
    X = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))
kw = dict(min_cluster_size=min_cluster_size, min_samples=min_samples,
          cluster_selection_method=selection, cluster_selection_epsilon=epsilon)
try:
    from sklearn.cluster import HDBSCAN
    model = HDBSCAN(**kw)
    backend = "scikit-learn"
except ImportError:
    import hdbscan
    model = hdbscan.HDBSCAN(**kw)
    backend = "hdbscan"
model.fit(X)
df["cluster"] = model.labels_.astype(int)
df["probability"] = np.round(model.probabilities_, 4)
n_clu = int(df["cluster"].max() + 1) if (df["cluster"] >= 0).any() else 0
n_noise = int((df["cluster"] < 0).sum())
print(f"{backend}: {n} samples, {npcs} PCs -> {n_clu} clusters, {n_noise} noise samples ({100 * n_noise / n:.1f}%)")
if n_clu == 0:
    print("WARNING: no clusters found; try a smaller MIN_CLUSTER_SIZE / MIN_SAMPLES or SELECTION=leaf")

# ---- outputs: per-sample table + cluster summary -------------------------------------------------
keep_cols = ["sid", "cluster", "probability"] + pcs_used + ["Population", "Super-population"]
keep_cols += [c for c in ("Source Publication(s)", "Empirical cluster (Study Pop only)", "Undergone Imputation")
              if c in df.columns]
out = df[keep_cols].rename(columns={"sid": "IID"}).sort_values(["cluster", "probability"], ascending=[True, False])
out.to_csv(f"{prefix}.hdbscan.tsv", sep="\t", index=False)

rows = []
for cl, g in df.groupby("cluster"):
    pc = g["pop"].value_counts()
    sc = g["superpop"].value_counts()
    rows.append({"cluster": cl, "n": len(g), "n_populations": len(pc),
                 "top_population": pc.index[0], "top_population_fraction": round(pc.iloc[0] / len(g), 3),
                 "top_super_population": sc.index[0],
                 "mean_probability": round(g["probability"].mean(), 3),
                 "populations": "; ".join(f"{p} ({k})" for p, k in pc.head(8).items())
                                + ("; ..." if len(pc) > 8 else "")})
summary = pd.DataFrame(rows).sort_values("cluster")
summary.to_csv(f"{prefix}.hdbscan_summary.tsv", sep="\t", index=False)
print(summary[["cluster", "n", "n_populations", "top_population", "top_population_fraction"]]
      .to_string(index=False))

# ---- plots: PC1 vs PC2 (and PC3 vs PC4), coloured by cluster, noise grey -------------------------
cmap = np.vstack([plt.cm.tab20(np.arange(20)), plt.cm.tab20b(np.arange(20))])
clusters = sorted(c for c in df["cluster"].unique() if c >= 0)
col = {c: cmap[i % len(cmap)] for i, c in enumerate(clusters)}

def plot(a, b):
    fig, ax = plt.subplots(figsize=(11, 8))
    nz = df[df["cluster"] < 0]
    if len(nz):
        ax.scatter(nz[a], nz[b], s=10, c="#9a9a9a", marker="x", linewidths=0.8, alpha=0.7,
                   label=f"noise (n={len(nz)})", zorder=1)
    for c in clusters:
        g = df[df["cluster"] == c]
        top = g["pop"].value_counts().index[0]
        ax.scatter(g[a], g[b], s=14, color=col[c], alpha=0.9, linewidths=0,
                   label=f"{c}: {top} (n={len(g)})", zorder=2)
    ax.set_xlabel(a)
    ax.set_ylabel(b)
    ax.set_title(f"{os.path.basename(prefix)}  |  HDBSCAN on {npcs} PCs: {n_clu} clusters, {n_noise} noise "
                 f"(min_cluster_size={min_cluster_size}, {selection})")
    ax.legend(markerscale=1.6, fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left",
              ncol=1 if len(clusters) <= 40 else 2, title="cluster: most common population")
    fig.tight_layout()
    for ext, kwargs in (("png", {"dpi": 200}), ("pdf", {})):
        fig.savefig(f"{prefix}.hdbscan_{a}_{b}.{ext}", **kwargs)
    plt.close(fig)

plot(pcs_used[0], pcs_used[1])
if npcs >= 4:
    plot(pcs_used[2], pcs_used[3])
print(f"wrote {prefix}.hdbscan.tsv, {prefix}.hdbscan_summary.tsv and the {prefix}.hdbscan_PC*.png/pdf plots")
PY
