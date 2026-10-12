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

# Plot a PCA coloured by the labels in the master sample table.
#
# Usage:   sbatch plot_pca.sh <prefix> [master_sample_table.tsv|.xlsx]
#          bash   plot_pca.sh <prefix> [master_sample_table.tsv|.xlsx]
#
#   <prefix>  expects  <prefix>.fam  <prefix>.eigenvec  <prefix>.eigenval
#             (PLINK 1.9 or PLINK 2 eigenvec; with or without a header line)
#   table     .tsv (preferred) or .xlsx (needs openpyxl; reads sheet "sample_table"); defaults to $MASTER_TABLE, else metadata/master_sample_table.tsv under
#             $SLURM_SUBMIT_DIR (or the current directory)
#
# Colour key:  ancient -> analysis_label_published (falls back to empirical_cluster when the study
#                         gives no analysis label, i.e. Mathieson 2015 samples)
#              modern  -> self_reported_label
# Colour/marker per label is fixed by the label's rank among ALL labels in the master table, so a
# label looks the same in every PCA you make from the same table.
#
# Outputs (next to the prefix):  <prefix>.pca_PC1_PC2.{png,pdf}, <prefix>.pca_PC3_PC4.{png,pdf},
#                                <prefix>.pca_labels.tsv (label, colour, marker, n, group)

set -euo pipefail

PREFIX="${ARGS[0]:-}"
if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <prefix> [master_sample_table.tsv|.xlsx]" >&2
    exit 1
fi
TABLE="${ARGS[1]:-${MASTER_TABLE:-${SLURM_SUBMIT_DIR:-$PWD}/metadata/master_sample_table.tsv}}"

for f in "${PREFIX}.fam" "${PREFIX}.eigenvec" "${PREFIX}.eigenval" "$TABLE"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas and matplotlib" >&2; exit 1; }

echo "prefix : $PREFIX"
echo "table  : $TABLE"

python3 - "$PREFIX" "$TABLE" <<'PY'
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

prefix, table_path = sys.argv[1], sys.argv[2]

# ---- read fam / eigenvec / eigenval ------------------------------------------------------------
fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
# Un-repaired fam where every IID is one constant string (ID stored in FID): fall back to FID.
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; using FID as the sample ID")
fam["sid"] = fam["FID"] if use_fid else fam["IID"]

with open(f"{prefix}.eigenvec") as fh:
    first = fh.readline().split()
has_header = first and first[0].startswith("#")
if has_header:                                   # PLINK 2: "#FID IID PC1..." or "#IID PC1..."
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", dtype=str)
    ev.columns = [c.lstrip("#") for c in ev.columns]
    pc_cols = [c for c in ev.columns if c.upper().startswith("PC")]
    id_cols = [c for c in ev.columns if c not in pc_cols]
    ev["FID"] = ev["FID"] if "FID" in ev else "0"
else:                                            # PLINK 1.9: FID IID PC1 PC2 ...
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", header=None, dtype=str)
    pc_cols = [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
    ev.columns = ["FID", "IID"] + pc_cols
ev[pc_cols] = ev[pc_cols].astype(float)

key = "FID" if use_fid else "IID"
ev = ev.rename(columns={key: "sid"})
if len(pc_cols) < 2:
    sys.exit("ERROR: eigenvec has fewer than 2 PCs")
eigval = np.loadtxt(f"{prefix}.eigenval", ndmin=1)
pct = 100 * eigval / eigval.sum()

# ---- fam <-> eigenvec consistency --------------------------------------------------------------
fam_ids, ev_ids = set(fam["sid"]), set(ev["sid"])
if fam_ids != ev_ids:
    print(f"WARNING: fam has {len(fam_ids)} samples, eigenvec {len(ev_ids)}; "
          f"{len(fam_ids - ev_ids)} only in fam, {len(ev_ids - fam_ids)} only in eigenvec")
if ev["sid"].duplicated().any():
    sys.exit("ERROR: duplicate sample IDs in eigenvec; cannot join to the sample table unambiguously")
ev = ev[ev["sid"].isin(fam_ids)]               # fam defines the sample set

# ---- labels from the master table --------------------------------------------------------------
if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="sample_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass master_sample_table.tsv instead")
else:
    t = pd.read_csv(table_path, sep="\t", dtype=str)
# Label for every table row: ancient -> analysis_label_published, modern -> self_reported_label.
def pick(r):
    if r.sample_type == "ancient":
        for c in ("analysis_label_published", "empirical_cluster"):
            if isinstance(getattr(r, c), str):
                return getattr(r, c), c
        return "Ancient_unlabelled", "none"
    if isinstance(r.self_reported_label, str):
        return r.self_reported_label, "self_reported_label"
    return "Modern_unreported", "none"

t[["label", "label_source"]] = t.apply(lambda r: pd.Series(pick(r)), axis=1)
t = t.drop_duplicates("sample_id")

df = ev.merge(t[["sample_id", "sample_type", "label", "label_source", "study", "inclusion_status"]],
              left_on="sid", right_on="sample_id", how="left")
unk = df["sample_type"].isna()
if unk.any():
    print(f"WARNING: {unk.sum()} samples not in the master table (plotted grey as 'Not_in_table'): "
          f"{', '.join(df.loc[unk, 'sid'].head(10))}{' ...' if unk.sum() > 10 else ''}")
df.loc[unk, ["sample_type", "label", "label_source"]] = ["unknown", "Not_in_table", "none"]
fb = (df["label_source"] == "empirical_cluster").sum()
if fb:
    print(f"NOTE: {fb} ancient samples have no analysis_label_published; used empirical_cluster instead")
df["group"] = df["sample_type"].map({"ancient": "Ancient", "modern": "Modern"}).fillna("Unknown")

# ---- stable colour / marker per label ----------------------------------------------------------
# Rank over ALL labels in the table (not just this PCA) so a label keeps its style across plots.
lab_tab = pd.concat([t.loc[t.sample_type == "ancient", ["label"]].assign(g="Ancient"),
                     t.loc[t.sample_type == "modern", ["label"]].assign(g="Modern")])
all_lab = sorted(set(map(tuple, lab_tab[["g", "label"]].drop_duplicates().values.tolist())))
rank = {k: i for i, k in enumerate(all_lab)}
n_col = 48
cols = np.vstack([plt.cm.tab20(np.arange(20)), plt.cm.tab20b(np.arange(20)),
                  plt.cm.Dark2(np.arange(8))])[:n_col]
markers = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "h", "p", "*"]
def style(g, lab):
    i = rank.get((g, lab), len(rank))
    return cols[i % n_col], markers[(i // n_col) % len(markers)]

# ---- plot --------------------------------------------------------------------------------------
rows = []
def plot_pair(a, b, out):
    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_axes([0.05, 0.07, 0.40, 0.88])
    for grp, z, size, edge in (("Modern", 1, 14, "none"), ("Unknown", 2, 14, "none"), ("Ancient", 3, 34, "black")):
        sub = df[df["group"] == grp]
        for lab, s in sub.groupby("label", sort=True):
            if grp == "Unknown":
                c, m = np.array([0.5, 0.5, 0.5, 1]), "x"
            else:
                c, m = style(grp, lab)
            ax.scatter(s[a], s[b], s=size, c=[c], marker=m, edgecolors=edge, linewidths=0.4,
                       alpha=0.75 if grp == "Modern" else 0.95, zorder=z)
    ia, ib = pc_cols.index(a), pc_cols.index(b)
    ax.set_xlabel(f"{a} ({pct[ia]:.2f}% of listed eigenvalues)" if ia < len(pct) else a)
    ax.set_ylabel(f"{b} ({pct[ib]:.2f}% of listed eigenvalues)" if ib < len(pct) else b)
    ax.set_title(f"{prefix.split('/')[-1]}  |  n={len(df)}  "
                 f"(ancient {int((df.group=='Ancient').sum())}, modern {int((df.group=='Modern').sum())})")
    ax.grid(alpha=0.2)

    # legends: one block per group, own columns, small font
    for grp, x0, w in (("Ancient", 0.46, 0.32), ("Modern", 0.78, 0.21)):
        labs = sorted(df.loc[df.group == grp, "label"].unique())
        if not labs:
            continue
        handles = []
        for lab in labs:
            c, m = style(grp, lab)
            n = int(((df.group == grp) & (df.label == lab)).sum())
            handles.append(Line2D([], [], marker=m, ls="", color=c, markersize=5,
                                  markeredgecolor="black" if grp == "Ancient" else "none",
                                  label=f"{lab if len(lab) <= 42 else lab[:41] + '~'} ({n})"))
        lax = fig.add_axes([x0, 0.02, w, 0.94]); lax.axis("off")
        leg = lax.legend(handles=handles, loc="upper left", ncol=int(np.ceil(len(labs) / 120)), fontsize=float(np.clip(700 / max(len(labs), 1), 5, 9)),
                   frameon=False, title=f"{grp} - " + ("analysis_label_published" if grp == "Ancient"
                                                      else "self_reported_label"),
                   title_fontsize=8, labelspacing=0.22, handletextpad=0.2, columnspacing=0.8, borderaxespad=0)
        leg._legend_box.align = "left"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    plt.close(fig)
    print("wrote", f"{out}.png", f"{out}.pdf")

plot_pair(pc_cols[0], pc_cols[1], f"{prefix}.pca_PC1_PC2")
if len(pc_cols) >= 4:
    plot_pair(pc_cols[2], pc_cols[3], f"{prefix}.pca_PC3_PC4")

# label key for reproducibility
for (grp, lab), s in df.groupby(["group", "label"]):
    c, m = style(grp, lab)
    rows.append((grp, lab, matplotlib.colors.to_hex(c), m, len(s)))
pd.DataFrame(rows, columns=["group", "label", "colour", "marker", "n"]).to_csv(
    f"{prefix}.pca_labels.tsv", sep="\t", index=False)
print("wrote", f"{prefix}.pca_labels.tsv")
PY
