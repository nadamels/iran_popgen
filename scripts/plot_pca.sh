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

# Plot a PCA from the mega sample table.
#   COLOUR = Super-population (colours from Plotting/pca_colors.txt)
#   SIZE   = age of the sample (older = bigger), ancient samples only
#   Modern samples = tiny dots (coloured by super-population); no marker shapes for now.
#   Modern IRANIANS (Super-population Iran, the target population) = large outlined SQUARES drawn
#   on top of everything else.
#
# Usage:   sbatch plot_pca.sh <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [PC1_range] [PC2_range]
#          bash   plot_pca.sh <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [PC1_range] [PC2_range]
#
#   Zoom example (PC1 from -0.01 to 0.01, PC2 from 0.0 to 0.02):
#          sbatch plot_pca.sh <prefix> <table> <colors> -0.01,0.01 0.0,0.02
#   A range is "min,max" (or "min:max"); use "auto" (or "-") to leave that axis unzoomed.
#   Defaults for the first three arguments can be skipped with "" (empty string).
#
#   <prefix>  expects  <prefix>.fam  <prefix>.eigenvec  <prefix>.eigenval
#             (PLINK 1.9 or PLINK 2 eigenvec; with or without a header line)
#   table     .xlsx (sheet "mega_table"; needs openpyxl) or .csv. Default: $MASTER_TABLE, else
#             metadata/mega_sample_table.xlsx under $SLURM_SUBMIT_DIR (or the current directory),
#             falling back to the .csv next to it if openpyxl is missing
#   colours   "name: #hex" lines. Default: $COLORS_FILE, else Plotting/pca_colors.txt (same base dir)
#
# Age: parsed from 'Archeological Date' (e.g. "5898-5531 calBCE", "600-800 CE", "1500-800 BCE")
#      as the range midpoint in years before present (1950). Unparseable dates are drawn at the
#      smallest ancient size and counted in the log.
#
#   When a PC1/PC2 range is given, ONLY the zoomed PC1 vs PC2 plot is made, written to
#   <prefix>.pca_PC1_PC2.zoom_x<min>_<max>_y<min>_<max>.{png,pdf} (the full plots are not overwritten).
#
# Interactive version: every run also writes <prefix>.pca_interactive.html (one self-contained file,
#   no internet needed): hover a point for its metadata, click a Super-population in the legend to
#   hide/show it, switch the PCs on each axis, wheel = zoom, drag = pan, double-click = reset.
#   With a zoom range the HTML opens at that range (button "Full extent" shows everything) and is
#   named <prefix>.pca_interactive.zoom_x<min>_<max>_y<min>_<max>.html.
#
# Outputs (next to the prefix):  <prefix>.pca_PC1_PC2.{png,pdf}, <prefix>.pca_PC3_PC4.{png,pdf},
#                                <prefix>.pca_samples.tsv (IID, super-population, age BP, colour)

set -euo pipefail

PREFIX="${ARGS[0]:-}"
if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [PC1_min,max] [PC2_min,max]" >&2
    exit 1
fi
BASE="${SLURM_SUBMIT_DIR:-$PWD}"
TABLE="${ARGS[1]:-${MASTER_TABLE:-$BASE/metadata/mega_sample_table.xlsx}}"
COLORS="${ARGS[2]:-${COLORS_FILE:-$BASE/Plotting/pca_colors.txt}}"
PC1_RANGE="${ARGS[3]:-auto}"
PC2_RANGE="${ARGS[4]:-auto}"

if [[ "$TABLE" == *.xlsx || "$TABLE" == *.xls ]] && ! python3 -c "import openpyxl" 2>/dev/null; then
    alt="${TABLE%.*}.csv"
    if [[ -s "$alt" ]]; then
        echo "NOTE: openpyxl not available; using $alt instead of $TABLE" >&2
        TABLE="$alt"
    fi
fi

for f in "${PREFIX}.fam" "${PREFIX}.eigenvec" "${PREFIX}.eigenval" "$TABLE" "$COLORS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy, pandas and matplotlib" >&2; exit 1; }

echo "prefix : $PREFIX"
echo "table  : $TABLE"
echo "colours: $COLORS"
echo "PC1 range: $PC1_RANGE   PC2 range: $PC2_RANGE"

python3 - "$PREFIX" "$TABLE" "$COLORS" "$PC1_RANGE" "$PC2_RANGE" <<'PY'
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

prefix, table_path, colors_path = sys.argv[1], sys.argv[2], sys.argv[3]

def parse_range(txt, name):
    """'min,max' or 'min:max' -> (min, max); 'auto' / '-' / '' -> None"""
    if txt.strip().lower() in ("", "auto", "-", "none"):
        return None
    parts = [x for x in re.split(r"[,:]", txt.strip()) if x != ""]
    try:
        lo, hi = (float(x) for x in parts)
    except ValueError:
        sys.exit(f"ERROR: {name} range {txt!r} must look like min,max (e.g. -0.01,0.01) or 'auto'")
    if not lo < hi:
        sys.exit(f"ERROR: {name} range {txt!r}: min must be smaller than max")
    return lo, hi
xr, yr = parse_range(sys.argv[4], "PC1"), parse_range(sys.argv[5], "PC2")
zoom = xr is not None or yr is not None

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

# ---- mega table: super-population, age ---------------------------------------------------------
if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass the .csv instead")
else:
    t = pd.read_csv(table_path, dtype=str)
need = ["IID", "Super-population", "Archeological Date", "Undergone Imputation"]
miss_cols = [c for c in need if c not in t.columns]
if miss_cols:
    sys.exit(f"ERROR: table lacks column(s) {miss_cols}")
t = t.drop_duplicates("IID")

DATE = re.compile(r"^\s*(?:(?P<a>\d+)\s*-\s*(?P<b>\d+)|(?P<c>\d+))\s*(?P<era>cal\s*BCE|BCE|cal\s*CE|CE|AD|cal\s*AD|calBP|BP)", re.I)
def age_bp(txt):
    """Midpoint of the date range in years before present (1950)."""
    m = DATE.match(txt) if isinstance(txt, str) else None
    if not m:
        return np.nan
    lo, hi = (m["a"], m["b"]) if m["a"] else (m["c"], m["c"])
    mid = (float(lo) + float(hi)) / 2
    era = re.sub(r"\s+", "", m["era"]).lower()
    if era.endswith("bce"):
        return mid + 1950
    if era.endswith("bp"):
        return mid
    return 1950 - mid                                   # CE / AD

t["is_ancient"] = t["Undergone Imputation"].str.strip().str.lower() == "yes"
t["age"] = [age_bp(d) if a else np.nan for d, a in zip(t["Archeological Date"], t["is_ancient"])]

df = ev.drop(columns=["IID"], errors="ignore").merge(t, left_on="sid", right_on="IID", how="left")
unk = df["IID"].isna()
if unk.any():
    print(f"WARNING: {unk.sum()} samples not in the mega table (plotted grey 'x'): "
          f"{', '.join(df.loc[unk, 'sid'].head(10))}{' ...' if unk.sum() > 10 else ''}")
df["is_ancient"] = df["is_ancient"].fillna(False).astype(bool)
df["sp"] = df["Super-population"].fillna("Unassigned")
n_bad = int((df["is_ancient"] & df["age"].isna()).sum())
if n_bad:
    print(f"NOTE: {n_bad} ancient samples have no parseable date; drawn at the smallest ancient size")

# ---- colours -----------------------------------------------------------------------------------
def norm(s):
    return re.sub(r"[^a-z]", "", s.lower())
file_cols = {}
with open(colors_path) as fh:
    for line in fh:
        if ":" in line:
            k, v = line.split(":", 1)
            if v.strip().startswith("#"):
                file_cols[norm(k)] = v.strip()
# table super-population name -> key in the colour file (names differ slightly)
ALIAS = {"iran": "iraniansamplesofthisstudy", "arabia": "arabian", "levant": "levantine",
         "northafrica": "africa"}
# not listed in pca_colors.txt: fallback colours (edit the file to override; add a line "Name: #hex")
FALLBACK = {"anatolia": "#8c564b", "americas": "#2ca02c", "oceania": "#17becf", "unassigned": "#7f7f7f"}
def colour_of(sp):
    k = norm(sp)
    for key in (k, ALIAS.get(k)):
        if key and key in file_cols:
            return file_cols[key], True
    return FALLBACK.get(k, "#7f7f7f"), False
sps = sorted(df["sp"].unique())
col_of, fallback_used = {}, []
for sp in sps:
    col_of[sp], ok = colour_of(sp)
    if not ok:
        fallback_used.append(sp)
if fallback_used:
    print("NOTE: no colour in pca_colors.txt for: " + ", ".join(fallback_used) +
          " (fallback colours used; add 'Name: #hex' lines to set them)")

# ---- size from age -----------------------------------------------------------------------------
S_MIN, S_MAX = 12.0, 130.0                      # marker area (pt^2) for the youngest / oldest ancient sample
AGE_MAX = float(np.nanmax(df.loc[df.is_ancient, "age"])) if df.is_ancient.any() else 1.0
AGE_MAX = max(AGE_MAX, 1.0)
def size_of(age):
    a = np.where(np.isfinite(age), np.clip(age, 0, AGE_MAX), 0.0)
    return S_MIN + (S_MAX - S_MIN) * a / AGE_MAX
df["size"] = size_of(df["age"].to_numpy(dtype=float))

# ---- plot --------------------------------------------------------------------------------------
IRAN_SIZE = 45.0                                # marker area (pt^2) for modern Iranians (target population)
anc = df[df.is_ancient]
mod_all = df[~df.is_ancient & df["IID"].notna()]
is_iran = mod_all["sp"].str.lower() == "iran"
mod, mod_iran = mod_all[~is_iran], mod_all[is_iran]
def plot_pair(a, b, out, xlim=None, ylim=None):
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_axes([0.06, 0.07, 0.66, 0.88])
    # modern: tiny dots, bottom layer
    for sp, s in mod.groupby("sp"):
        ax.scatter(s[a], s[b], s=3, c=[col_of[sp]], marker="o", edgecolors="none", alpha=0.7, zorder=1)
    u = df[df["IID"].isna()]
    if len(u):
        ax.scatter(u[a], u[b], s=14, c="grey", marker="x", zorder=2)
    # ancient: colour = super-population, size = age; younger on top so big old points stay visible
    s = anc.sort_values("age", ascending=False, na_position="last")
    ax.scatter(s[a], s[b], s=s["size"], c=[col_of[x] for x in s["sp"]], marker="o",
               edgecolors="black", linewidths=0.5, alpha=0.9, zorder=3)
    # modern Iranians: large, on top of ancient points too
    if len(mod_iran):
        ax.scatter(mod_iran[a], mod_iran[b], s=IRAN_SIZE, c=[col_of[x] for x in mod_iran["sp"]], marker="s",
                   edgecolors="black", linewidths=0.8, alpha=0.95, zorder=5)
    ia, ib = pc_cols.index(a), pc_cols.index(b)
    ax.set_xlabel(f"{a} ({pct[ia]:.2f}% of listed eigenvalues)" if ia < len(pct) else a)
    ax.set_ylabel(f"{b} ({pct[ib]:.2f}% of listed eigenvalues)" if ib < len(pct) else b)
    title = f"{prefix.split('/')[-1]}  |  n={len(df)}  (ancient {len(anc)}, modern {len(mod) + len(mod_iran)}, of which Iranian {len(mod_iran)})"
    if xlim is not None or ylim is not None:
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        sel = np.ones(len(df), bool)
        if xlim is not None:
            sel &= ((df[a] >= xlim[0]) & (df[a] <= xlim[1])).to_numpy()
        if ylim is not None:
            sel &= ((df[b] >= ylim[0]) & (df[b] <= ylim[1])).to_numpy()
        title += f"\nZOOM: {a} {xlim if xlim is not None else 'auto'}, {b} {ylim if ylim is not None else 'auto'}  -  {int(sel.sum())} samples in view"
        print(f"samples in the zoom window: {int(sel.sum())} of {len(df)}")
    ax.set_title(title)
    ax.grid(alpha=0.2)

    # legend 1: colour = super-population (ancient / modern counts)
    h1 = []
    for sp in sorted(sps, key=lambda x: (x == "Unassigned", x)):
        na, nm = int(((anc.sp == sp)).sum()), int(((mod_all.sp == sp)).sum())
        h1.append(Line2D([], [], marker="o", ls="", color=col_of[sp], markersize=8,
                         markeredgecolor="black", markeredgewidth=0.4,
                         label=f"{sp} (ancient {na}, modern {nm})"))
    # legend 2: size = age
    ticks = [t_ for t_ in (500, 2000, 4000, 6000, 8000, 10000, 12000) if t_ <= AGE_MAX * 1.05]
    h2 = [Line2D([], [], marker="o", ls="", color="0.75", markeredgecolor="black", markeredgewidth=0.5,
                 markersize=float(np.sqrt(size_of(np.array([t_]))[0])), label=f"{t_:,} years BP")
          for t_ in ticks]
    h2.append(Line2D([], [], marker="o", ls="", color="0.5", markersize=2, label="Modern (tiny dots)"))
    if len(mod_iran):
        h2.append(Line2D([], [], marker="s", ls="", color=col_of[mod_iran["sp"].iloc[0]], markeredgecolor="black",
                         markersize=float(np.sqrt(IRAN_SIZE)), label=f"Modern Iranians (n={len(mod_iran)})"))
    lax = fig.add_axes([0.73, 0.02, 0.26, 0.94]); lax.axis("off")
    l1 = lax.legend(handles=h1, loc="upper left", bbox_to_anchor=(0, 1), fontsize=8, frameon=False,
                    title="Colour: super-population", title_fontsize=9, labelspacing=0.45)
    l1._legend_box.align = "left"
    lax.add_artist(l1)
    l2 = lax.legend(handles=h2, loc="lower left", bbox_to_anchor=(0, 0), fontsize=8, frameon=False,
                    title="Size: age of ancient sample", title_fontsize=9, labelspacing=1.1)
    l2._legend_box.align = "left"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    plt.close(fig)
    print("wrote", f"{out}.png", f"{out}.pdf")

if zoom:
    def tag(r, pre):
        return f"{pre}{r[0]:g}_{r[1]:g}" if r is not None else f"{pre}auto"
    plot_pair(pc_cols[0], pc_cols[1], f"{prefix}.pca_PC1_PC2.zoom_{tag(xr, 'x')}_{tag(yr, 'y')}", xr, yr)
else:
    plot_pair(pc_cols[0], pc_cols[1], f"{prefix}.pca_PC1_PC2")
    if len(pc_cols) >= 4:
        plot_pair(pc_cols[2], pc_cols[3], f"{prefix}.pca_PC3_PC4")

# ---- interactive HTML ----------------------------------------------------------------------------
import json
def _clean(v):
    return None if (v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip() in ("", "nan")) else str(v)

skip_cols = {"IID", "Super-population", "Undergone Imputation", "is_ancient", "age", "sid", "sp"}
info_cols = [c for c in t.columns if c not in skip_cols]
npc = min(len(pc_cols), 10)
recs = []
for _, r in df.iterrows():
    info = {k: _clean(r.get(k)) for k in info_cols}
    recs.append({
        "id": r["sid"],
        "sp": r["sp"],
        "anc": bool(r["is_ancient"]),
        "iran": bool((not r["is_ancient"]) and str(r["sp"]).lower() == "iran" and not pd.isna(r.get("IID"))),
        "ok": bool(not pd.isna(r.get("IID"))),
        "age": None if pd.isna(r["age"]) else round(float(r["age"])),
        "p": [round(float(r[c]), 6) for c in pc_cols[:npc]],
        "info": {k: v for k, v in info.items() if v is not None},
    })
payload = {"title": f"{prefix.split('/')[-1]}  |  n={len(recs)}", "pcs": pc_cols[:npc],
           "pct": [round(float(x), 2) for x in pct[:npc]], "colors": {k: str(v) if isinstance(v, str) else matplotlib.colors.to_hex(v) for k, v in col_of.items()},
           "zoom": ([list(xr) if xr else None, list(yr) if yr else None] if zoom else None), "rows": recs}

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Interactive PCA</title>
<style>
:root{--bg:#fff;--fg:#1d1d1f;--mut:#6b6b70;--grid:#e6e6ea;--card:#fff;--bd:#cfcfd6}
@media (prefers-color-scheme:dark){:root{--bg:#16161a;--fg:#e8e8ec;--mut:#9a9aa3;--grid:#2c2c33;--card:#202026;--bd:#3a3a44}}
html,body{margin:0;height:100%;background:var(--bg);color:var(--fg);font:13px/1.35 system-ui,sans-serif}
#app{display:flex;height:100%}#main{flex:1;min-width:0;display:flex;flex-direction:column}
#bar{display:flex;gap:12px;align-items:center;padding:10px 16px;flex-wrap:wrap}#bar h1{font-size:15px;margin:0;font-weight:600}
select,button{background:var(--card);color:var(--fg);border:1px solid var(--bd);border-radius:6px;padding:3px 6px;font:inherit}
#plot{flex:1;min-height:0;position:relative}canvas{position:absolute;inset:0;width:100%;height:100%;cursor:grab}canvas:active{cursor:grabbing}
#side{width:300px;border-left:1px solid var(--bd);overflow:auto;padding:10px 12px}
#side h2{font-size:12px;margin:8px 0 6px;color:var(--mut);font-weight:600;text-transform:uppercase;letter-spacing:.04em}
.it{display:flex;align-items:center;gap:7px;padding:2px 0;cursor:pointer;user-select:none}.it.off{opacity:.35}
.sw{width:11px;height:11px;border-radius:50%;flex:none;border:1px solid #000}.nm{flex:1}.n{color:var(--mut)}
#tip{position:fixed;pointer-events:none;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px 10px;box-shadow:0 4px 18px rgba(0,0,0,.18);display:none;max-width:360px;z-index:5}
#tip b{font-size:14px}#tip .k{color:var(--mut)}#note{padding:0 16px 8px;color:var(--mut);font-size:12px}
@media (max-width:720px){#app{flex-direction:column}#side{width:auto;max-height:34%;border-left:0;border-top:1px solid var(--bd)}}
</style></head><body>
<div id="app"><div id="main">
 <div id="bar"><h1 id="ttl"></h1><label>x <select id="sx"></select></label><label>y <select id="sy"></select></label>
  <label><input type="checkbox" id="cA" checked> ancient</label><label><input type="checkbox" id="cM" checked> modern</label>
  <button id="rst">Reset view</button><button id="full">Full extent</button><button id="all">Show all</button><button id="none">Hide all</button></div>
 <div id="plot"><canvas id="cv"></canvas></div>
 <div id="note">Hover a point for details &middot; wheel = zoom, drag = pan, double-click = reset &middot; click a Super-population to hide/show it &middot; colour = Super-population, ancient: outlined, larger = older &middot; modern Iranians: squares</div>
</div>
<div id="side"><h2>Super-population</h2><div id="leg"></div></div></div>
<div id="tip"></div>
<script>
const D=__DATA__, cv=document.getElementById("cv"), ctx=cv.getContext("2d"), tip=document.getElementById("tip");
const rows=D.rows, M={l:62,r:18,t:12,b:46}; let ix=0, iy=1, view=null, plotted=[];
const hidden=new Set(); document.getElementById("ttl").textContent=D.title;
const ages=rows.filter(r=>r.anc&&r.age!==null).map(r=>r.age), AMAX=Math.max(1,...ages);
const rad=r=>r.iran?4.5:(r.anc?3+7*Math.min(Math.max(r.age||0,0),AMAX)/AMAX:1.7);
const order=[...rows.filter(r=>!r.anc&&!r.iran),...rows.filter(r=>r.anc).sort((a,b)=>(b.age||0)-(a.age||0)),...rows.filter(r=>r.iran)];
const col=sp=>D.colors[sp]||"#7f7f7f";
for(const id of ["sx","sy"]){const s=document.getElementById(id);D.pcs.forEach((c,i)=>{const o=document.createElement("option");o.value=i;o.textContent=c+" ("+D.pct[i]+"%)";s.appendChild(o)})}
sy.value=1; sx.onchange=()=>{ix=+sx.value;home()}; sy.onchange=()=>{iy=+sy.value;home()};
const sps=[...new Set(rows.map(r=>r.sp))].sort(); const cnt={};
rows.forEach(r=>{const c=cnt[r.sp]||(cnt[r.sp]={a:0,m:0});r.anc?c.a++:c.m++});
const leg=document.getElementById("leg");
sps.forEach(sp=>{const d=document.createElement("div");d.className="it";
 d.innerHTML='<span class="sw" style="background:'+col(sp)+'"></span><span class="nm">'+sp+'</span><span class="n">'+cnt[sp].a+' anc / '+cnt[sp].m+' mod</span>';
 d.onclick=()=>{hidden.has(sp)?hidden.delete(sp):hidden.add(sp);d.classList.toggle("off");draw()};leg.appendChild(d)});
all.onclick=()=>{hidden.clear();[...leg.children].forEach(c=>c.classList.remove("off"));draw()};
none.onclick=()=>{sps.forEach(s=>hidden.add(s));[...leg.children].forEach(c=>c.classList.add("off"));draw()};
cA.onchange=cM.onchange=draw;
function bounds(){let a=[1e9,-1e9,1e9,-1e9];rows.forEach(r=>{const x=r.p[ix],y=r.p[iy];if(x<a[0])a[0]=x;if(x>a[1])a[1]=x;if(y<a[2])a[2]=y;if(y>a[3])a[3]=y});
 const px=(a[1]-a[0])*.05||1,py=(a[3]-a[2])*.05||1;return [a[0]-px,a[1]+px,a[2]-py,a[3]+py]}
function home(){const b=bounds();if(D.zoom&&ix==0&&iy==1){const z=D.zoom;if(z[0]){b[0]=z[0][0];b[1]=z[0][1]}if(z[1]){b[2]=z[1][0];b[3]=z[1][1]}}view=b;draw()}
rst.onclick=home; cv.ondblclick=home; full.onclick=()=>{view=bounds();draw()};
function nice(lo,hi,n){const span=hi-lo,st=Math.pow(10,Math.floor(Math.log10(span/n)));const e=[1,2,5,10].find(m=>span/(st*m)<=n)*st;const o=[];for(let v=Math.ceil(lo/e)*e;v<=hi;v+=e)o.push(+v.toFixed(10));return o}
function vis(r){return !hidden.has(r.sp)&&(r.anc?cA.checked:cM.checked)}
function draw(){const b=cv.getBoundingClientRect(),dpr=window.devicePixelRatio||1,w=b.width,h=b.height;
 cv.width=w*dpr;cv.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
 const cs=getComputedStyle(document.documentElement),grid=cs.getPropertyValue("--grid"),mut=cs.getPropertyValue("--mut"),fg=cs.getPropertyValue("--fg");
 const sx_=v=>M.l+(v-view[0])/(view[1]-view[0])*(w-M.l-M.r),sy_=v=>h-M.b-(v-view[2])/(view[3]-view[2])*(h-M.t-M.b);
 ctx.font="11px system-ui";ctx.lineWidth=1;
 ctx.textAlign="center";nice(view[0],view[1],8).forEach(v=>{const x=sx_(v);ctx.strokeStyle=grid;ctx.beginPath();ctx.moveTo(x,M.t);ctx.lineTo(x,h-M.b);ctx.stroke();ctx.fillStyle=mut;ctx.fillText(v,x,h-M.b+16)});
 ctx.textAlign="right";nice(view[2],view[3],6).forEach(v=>{const y=sy_(v);ctx.strokeStyle=grid;ctx.beginPath();ctx.moveTo(M.l,y);ctx.lineTo(w-M.r,y);ctx.stroke();ctx.fillStyle=mut;ctx.fillText(v,M.l-8,y+4)});
 ctx.fillStyle=fg;ctx.font="12px system-ui";ctx.textAlign="center";ctx.fillText(D.pcs[ix]+" ("+D.pct[ix]+"%)",(M.l+w-M.r)/2,h-8);
 ctx.save();ctx.translate(14,(M.t+h-M.b)/2);ctx.rotate(-Math.PI/2);ctx.fillText(D.pcs[iy]+" ("+D.pct[iy]+"%)",0,0);ctx.restore();
 ctx.save();ctx.beginPath();ctx.rect(M.l,M.t,w-M.l-M.r,h-M.t-M.b);ctx.clip();plotted=[];
 for(const r of order){if(!vis(r))continue;const x=sx_(r.p[ix]),y=sy_(r.p[iy]);if(x<M.l-12||x>w-M.r+12||y<M.t-12||y>h-M.b+12)continue;
  const q=rad(r);plotted.push([r,x,y,q]);ctx.fillStyle=r.ok?col(r.sp):"#888";
  if(r.iran){ctx.globalAlpha=.95;ctx.fillRect(x-q,y-q,2*q,2*q);ctx.lineWidth=.8;ctx.strokeStyle="#000";ctx.strokeRect(x-q,y-q,2*q,2*q)}
  else{ctx.globalAlpha=r.anc?.9:.7;ctx.beginPath();ctx.arc(x,y,q,0,6.2832);ctx.fill();if(r.anc){ctx.lineWidth=.6;ctx.strokeStyle="#000";ctx.stroke()}}}
 ctx.restore();ctx.globalAlpha=1}
function pick(e){const b=cv.getBoundingClientRect(),mx=e.clientX-b.left,my=e.clientY-b.top;
 for(let i=plotted.length-1;i>=0;i--){const [r,x,y,q]=plotted[i],R=Math.max(q,4)+1;if(Math.abs(mx-x)<=R&&Math.abs(my-y)<=R&&(mx-x)**2+(my-y)**2<=R*R*1.6)return r}return null}
cv.addEventListener("mousemove",e=>{if(drag)return;const r=pick(e);if(!r){tip.style.display="none";return}
 const rows_=[["Sample ID","<b>"+r.id+"</b>"],["Super-population",r.sp],["Type",r.ok?(r.anc?"Ancient":"Modern"):"not in mega table"]];
 Object.keys(r.info).forEach(k=>{let v=r.info[k];if(/^(lat|lon)/i.test(k)&&!isNaN(+v))v=(+v).toFixed(4);rows_.push([k,v])});
 if(r.age!==null)rows_.push(["Age (yr BP)",r.age.toLocaleString()]);
 rows_.push([D.pcs[ix],r.p[ix].toFixed(4)],[D.pcs[iy],r.p[iy].toFixed(4)]);
 tip.innerHTML=rows_.map(x=>'<div><span class="k">'+x[0]+':</span> '+x[1]+'</div>').join("");tip.style.display="block";
 const w=tip.offsetWidth,h=tip.offsetHeight;let X=e.clientX+14,Y=e.clientY+14;if(X+w>innerWidth)X=e.clientX-w-14;if(Y+h>innerHeight)Y=e.clientY-h-14;tip.style.left=X+"px";tip.style.top=Y+"px"});
cv.addEventListener("mouseleave",()=>tip.style.display="none");
cv.addEventListener("wheel",e=>{e.preventDefault();const b=cv.getBoundingClientRect(),w=b.width,h=b.height;
 const fx=(e.clientX-b.left-M.l)/(w-M.l-M.r),fy=1-(e.clientY-b.top-M.t)/(h-M.t-M.b),k=e.deltaY<0?.85:1/.85;
 const cx=view[0]+fx*(view[1]-view[0]),cy=view[2]+fy*(view[3]-view[2]);view=[cx-(cx-view[0])*k,cx+(view[1]-cx)*k,cy-(cy-view[2])*k,cy+(view[3]-cy)*k];draw()},{passive:false});
let drag=null;cv.onmousedown=e=>{drag={x:e.clientX,y:e.clientY,v:view.slice()};tip.style.display="none"};addEventListener("mouseup",()=>drag=null);
addEventListener("mousemove",e=>{if(!drag)return;const b=cv.getBoundingClientRect(),dx=(e.clientX-drag.x)/(b.width-M.l-M.r)*(drag.v[1]-drag.v[0]),dy=(e.clientY-drag.y)/(b.height-M.t-M.b)*(drag.v[3]-drag.v[2]);view=[drag.v[0]-dx,drag.v[1]-dx,drag.v[2]+dy,drag.v[3]+dy];draw()});
addEventListener("resize",draw);home();
</script></body></html>"""

if zoom:
    html_out = f"{prefix}.pca_interactive.zoom_{tag(xr, 'x')}_{tag(yr, 'y')}.html"
else:
    html_out = f"{prefix}.pca_interactive.html"
with open(html_out, "w") as fh:
    fh.write(HTML.replace("__DATA__", json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")))
print("wrote", html_out)

out = df[["sid", "sp", "is_ancient", "age"]].copy()
out["colour"] = out["sp"].map(col_of)
out.rename(columns={"sid": "IID", "sp": "super_population", "age": "age_BP"}).to_csv(
    f"{prefix}.pca_samples.tsv", sep="\t", index=False)
print("wrote", f"{prefix}.pca_samples.tsv")
PY
