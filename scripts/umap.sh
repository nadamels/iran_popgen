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
conda activate /scratch/nme3923/conda-envs/glimpse2

# UMAP of a PLINK binary fileset, restricted to the samples that are in BOTH the fileset and the
# mega sample table. The input fileset can therefore be any subset of the samples: only the IIDs that
# are in <prefix>.fam AND in the table's IID column are used (the rest are reported and ignored).
#   1. intersect fam IIDs with the table IIDs -> keep list
#   2. PCA of those samples (no separate subsetting step)              plink --keep --pca
#   3. UMAP on the top PCs                                             umap-learn
#   4. plot, coloured by Population (colours from Plotting/pca_colors.txt if available)
# The input fileset is assumed to be LD-pruned already: no pruning is done here.
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
#           <out>.interactive.html   interactive UMAP (hover for metadata, colour-by selector, legend toggle/filter,
#                            wheel zoom, drag pan); one self-contained file, no internet needed
#
# Settings (environment variables):
#   NPCS=20           number of PCs fed into UMAP
#   N_NEIGHBORS=15    UMAP n_neighbors
#   MIN_DIST=0.1      UMAP min_dist
#   SEED=42           random seed
#   KEEP_PCA=0        1 = also keep the PLINK PCA (<out>.pca.eigenvec / .eigenval); by default it is
#                     written to a temporary directory and deleted, so it cannot clash with other eigenvecs
#   COLORS_FILE       "name: #hex" lines (default Plotting/pca_colors.txt under $SLURM_SUBMIT_DIR or cwd)
#   COLOR_BY="Population"   table column used for the colours (478 populations in the current table:
#                     colours repeat unless a population is listed in the colours file; the legend shows
#                     only the LEGEND_N=25 largest populations)
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
COLOR_BY="${COLOR_BY:-Population}"
LEGEND_N=${LEGEND_N:-25}
NPCS=${NPCS:-4}
N_NEIGHBORS=${N_NEIGHBORS:-5}
MIN_DIST=${MIN_DIST:-0.1}
SEED=${SEED:-42}
KEEP_PCA=${KEEP_PCA:-0}
THREADS=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-4}}

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
echo "output : ${OUT}.*   PCs: $NPCS   n_neighbors: $N_NEIGHBORS   min_dist: $MIN_DIST"

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

# ---- 2. PCA of the shared samples --------------------------------------------------------------
# The fileset is assumed to be LD-pruned already, so no pruning (and no separate subsetting) is done.
plink --bfile "$IN" --keep "${OUT}.keep.txt" --allow-no-sex --threads "$THREADS" --pca "$NPCS" --out "${TMP}/pca"

if [ "$KEEP_PCA" = "1" ]; then
	cp "${TMP}/pca.eigenvec" "${OUT}.pca.eigenvec"
	cp "${TMP}/pca.eigenval" "${OUT}.pca.eigenval"
fi

# ---- 3./4. UMAP + plot -------------------------------------------------------------------------
python3 -I - "${TMP}/pca.eigenvec" "$TABLE" "${IN}.fam" "$COLORS" "$COLOR_BY" "$NPCS" "$N_NEIGHBORS" "$MIN_DIST" "$SEED" "$OUT" "$LEGEND_N" <<'PY'
import json
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
legend_n = int(sys.argv[11])

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
# Colours: the colours file if it lists the group, else a cycling 60-colour palette (many groups => repeats).
palette = [c for name in ("tab20", "tab20b", "tab20c") for c in plt.get_cmap(name).colors]
fig, ax = plt.subplots(figsize=(11, 8))
for i, g in enumerate(groups):
    sub = df[df[color_by] == g]
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=20, alpha=1, linewidths=0,
               color=col_of.get(g, palette[i % len(palette)]), label=f"{g} (n={len(sub)})")
if len(groups) > legend_n:
    print(f"NOTE: {len(groups)} {color_by} groups; legend shows the {legend_n} largest (all are in the coords.tsv)")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_title(f"{os.path.basename(out)}  |  n={len(df)}, {npcs} PCs, n_neighbors={n_neighbors}, min_dist={min_dist}")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:legend_n], labels[:legend_n], markerscale=3, fontsize=7, frameon=False,
          bbox_to_anchor=(1.01, 1), loc="upper left", title=f"{color_by} (largest {min(legend_n, len(groups))})")
fig.tight_layout()
fig.savefig(f"{out}.png", dpi=200)
fig.savefig(f"{out}.pdf")
print(f"wrote {out}.coords.tsv, {out}.png, {out}.pdf")

# ---- interactive HTML (one self-contained file, no internet needed) ----------------------------------
skip = {"UMAP1", "UMAP2", "Latitude", "Longitude"}
columns = {}
for c in df.columns:
    if c in skip:
        continue
    vals = df[c].fillna("").astype(str)
    uniq = list(dict.fromkeys(vals))
    index = {u: k for k, u in enumerate(uniq)}
    columns[c] = {"u": uniq, "i": [index[v] for v in vals]}
colorable = [c for c, d in columns.items()
             if c != "IID" and 2 <= len(d["u"]) <= 600 and "Date" not in c and "Coverage" not in c]
payload = {
    "title": f"{os.path.basename(out)}  |  n={len(df)}, {npcs} PCs, n_neighbors={n_neighbors}, min_dist={min_dist}",
    "x": [round(float(v), 4) for v in df["UMAP1"]],
    "y": [round(float(v), 4) for v in df["UMAP2"]],
    "cols": columns,
    "colorable": colorable,
    "colorBy": color_by if color_by in colorable else (colorable[0] if colorable else "IID"),
    "fileColors": col_of,
    "colors": {g: matplotlib.colors.to_hex(col_of.get(g, palette[i % len(palette)]))
               for i, g in enumerate(groups)},
}

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Interactive UMAP</title>
<style>
:root{--bg:#fff;--fg:#1d1d1f;--mut:#6b6b70;--grid:#e6e6ea;--card:#fff;--bd:#cfcfd6}
@media (prefers-color-scheme:dark){:root{--bg:#16161a;--fg:#e8e8ec;--mut:#9a9aa3;--grid:#2c2c33;--card:#202026;--bd:#3a3a44}}
html,body{margin:0;height:100%;background:var(--bg);color:var(--fg);font:13px/1.35 system-ui,sans-serif}
#app{display:flex;height:100%}#main{flex:1;min-width:0;display:flex;flex-direction:column}
#bar{display:flex;gap:10px;align-items:center;padding:10px 16px;flex-wrap:wrap}#bar h1{font-size:15px;margin:0;font-weight:600}
select,button,input{background:var(--card);color:var(--fg);border:1px solid var(--bd);border-radius:6px;padding:3px 6px;font:inherit}
#plot{flex:1;min-height:0;position:relative}canvas{position:absolute;inset:0;width:100%;height:100%;cursor:grab}canvas:active{cursor:grabbing}
#side{width:310px;border-left:1px solid var(--bd);display:flex;flex-direction:column;padding:10px 12px;min-height:0}
#side h2{font-size:12px;margin:0 0 6px;color:var(--mut);font-weight:600;text-transform:uppercase;letter-spacing:.04em}
#q{width:100%;box-sizing:border-box;margin-bottom:6px}#leg{flex:1;overflow:auto;min-height:0}
.it{display:flex;align-items:center;gap:7px;padding:2px 0;cursor:pointer;user-select:none}.it.off{opacity:.35}
.sw{width:11px;height:11px;border-radius:50%;flex:none}.nm{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.n{color:var(--mut)}
#tip{position:fixed;pointer-events:none;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px 10px;box-shadow:0 4px 18px rgba(0,0,0,.18);display:none;max-width:380px;z-index:5}
#tip b{font-size:14px}#tip .k{color:var(--mut)}#note{padding:0 16px 8px;color:var(--mut);font-size:12px}
@media (max-width:720px){#app{flex-direction:column}#side{width:auto;max-height:36%;border-left:0;border-top:1px solid var(--bd)}}
</style></head><body>
<div id="app"><div id="main">
 <div id="bar"><h1 id="ttl"></h1><label>colour by <select id="by"></select></label>
  <button id="rst">Reset view</button><button id="all">Show all</button><button id="none">Hide all</button></div>
 <div id="plot"><canvas id="cv"></canvas></div>
 <div id="note">Hover a point for details &middot; wheel = zoom, drag = pan, double-click = reset &middot; click a legend entry to hide/show it</div>
</div>
<div id="side"><h2 id="lh"></h2><input id="q" type="search" placeholder="filter legend..."><div id="leg"></div></div></div>
<div id="tip"></div>
<script>
const D=__DATA__, N=D.x.length, cv=document.getElementById("cv"), ctx=cv.getContext("2d"), tip=document.getElementById("tip");
const PAL=["#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f","#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
 "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"];
const M={l:48,r:16,t:10,b:36}; let by=D.colorBy, groups=[], colorOf=[], hidden=new Set(), view=null, vis=[], W=0, H=0, dpr=1;
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
document.getElementById("ttl").textContent=D.title;
const xs=D.x, ys=D.y; let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;
for(let i=0;i<N;i++){x0=Math.min(x0,xs[i]);x1=Math.max(x1,xs[i]);y0=Math.min(y0,ys[i]);y1=Math.max(y1,ys[i]);}
const padx=(x1-x0)*0.04||1, pady=(y1-y0)*0.04||1, FULL=[x0-padx,x1+padx,y0-pady,y1+pady]; view=FULL.slice();
const sel=document.getElementById("by"); D.colorable.forEach(c=>{const o=document.createElement("option");o.value=o.textContent=c;sel.appendChild(o);}); sel.value=by;

function setColumn(c){
  by=c; hidden=new Set(); const col=D.cols[c], cnt=new Array(col.u.length).fill(0);
  col.i.forEach(k=>cnt[k]++);
  groups=cnt.map((n,k)=>({k,n,name:col.u[k]||"(blank)"})).sort((a,b)=>b.n-a.n);
  colorOf=new Array(col.u.length);
  groups.forEach((g,r)=>{colorOf[g.k]=D.fileColors[col.u[g.k]]||(c===D.colorBy&&D.colors[col.u[g.k]])||PAL[r%PAL.length];});
  document.getElementById("lh").textContent=c+" ("+groups.length+")"; buildLegend(); draw();
}
function buildLegend(){
  const f=document.getElementById("q").value.toLowerCase(), box=document.getElementById("leg"); box.innerHTML="";
  groups.filter(g=>g.name.toLowerCase().includes(f)).forEach(g=>{
    const d=document.createElement("div"); d.className="it"+(hidden.has(g.k)?" off":"");
    d.innerHTML='<span class="sw" style="background:'+colorOf[g.k]+'"></span><span class="nm" title="'+g.name.replace(/"/g,"&quot;")+'"></span><span class="n">'+g.n+'</span>';
    d.querySelector(".nm").textContent=g.name;
    d.onclick=()=>{hidden.has(g.k)?hidden.delete(g.k):hidden.add(g.k); d.classList.toggle("off"); draw();};
    box.appendChild(d);
  });
}
const sx=x=>M.l+(x-view[0])/(view[1]-view[0])*(W-M.l-M.r), sy=y=>H-M.b-(y-view[2])/(view[3]-view[2])*(H-M.t-M.b);
const ux=px=>view[0]+(px-M.l)/(W-M.l-M.r)*(view[1]-view[0]), uy=py=>view[2]+(H-M.b-py)/(H-M.t-M.b)*(view[3]-view[2]);
function ticks(lo,hi){const s=(hi-lo)/6, p=Math.pow(10,Math.floor(Math.log10(s))), m=s/p, st=(m<1.5?1:m<3.5?2:m<7.5?5:10)*p, o=[];
  for(let v=Math.ceil(lo/st)*st; v<=hi; v+=st)o.push(+v.toFixed(10)); return o;}
function resize(){const r=cv.getBoundingClientRect(); dpr=window.devicePixelRatio||1; W=r.width; H=r.height; cv.width=W*dpr; cv.height=H*dpr; ctx.setTransform(dpr,0,0,dpr,0,0); draw();}
let hover=-1;
function draw(){
  ctx.clearRect(0,0,W,H); ctx.font="11px system-ui"; ctx.fillStyle=css("--mut"); ctx.strokeStyle=css("--grid"); ctx.lineWidth=1;
  ticks(view[0],view[1]).forEach(v=>{const p=sx(v); ctx.beginPath(); ctx.moveTo(p,M.t); ctx.lineTo(p,H-M.b); ctx.stroke(); ctx.textAlign="center"; ctx.fillText(v,p,H-M.b+14);});
  ticks(view[2],view[3]).forEach(v=>{const p=sy(v); ctx.beginPath(); ctx.moveTo(M.l,p); ctx.lineTo(W-M.r,p); ctx.stroke(); ctx.textAlign="right"; ctx.fillText(v,M.l-6,p+4);});
  ctx.textAlign="center"; ctx.fillText("UMAP1",(M.l+W-M.r)/2,H-4); ctx.save(); ctx.translate(12,(M.t+H-M.b)/2); ctx.rotate(-Math.PI/2); ctx.fillText("UMAP2",0,0); ctx.restore();
  ctx.save(); ctx.beginPath(); ctx.rect(M.l,M.t,W-M.l-M.r,H-M.t-M.b); ctx.clip();
  const ci=D.cols[by].i; vis=[];
  for(let i=0;i<N;i++){const k=ci[i]; if(hidden.has(k))continue; const px=sx(xs[i]), py=sy(ys[i]);
    if(px<M.l-5||px>W-M.r+5||py<M.t-5||py>H-M.b+5)continue;
    vis.push(i,px,py); ctx.globalAlpha=0.85; ctx.fillStyle=colorOf[k]; ctx.beginPath(); ctx.arc(px,py,3.5,0,6.2832); ctx.fill();}
  ctx.globalAlpha=1;
  if(hover>=0){ctx.strokeStyle=css("--fg"); ctx.lineWidth=2; ctx.beginPath(); ctx.arc(sx(xs[hover]),sy(ys[hover]),6,0,6.2832); ctx.stroke();}
  ctx.restore();
}
function nearest(mx,my){let b=-1,bd=100; for(let j=0;j<vis.length;j+=3){const dx=vis[j+1]-mx,dy=vis[j+2]-my,d=dx*dx+dy*dy; if(d<bd){bd=d;b=vis[j];}} return b;}
cv.addEventListener("mousemove",e=>{
  const r=cv.getBoundingClientRect(), mx=e.clientX-r.left, my=e.clientY-r.top;
  if(drag){view=[drag.v[0]-(mx-drag.x)*drag.sx, drag.v[1]-(mx-drag.x)*drag.sx, drag.v[2]+(my-drag.y)*drag.sy, drag.v[3]+(my-drag.y)*drag.sy]; hover=-1; tip.style.display="none"; draw(); return;}
  const h=nearest(mx,my); if(h!==hover){hover=h; draw();}
  if(h<0){tip.style.display="none"; return;}
  let s=""; Object.keys(D.cols).forEach(c=>{const v=D.cols[c].u[D.cols[c].i[h]]; if(v!=="")s+=(c==="IID"?'<b>'+v+'</b><br>':'<span class="k">'+c+':</span> '+v.replace(/</g,"&lt;")+'<br>');});
  tip.innerHTML=s; tip.style.display="block"; tip.style.left=Math.min(e.clientX+14,innerWidth-tip.offsetWidth-8)+"px"; tip.style.top=Math.min(e.clientY+14,innerHeight-tip.offsetHeight-8)+"px";
});
let drag=null;
cv.addEventListener("mousedown",e=>{const r=cv.getBoundingClientRect(); drag={x:e.clientX-r.left,y:e.clientY-r.top,v:view.slice(),sx:(view[1]-view[0])/(W-M.l-M.r),sy:(view[3]-view[2])/(H-M.t-M.b)};});
addEventListener("mouseup",()=>{drag=null;});
cv.addEventListener("mouseleave",()=>{hover=-1; tip.style.display="none"; draw();});
cv.addEventListener("wheel",e=>{e.preventDefault(); const r=cv.getBoundingClientRect(), mx=e.clientX-r.left, my=e.clientY-r.top, f=e.deltaY<0?0.8:1.25, cx=ux(mx), cy=uy(my);
  view=[cx-(cx-view[0])*f,cx+(view[1]-cx)*f,cy-(cy-view[2])*f,cy+(view[3]-cy)*f]; draw();},{passive:false});
const reset=()=>{view=FULL.slice(); draw();};
cv.addEventListener("dblclick",reset); document.getElementById("rst").onclick=reset;
document.getElementById("all").onclick=()=>{hidden=new Set(); buildLegend(); draw();};
document.getElementById("none").onclick=()=>{hidden=new Set(groups.map(g=>g.k)); buildLegend(); draw();};
document.getElementById("q").oninput=buildLegend; sel.onchange=()=>setColumn(sel.value);
addEventListener("resize",resize); setColumn(by); resize();
</script></body></html>
"""
with open(f"{out}.interactive.html", "w") as fh:
    fh.write(HTML.replace("__DATA__", json.dumps(payload).replace("</", "<\\/")))
print(f"wrote {out}.interactive.html")
PY
