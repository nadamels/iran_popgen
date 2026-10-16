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

# Interactive PCA of the MODERN Iranian samples only (Super-population == Iran and not ancient,
# i.e. "Undergone Imputation" != Yes in the mega table), as ONE self-contained HTML file.
#   - hover a point to see its Population and every other metadata column of the mega table
#     that has a value for that sample (study cluster, source, sex, latitude/longitude, ...)
#   - colour = Population. Markers are ~4x the size used in plot_pca.sh
#   - legend: click a Population to hide/show it; dropdowns switch the PCs on each axis;
#     mouse wheel zooms, drag pans, double-click resets
# Needs only numpy + pandas (no plotly / internet): the HTML opens in any browser.
#
# Usage:   sbatch plot_pca_iran.sh <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [Super-population=Iran]
#          bash   plot_pca_iran.sh <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [Super-population=Iran]
#
#   <prefix>  expects  <prefix>.fam  <prefix>.eigenvec  <prefix>.eigenval
#             (PLINK 1.9 or PLINK 2 eigenvec; with or without a header line)
#   table     defaults to $MASTER_TABLE, else metadata/mega_sample_table.xlsx under $SLURM_SUBMIT_DIR
#             (or the current directory); the .csv next to it is used if openpyxl is missing
#   3rd arg   pca_colors.txt - accepted so the command line matches plot_pca.sh; NOT used here
#             (points are coloured by Population). If the 3rd argument is not a path it is
#             taken as the Super-population name (old behaviour); anything containing "/" or ending
#             in .txt is always treated as the colour file.
#   4th arg   which Super-population to plot (default Iran), or env SUPERPOP=<name>. Ancient samples
#             are always left out (set INCLUDE_ANCIENT=1 to put them back).
#   Set MARKER_SCALE=<x> to resize points further (default 1 = the 4x size described above).
#
# Output: <prefix>.pca_iran.html   (+ <prefix>.pca_iran.tsv with the plotted rows)

set -euo pipefail

PREFIX="${ARGS[0]:-}"
if [[ -z "$PREFIX" ]]; then
    echo "Usage: $0 <prefix> [mega_sample_table.xlsx|.csv] [pca_colors.txt] [Super-population=Iran]" >&2
    exit 1
fi
BASE="${SLURM_SUBMIT_DIR:-$PWD}"
TABLE="${ARGS[1]:-${MASTER_TABLE:-$BASE/metadata/mega_sample_table.xlsx}}"
ARG3="${ARGS[2]:-}"
if [[ -n "$ARG3" && ( "$ARG3" == */* || "$ARG3" == *.txt ) ]]; then
    GROUP="${ARGS[3]:-${SUPERPOP:-Iran}}"           # looks like a path: the (unused) colour file
    [[ -f "$ARG3" ]] || echo "NOTE: $ARG3 not found; ignored (colours file is not used by this script)" >&2
elif [[ -n "$ARG3" && -z "${ARGS[3]:-}" ]]; then
    GROUP="$ARG3"                                   # 3rd argument is a group name, not a colour file
else
    GROUP="${ARGS[3]:-${SUPERPOP:-Iran}}"
fi
MARKER_SCALE="${MARKER_SCALE:-1}"

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

python3 -c "import numpy, pandas" 2>/dev/null || {
    echo "ERROR: python3 in the active environment needs numpy and pandas" >&2; exit 1; }

echo "prefix : $PREFIX"
echo "table  : $TABLE"
echo "group  : Super-population == $GROUP (modern samples only)"

python3 - "$PREFIX" "$TABLE" "$GROUP" "$MARKER_SCALE" "${INCLUDE_ANCIENT:-0}" <<'PY'
import json
import os
import re
import sys
import numpy as np
import pandas as pd

prefix, table_path, group, mscale = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
include_ancient = sys.argv[5] == "1"

# ---- read fam / eigenvec / eigenval ------------------------------------------------------------
fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
use_fid = bool(fam["IID"].duplicated().any() and not fam["FID"].duplicated().any())
if use_fid:
    print("WARNING: IIDs in the fam are not unique but FIDs are; using FID as the sample ID")
fam["sid"] = fam["FID"] if use_fid else fam["IID"]

with open(f"{prefix}.eigenvec") as fh:
    first = fh.readline().split()
if first and first[0].startswith("#"):          # PLINK 2: "#FID IID PC1..." or "#IID PC1..."
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", dtype=str)
    ev.columns = [c.lstrip("#") for c in ev.columns]
    pc_cols = [c for c in ev.columns if c.upper().startswith("PC")]
    if "FID" not in ev:
        ev["FID"] = "0"
else:                                            # PLINK 1.9: FID IID PC1 PC2 ...
    ev = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", header=None, dtype=str)
    pc_cols = [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
    ev.columns = ["FID", "IID"] + pc_cols
ev[pc_cols] = ev[pc_cols].astype(float)
ev = ev.rename(columns={"FID" if use_fid else "IID": "sid"})
if len(pc_cols) < 2:
    sys.exit("ERROR: eigenvec has fewer than 2 PCs")
eigval = np.loadtxt(f"{prefix}.eigenval", ndmin=1)
pct = (100 * eigval / eigval.sum()).tolist()
if ev["sid"].duplicated().any():
    sys.exit("ERROR: duplicate sample IDs in eigenvec; cannot join to the sample table unambiguously")
ev = ev[ev["sid"].isin(set(fam["sid"]))]         # fam defines the sample set

# ---- mega table ---------------------------------------------------------------------------------
if table_path.lower().endswith((".xlsx", ".xls")):
    try:
        t = pd.read_excel(table_path, sheet_name="mega_table", dtype=str)
    except ImportError:
        sys.exit("ERROR: reading .xlsx needs openpyxl; pass the .csv instead")
else:
    t = pd.read_csv(table_path, dtype=str)
need = ["IID", "Population", "Super-population", "Undergone Imputation"]
bad = [c for c in need if c not in t.columns]
if bad:
    sys.exit(f"ERROR: table {table_path} lacks column(s) {bad}; its columns are: {', '.join(t.columns)}")
t = t.drop_duplicates("IID")
sel = t[t["Super-population"] == group]
if sel.empty:
    sys.exit(f"ERROR: no samples with Super-population == {group!r} in the table "
             f"(available: {', '.join(sorted(t['Super-population'].dropna().unique()))})")
is_anc = sel["Undergone Imputation"].str.strip().str.lower() == "yes"
if not include_ancient:
    print(f"leaving out {int(is_anc.sum())} ancient {group} samples (modern only)")
    sel = sel[~is_anc]
    # The table also calls a few reference samples (e.g. SGDP) Super-population Iran. Keep only the
    # study's own Iranians (Source Publication(s) contains IRAN_SOURCE, default "Iranian cohort");
    # IRAN_SOURCE="" keeps every modern Iran-labelled sample.
    iran_src = os.environ.get("IRAN_SOURCE", "Iranian cohort")
    if group.lower() == "iran" and iran_src and "Source Publication(s)" in sel.columns:
        keep = sel["Source Publication(s)"].fillna("").str.contains(iran_src, case=False, regex=False)
        if keep.any() and (~keep).any():
            print(f"leaving out {int((~keep).sum())} other modern Iran-labelled reference samples (Source lacks {iran_src!r}): "
                  + ", ".join(f"{i} [{src}]" for i, src in zip(sel.loc[~keep, 'IID'], sel.loc[~keep, 'Source Publication(s)'])))
        if keep.any():
            sel = sel[keep]
df = ev.merge(sel, left_on="sid", right_on="IID", how="inner")
print(f"{group}: {len(sel)} {'' if include_ancient else 'modern '}samples in the table, {len(df)} present in the PCA")
if df.empty:
    sys.exit("ERROR: none of those samples are in the eigenvec")
df["ancient"] = df["Undergone Imputation"].str.strip().str.lower() == "yes"
df["r"] = 4.0 * mscale                           # ~4x the plot_pca.sh marker radius (px)

def clean(v):
    return None if (v is None or (isinstance(v, float) and np.isnan(v)) or str(v).strip() == "" or str(v) == "nan") else str(v)

# every other table column goes in the tooltip when it has a value, so new columns show up automatically
skip = {"IID", "Population", "Super-population", "Undergone Imputation"}
info_cols = [c for c in t.columns if c not in skip]
empty_cols = [c for c in info_cols if df[c].isna().all()]
if empty_cols:
    print("NOTE: no values for these samples in: " + ", ".join(empty_cols) + " (not shown in the tooltip)")

pops = sorted(df["Population"].fillna("Population not recorded").unique())
rows = []
for _, r in df.iterrows():
    rows.append({
        "id": r["sid"],
        "pop": clean(r["Population"]) or "Population not recorded",
        "pcs": [round(float(r[c]), 6) for c in pc_cols],
        "r": round(float(r["r"]), 2),
        "anc": bool(r["ancient"]),
        "info": {k: clean(r[k]) for k in info_cols if clean(r[k]) is not None},
    })
payload = {"title": f"{prefix.split('/')[-1]}  |  {group}  |  n={len(rows)}",
           "pcs": pc_cols, "pct": [round(x, 2) for x in pct], "pops": pops, "rows": rows,
           "group": group}

HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Iranian PCA</title>
<style>
:root{--bg:#fff;--fg:#1d1d1f;--mut:#6b6b70;--grid:#e6e6ea;--card:#fff;--bd:#cfcfd6}
@media (prefers-color-scheme:dark){:root{--bg:#16161a;--fg:#e8e8ec;--mut:#9a9aa3;--grid:#2c2c33;--card:#202026;--bd:#3a3a44}}
html,body{margin:0;height:100%;background:var(--bg);color:var(--fg);font:13px/1.35 system-ui,sans-serif}
#app{display:flex;height:100%;gap:0}
#main{flex:1;min-width:0;display:flex;flex-direction:column}
#bar{display:flex;gap:14px;align-items:center;padding:10px 16px;flex-wrap:wrap}
#bar h1{font-size:15px;margin:0;font-weight:600}
select,button{background:var(--card);color:var(--fg);border:1px solid var(--bd);border-radius:6px;padding:3px 6px;font:inherit}
#plot{flex:1;min-height:0;position:relative}
svg{width:100%;height:100%;display:block;cursor:grab}
svg:active{cursor:grabbing}
#side{width:290px;border-left:1px solid var(--bd);overflow:auto;padding:10px 12px}
#side h2{font-size:12px;margin:2px 0 8px;color:var(--mut);font-weight:600;text-transform:uppercase;letter-spacing:.04em}
.it{display:flex;align-items:center;gap:7px;padding:2px 0;cursor:pointer;user-select:none}
.it.off{opacity:.35}
.sw{width:11px;height:11px;flex:none;border:1px solid #000}
.nm{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.n{color:var(--mut)}
#tip{position:fixed;pointer-events:none;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px 10px;box-shadow:0 4px 18px rgba(0,0,0,.18);display:none;max-width:340px;z-index:5}
#tip b{font-size:14px}#tip .k{color:var(--mut)}
#note{padding:0 16px 8px;color:var(--mut);font-size:12px}
@media (max-width:720px){#app{flex-direction:column}#side{width:auto;max-height:34%;border-left:0;border-top:1px solid var(--bd)}}
</style></head><body>
<div id="app"><div id="main">
 <div id="bar"><h1 id="ttl"></h1>
  <label>x <select id="sx"></select></label><label>y <select id="sy"></select></label>
  <button id="rst">Reset view</button><button id="all">Show all</button><button id="none">Hide all</button></div>
 <div id="plot"><svg id="svg"></svg></div>
 <div id="note">Hover a point for details &middot; wheel = zoom, drag = pan, double-click = reset &middot; click a Population to hide/show it &middot; colour = Population</div>
</div>
<div id="side"><h2>Population</h2><div id="leg"></div></div></div>
<div id="tip"></div>
<script>
const D = __DATA__;
const NS="http://www.w3.org/2000/svg", svg=document.getElementById("svg"), tip=document.getElementById("tip");
const PAL=["#4e79a7","#f28e2b","#e15759","#76b7b2","#59a14f","#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
 "#1f77b4","#d62728","#2ca02c","#9467bd","#8c564b","#e377c2","#17becf","#bcbd22","#7f7f7f","#ff7f0e"];
const col={}; D.pops.forEach((p,i)=>col[p]=PAL[i%PAL.length]);
const hidden=new Set(); let ix=0, iy=1, view=null, M={l:60,r:20,t:14,b:46};
document.getElementById("ttl").textContent=D.title;
for(const id of ["sx","sy"]){const s=document.getElementById(id);
  D.pcs.forEach((c,i)=>{const o=document.createElement("option");o.value=i;o.textContent=c+" ("+D.pct[i]+"%)";s.appendChild(o)});}
document.getElementById("sy").value=1;
sx.onchange=()=>{ix=+sx.value;resetView()};sy.onchange=()=>{iy=+sy.value;resetView()};
const counts={}; D.rows.forEach(r=>counts[r.pop]=(counts[r.pop]||0)+1);
const leg=document.getElementById("leg");
D.pops.forEach(p=>{const d=document.createElement("div");d.className="it";
  d.innerHTML='<span class="sw" style="background:'+col[p]+'"></span><span class="nm" title="'+p+'">'+p+'</span><span class="n">'+counts[p]+'</span>';
  d.onclick=()=>{hidden.has(p)?hidden.delete(p):hidden.add(p);d.classList.toggle("off");draw()};leg.appendChild(d)});
all.onclick=()=>{hidden.clear();[...leg.children].forEach(c=>c.classList.remove("off"));draw()};
none.onclick=()=>{D.pops.forEach(p=>hidden.add(p));[...leg.children].forEach(c=>c.classList.add("off"));draw()};
function bounds(){let a=[1e9,-1e9,1e9,-1e9];D.rows.forEach(r=>{const x=r.pcs[ix],y=r.pcs[iy];a[0]=Math.min(a[0],x);a[1]=Math.max(a[1],x);a[2]=Math.min(a[2],y);a[3]=Math.max(a[3],y)});
  const px=(a[1]-a[0])*.06||1,py=(a[3]-a[2])*.06||1;return [a[0]-px,a[1]+px,a[2]-py,a[3]+py]}
function resetView(){view=bounds();draw()}
rst.onclick=resetView;svg.ondblclick=resetView;
function dims(){const b=svg.getBoundingClientRect();return {w:b.width,h:b.height}}
function nice(lo,hi,n){const span=hi-lo,st=Math.pow(10,Math.floor(Math.log10(span/n)));const e=[1,2,5,10].find(m=>span/(st*m)<=n)*st;
  const o=[];for(let v=Math.ceil(lo/e)*e;v<=hi;v+=e)o.push(+v.toFixed(10));return o}
function el(n,a,p){const e=document.createElementNS(NS,n);for(const k in a)e.setAttribute(k,a[k]);(p||svg).appendChild(e);return e}
function draw(){
  const {w,h}=dims();svg.setAttribute("viewBox","0 0 "+w+" "+h);svg.innerHTML="";
  const sx_=v=>M.l+(v-view[0])/(view[1]-view[0])*(w-M.l-M.r), sy_=v=>h-M.b-(v-view[2])/(view[3]-view[2])*(h-M.t-M.b);
  const cs=getComputedStyle(document.documentElement),grid=cs.getPropertyValue("--grid"),mut=cs.getPropertyValue("--mut"),fg=cs.getPropertyValue("--fg");
  nice(view[0],view[1],8).forEach(v=>{const x=sx_(v);el("line",{x1:x,x2:x,y1:M.t,y2:h-M.b,stroke:grid});
    el("text",{x:x,y:h-M.b+16,"text-anchor":"middle",fill:mut,"font-size":11}).textContent=v});
  nice(view[2],view[3],6).forEach(v=>{const y=sy_(v);el("line",{x1:M.l,x2:w-M.r,y1:y,y2:y,stroke:grid});
    el("text",{x:M.l-8,y:y+4,"text-anchor":"end",fill:mut,"font-size":11}).textContent=v});
  el("text",{x:(M.l+w-M.r)/2,y:h-8,"text-anchor":"middle",fill:fg,"font-size":12}).textContent=D.pcs[ix]+" ("+D.pct[ix]+"%)";
  const yl=el("text",{x:14,y:(M.t+h-M.b)/2,"text-anchor":"middle",fill:fg,"font-size":12,transform:"rotate(-90 14 "+(M.t+h-M.b)/2+")"});yl.textContent=D.pcs[iy]+" ("+D.pct[iy]+"%)";
  const cp=el("clipPath",{id:"cp"});el("rect",{x:M.l,y:M.t,width:w-M.l-M.r,height:h-M.t-M.b},cp);
  const g=el("g",{"clip-path":"url(#cp)"});
  // modern first (below), ancient on top, big old ancient points below small ones
  const order=D.rows.filter(r=>!hidden.has(r.pop)).sort((a,b)=>(a.anc-b.anc)||(b.r-a.r));
  // modern Iranians = squares with a thick border; ancient samples (only with INCLUDE_ANCIENT=1) = circles, thin outline
  order.forEach(r=>{const x=sx_(r.pcs[ix]),y=sy_(r.pcs[iy]),w0=r.anc?.6:1.6;
    const c=r.anc?el("circle",{cx:x,cy:y,r:r.r,fill:col[r.pop],"fill-opacity":.85,stroke:"#000","stroke-width":w0},g)
                 :el("rect",{x:x-r.r,y:y-r.r,width:2*r.r,height:2*r.r,fill:col[r.pop],"fill-opacity":.85,stroke:"#000","stroke-width":w0},g);
    c.onmouseenter=e=>show(r,e);c.onmousemove=e=>move(e);c.onmouseleave=()=>tip.style.display="none";
    c.onmouseover=()=>c.setAttribute("stroke-width",r.anc?1.6:3);c.onmouseout=()=>c.setAttribute("stroke-width",w0)});
}
function show(r,e){const rows=[["Population","<b>"+r.pop+"</b>"],["Sample ID",r.id]];
  Object.keys(r.info).forEach(k=>{let v=r.info[k];if(/^(lat|lon)/i.test(k)&&!isNaN(+v))v=(+v).toFixed(4);rows.push([k,v])});
  rows.push([D.pcs[ix],r.pcs[ix].toFixed(4)],[D.pcs[iy],r.pcs[iy].toFixed(4)]);
  tip.innerHTML=rows.map(x=>'<div><span class="k">'+x[0]+':</span> '+x[1]+'</div>').join("");tip.style.display="block";move(e)}
function move(e){const w=tip.offsetWidth,h=tip.offsetHeight;let x=e.clientX+14,y=e.clientY+14;
  if(x+w>innerWidth)x=e.clientX-w-14;if(y+h>innerHeight)y=e.clientY-h-14;tip.style.left=x+"px";tip.style.top=y+"px"}
svg.addEventListener("wheel",e=>{e.preventDefault();const b=svg.getBoundingClientRect(),w=b.width,h=b.height;
  const fx=(e.clientX-b.left-M.l)/(w-M.l-M.r),fy=1-(e.clientY-b.top-M.t)/(h-M.t-M.b),k=e.deltaY<0?.85:1/.85;
  const cx=view[0]+fx*(view[1]-view[0]),cy=view[2]+fy*(view[3]-view[2]);
  view=[cx-(cx-view[0])*k,cx+(view[1]-cx)*k,cy-(cy-view[2])*k,cy+(view[3]-cy)*k];draw()},{passive:false});
let drag=null;svg.onmousedown=e=>{drag={x:e.clientX,y:e.clientY,v:view.slice()}};
addEventListener("mouseup",()=>drag=null);
addEventListener("mousemove",e=>{if(!drag)return;const b=svg.getBoundingClientRect();
  const dx=(e.clientX-drag.x)/(b.width-M.l-M.r)*(drag.v[1]-drag.v[0]),dy=(e.clientY-drag.y)/(b.height-M.t-M.b)*(drag.v[3]-drag.v[2]);
  view=[drag.v[0]-dx,drag.v[1]-dx,drag.v[2]+dy,drag.v[3]+dy];draw()});
addEventListener("resize",draw);resetView();
</script></body></html>"""

out = f"{prefix}.pca_iran.html"
with open(out, "w") as fh:
    fh.write(HTML.replace("__DATA__", json.dumps(payload, ensure_ascii=False).replace("</", "<\\/")))
df[["sid", "Population", "Undergone Imputation"] + pc_cols[:4]].to_csv(f"{prefix}.pca_iran.tsv", sep="\t", index=False)
print("wrote", out, f"({len(rows)} points, {len(pops)} populations)")
print("wrote", f"{prefix}.pca_iran.tsv")
PY
