#!/usr/bin/env python3
"""Build the single definitive sample table for the Iran popgen project.

Spine      : fam_files/ancient.modern.iran.fam  (every row = one genotype record)
             + ancient samples that appear in the working workbook but NOT in the
               fam file (listed as Excluded so nothing silently disappears).
Group defs : metadata/full_sample_metadata.xlsx (working workbook)
Harmonised : metadata/laz-raw-metadata.xlsx   (Lazaridis et al. 2022, Science abm4247)
             metadata/nar-raw-metadata.xlsx   (Narasimhan et al. 2019, Science aat7487)
             metadata/matheison-raw-metadata.xlsx (Mathieson et al. 2015, Nature 16152)

Outputs    : metadata/master_sample_table.tsv
             metadata/master_sample_table.xlsx (table + data dictionary + checks)

Run: python3 scripts/build_sample_table.py
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"

# ---- tunable policy ---------------------------------------------------------
LOW_COVERAGE_FLAG = 0.1   # x on 1240k targets; flagged, NOT excluded
STUDY = {"L": "Lazaridis2022", "N": "Narasimhan2019", "M": "Mathieson2015"}
PRECEDENCE = ["L", "N", "M"]  # tie-break when coverage cannot pick a source
MISSING = {"", "..", "?", "-", "n/a", "na", "nan", "none", "n/a (female)"}


def clean(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    s = str(v).strip()
    return np.nan if s.lower() in MISSING else s


def num(v):
    v = clean(v)
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def norm_sex(v):
    v = clean(v)
    if pd.isna(v):
        return "U"
    v = v.upper()
    return v if v in ("M", "F") else "U"


# ---- load -------------------------------------------------------------------
fam = pd.read_csv(META.parent / "fam_files/ancient.modern.iran.fam", sep=r"\s+",
                  header=None, dtype=str, names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
fam["fam_row"] = np.arange(1, len(fam) + 1)
wb = pd.read_excel(META / "full_sample_metadata.xlsx", sheet_name=None, dtype=str)
laz = pd.read_excel(META / "laz-raw-metadata.xlsx", dtype=str)
nar = pd.read_excel(META / "nar-raw-metadata.xlsx", sheet_name=0, header=2, dtype=str)
mat = pd.read_excel(META / "matheison-raw-metadata.xlsx", dtype=str)

# ---- harmonise each raw study into one common schema -------------------------
COMMON = ["instance_id", "master_id", "analysis_label_published", "split_label", "archaeological_culture",
          "archaeological_period", "location", "country", "broad_region", "latitude", "longitude", "sex",
          "coverage", "snps_hit", "date_mean_calBP", "date_description", "skeletal_element", "skeletal_code",
          "age_at_death", "data_type", "mtdna_haplogroup", "y_haplogroup", "library_ids", "original_publication",
          "source_qc_status", "source_qc_detail"]


def qc_from_text(t):
    t = clean(t)
    if pd.isna(t):
        return np.nan
    u = t.upper()
    if u.startswith("FAIL") or "QUESTIONABLE_CRITICAL" in u:
        return "FAIL"
    if "QUESTIONABLE" in u:
        return "QUESTIONABLE"
    return "PASS"


def from_laz(d):
    c = list(d.columns)
    o = pd.DataFrame(index=d.index)
    o["instance_id"] = d[c[0]].map(clean)
    o["master_id"] = d["Master ID"].map(clean)
    o["analysis_label_published"] = d["Analysis_Label"].map(clean)
    o["location"] = d["Locality"].map(clean)
    o["country"] = d["Country"].map(clean)
    o["latitude"] = d["Lat."].map(num)
    o["longitude"] = d["Long."].map(num)
    o["sex"] = d["Sex"].map(norm_sex)
    o["coverage"] = d["Coverage on autosomal targets"].map(num)
    o["snps_hit"] = d["SNPs hit on autosomal targets"].map(num)
    o["date_mean_calBP"] = d[c[8]].map(num)
    o["date_description"] = d[c[10]].map(clean)
    o["skeletal_element"] = d["Skeletal element"].map(clean)
    o["skeletal_code"] = d["Skeletal code"].map(clean)
    o["age_at_death"] = d["Age at Death"].map(clean)
    o["data_type"] = d["Data type"].map(clean)
    o["mtdna_haplogroup"] = d[c[24]].map(clean)
    o["y_haplogroup"] = d[c[23]].map(clean)
    o["library_ids"] = d["LibraryID(s)"].map(clean)
    o["original_publication"] = d["Publication"].map(clean)
    o["source_qc_detail"] = d[c[-1]].map(clean)
    o["source_qc_status"] = d[c[-1]].map(qc_from_text)
    return o


def from_nar(d):
    c = list(d.columns)
    o = pd.DataFrame(index=d.index)
    o["instance_id"] = d[c[0]].map(clean)
    o["master_id"] = d["Master ID"].map(clean)
    o["analysis_label_published"] = d[c[11]].map(clean)
    o["split_label"] = d[c[12]].map(clean)
    o["archaeological_period"] = d[c[10]].map(clean)
    o["location"] = d["Location"].map(clean)
    o["country"] = d["Country"].map(clean)
    o["broad_region"] = d[c[16]].map(clean)
    o["latitude"] = d["Latitude"].map(num)
    o["longitude"] = d["Longitude"].map(num)
    o["sex"] = d["Sex"].map(norm_sex)
    o["coverage"] = d[c[23]].map(num)
    o["snps_hit"] = d[c[24]].map(num)
    o["date_mean_calBP"] = d[c[8]].map(num)
    o["date_description"] = d[c[9]].map(clean)
    o["skeletal_element"] = d[c[3]].map(clean)
    o["skeletal_code"] = d["Skeletal code"].map(clean)
    o["data_type"] = d["Data type"].map(clean)
    o["mtdna_haplogroup"] = d["mtDNA haplogroup"].map(clean)
    o["y_haplogroup"] = d[c[21]].map(clean)
    o["library_ids"] = d[c[27]].map(clean)
    o["original_publication"] = d[c[6]].map(clean)
    passed = d[c[13]].map(clean)
    o["source_qc_detail"] = d[c[-1]].map(clean)
    st = d[c[-1]].map(qc_from_text)
    st[passed == "No"] = "FAIL"   # study's own "Passed all analysis filters" column wins
    o["source_qc_status"] = st
    return o


def from_mat(d):
    o = pd.DataFrame(index=d.index)
    o["instance_id"] = d["Unique ID"].map(clean)
    o["master_id"] = o["instance_id"]
    o["archaeological_culture"] = d["Archaeological culture"].map(clean)
    # Mathieson's "Selection label" columns are broad region/period tags (e.g. HG, SA), not analysis
    # clusters, so analysis_label_published stays empty for this study.
    o["split_label"] = np.nan
    o["location"] = d["Location"].map(clean)
    o["country"] = d["Country"].map(clean)
    o["latitude"] = d["Latitude"].map(num)
    o["longitude"] = d["Longitude"].map(num)
    o["sex"] = d["Sex (genetically determined)"].map(norm_sex)
    o["coverage"] = d["Coverage"].map(num)
    o["snps_hit"] = d["SNPs"].map(num)
    o["date_mean_calBP"] = (d["Min Date"].map(num) + d["Max Date"].map(num)) / 2
    o["date_description"] = d["Date (2-sigma)"].map(clean)
    o["skeletal_element"] = d["Skeletal element"].map(clean)
    o["mtdna_haplogroup"] = d["mtDNA haplogroup"].map(clean)
    o["y_haplogroup"] = d["Y haplogroup"].map(clean)
    o["library_ids"] = d["Libraries"].map(clean)
    o["original_publication"] = d["Reference for first report of data from this sample"].map(clean)
    o["data_type"] = d["Data source"].map(clean)
    bad = (d["Selection"].map(clean) == "N") | (d["Popgen"].map(clean) == "N")
    o["source_qc_status"] = np.where(bad, "FAIL", "PASS")
    o["source_qc_detail"] = ("Popgen=" + d["Popgen"].fillna("?") + "; Selection=" + d["Selection"].fillna("?")
                             + "; Warnings=" + d["Warnings"].fillna("-"))
    return o


H = {"L": from_laz(laz), "N": from_nar(nar), "M": from_mat(mat)}
for k, df in H.items():
    for col in COMMON:
        if col not in df:
            df[col] = np.nan
    H[k] = df[COMMON].dropna(subset=["instance_id"]).drop_duplicates("instance_id").set_index("instance_id", drop=False)


def lookup(sid):
    """Return {study_key: row} for every raw study listing this sample (instance ID, else master ID)."""
    hits = {}
    for k, df in H.items():
        if sid in df.index:
            hits[k] = df.loc[sid]
        else:
            m = df[df["master_id"] == sid]
            if len(m):
                hits[k] = m.iloc[0]
    return hits


# ---- working-workbook ancient labels ----------------------------------------
WB_ANCIENT = {
    "European Samples": dict(label="Archaeological.culture", cov="Coverage", country="Country", period="Archaeological Period"),
    "Near Eastern Samples": dict(label="Analysis.Label", cov="Coverage", country="Country"),
    "Central Asian and Iranian Sampl": dict(period="Broad chronological / archaeological period", country="Country",
                                            region="Broad Geographic Region (Genetic clustering may sometimes differ)", cov="Coverage"),
    "Southern Arc Samples": dict(label="Analysis_Label", cov="Coverage.on.autosomal targets".replace(" ", "."), country="Country",
                                 period="Archeological_Period"),
}
wbrec = {}
for sheet, m in WB_ANCIENT.items():
    d = wb[sheet]
    for _, r in d.iterrows():
        sid = r.iloc[0]
        rec = wbrec.setdefault(sid, {"sheets": [], "label": [], "period": [], "region": [], "country": [], "cov": []})
        rec["sheets"].append(sheet)
        for key in ("label", "period", "region", "country", "cov"):
            if key in m and clean(r.get(m[key])) is not np.nan:
                rec[key].append(clean(r.get(m[key])))

# ---- modern ------------------------------------------------------------------
mod = wb["Modern Samples"].merge(wb["Sheet1"][["IID", "Ethnicity"]], on="IID", how="left")
modrec = mod.set_index("IID")

# ---- assemble ----------------------------------------------------------------
hc_list = set(wb["High Coverage Samples"].iloc[:, 0])
fam_ids = set(fam["IID"])
dup_iids = set(fam.loc[fam["IID"].duplicated(keep=False), "IID"])

rows = []
spine = [(r.IID, r.FID, r.fam_row) for r in fam.itertuples()]
spine += [(sid, np.nan, np.nan) for sid in wbrec if sid not in fam_ids]  # workbook-only ancients


def first(*vals):
    for v in vals:
        if not (v is None or (isinstance(v, float) and np.isnan(v))):
            return v
    return np.nan


for sid, fid, frow in spine:
    in_fam = not pd.isna(frow)
    r = dict(sample_id=sid, fam_FID=fid, fam_row=frow, in_fam="Y" if in_fam else "N")
    flags = []
    notes = []
    if sid in modrec.index and sid not in wbrec:
        m = modrec.loc[sid]
        pop = m["Population"]
        r["sample_type"] = "modern"
        r["study"] = "not recorded in repository metadata"
        r["all_matching_studies"] = np.nan
        r["original_publication"] = np.nan
        if re.fullmatch(r"[0-4]", str(pop)):
            r["empirical_cluster"] = f"Iran_cluster_{pop}"
            r["empirical_cluster_basis"] = "numeric cluster in 'Modern Samples'/Sheet1 (clustering method not documented in repo)"
        elif m["Continent"] == "Iran":
            r["empirical_cluster"] = "Iran_unclustered"
            r["empirical_cluster_basis"] = "Iranian-population individual without a numeric cluster"
            flags.append("no_cluster_assigned")
        else:
            r["empirical_cluster"] = pop
            r["empirical_cluster_basis"] = "reference population label (no empirical clustering applied)"
        r["self_reported_label"] = first(clean(m["Ethnicity"]), pop if not re.fullmatch(r"[0-4]", str(pop)) else np.nan)
        r["self_reported_label_type"] = ("ethnicity (Sheet1)" if not pd.isna(clean(m["Ethnicity"]))
                                         else "population label" if not re.fullmatch(r"[0-4]", str(pop)) else "not reported")
        r["sampling_location"] = np.nan
        r["country"] = "Iran" if m["Continent"] == "Iran" else np.nan
        r["broad_region"] = m["Continent"]
        r["sex"] = "U"
        r["coverage"] = np.nan
        r["source_qc_status"] = "not available"
    else:
        hits = lookup(sid)
        w = wbrec.get(sid, {})
        r["sample_type"] = "ancient"
        # choose primary study: the one whose coverage equals the working-workbook coverage
        wcov = [num(c) for c in w.get("cov", [])]
        match = [k for k in PRECEDENCE if k in hits and any(abs(hits[k]["coverage"] - c) < 2e-3 for c in wcov
                                                           if not np.isnan(c) and not np.isnan(hits[k]["coverage"]))]
        pk = (match or [k for k in PRECEDENCE if k in hits] or [None])[0]
        h = hits.get(pk)
        r["study"] = STUDY.get(pk, "UNMATCHED")
        r["all_matching_studies"] = ";".join(STUDY[k] for k in PRECEDENCE if k in hits) or np.nan
        if len(hits) > 1:
            notes.append("in_multiple_studies")
        if pk is None:
            flags.append("no_raw_metadata")
            h = pd.Series({c: np.nan for c in COMMON})
        wl = sorted(set(w.get("label", [])))
        if len(wl) > 1:
            flags.append("conflicting_workbook_labels")
        label = first(wl[0] if wl else np.nan, h["analysis_label_published"])
        r["empirical_cluster"] = label
        r["empirical_cluster_basis"] = ("working-workbook analysis label" if wl else
                                        f"analysis label published by {STUDY.get(pk)}" if not pd.isna(label) else "none")
        pub_label = clean(h["analysis_label_published"])
        if wl and not pd.isna(pub_label) and wl[0] != pub_label:
            notes.append("workbook_label_differs_from_published")
        r["analysis_label_published"] = pub_label
        r["self_reported_label"] = first(h["archaeological_culture"], h["split_label"], h["archaeological_period"],
                                         (w.get("period") or [np.nan])[0])
        r["self_reported_label_type"] = ("archaeological culture" if not pd.isna(h["archaeological_culture"]) else
                                         "site/split label" if not pd.isna(h["split_label"]) else
                                         "archaeological period" if not pd.isna(r["self_reported_label"]) else "not reported")
        r["sampling_location"] = h["location"]
        r["country"] = first(h["country"], (w.get("country") or [np.nan])[0])
        r["broad_region"] = first(h["broad_region"], (w.get("region") or [np.nan])[0])
        r["sex"] = h["sex"] if pk else "U"
        r["coverage"] = first(h["coverage"], next((c for c in wcov if not np.isnan(c)), np.nan))
        for c in ["master_id", "archaeological_culture", "archaeological_period", "latitude", "longitude", "snps_hit",
                  "date_mean_calBP", "date_description", "skeletal_element", "skeletal_code", "age_at_death",
                  "data_type", "mtdna_haplogroup", "y_haplogroup", "library_ids", "original_publication",
                  "source_qc_status", "source_qc_detail"]:
            r[c] = h[c]
        if pd.isna(r["archaeological_period"]):
            r["archaeological_period"] = (w.get("period") or [np.nan])[0]
        r["workbook_sheets"] = ";".join(w.get("sheets", [])) or np.nan
        if sid.endswith("_d"):
            notes.append("damage_restricted_instance")

    cov = r.get("coverage", np.nan)
    if r["sample_type"] == "ancient" and not pd.isna(cov) and cov < LOW_COVERAGE_FLAG:
        flags.append(f"low_coverage_<{LOW_COVERAGE_FLAG}")
    if sid in dup_iids:
        flags.append("duplicate_IID_in_fam(FID 0 and 1)")
    r["in_high_coverage_list"] = "Y" if sid in hc_list else "N"

    # inclusion policy
    reasons = []
    if not in_fam:
        reasons.append("not in fam / genotype dataset")
    if r.get("source_qc_status") == "FAIL":
        reasons.append("failed source-study QC")
    r["inclusion_status"] = "Excluded" if reasons else ("Included_flagged" if flags else "Included")
    r["exclusion_reason"] = "; ".join(reasons) if reasons else np.nan
    r["flags"] = ";".join(flags) if flags else np.nan
    r["notes"] = ";".join(notes) if notes else np.nan
    rows.append(r)

out = pd.DataFrame(rows)
order = ["sample_id", "fam_FID", "fam_row", "sample_type", "study", "all_matching_studies", "original_publication",
         "empirical_cluster", "empirical_cluster_basis", "self_reported_label", "self_reported_label_type",
         "analysis_label_published", "sampling_location", "country", "broad_region", "latitude", "longitude", "sex",
         "coverage", "snps_hit", "archaeological_period", "date_mean_calBP", "date_description", "skeletal_element",
         "skeletal_code", "age_at_death", "data_type", "mtdna_haplogroup", "y_haplogroup", "library_ids",
         "master_id", "source_qc_status", "source_qc_detail", "in_fam", "in_high_coverage_list", "workbook_sheets",
         "inclusion_status", "exclusion_reason", "flags", "notes"]
out = out[[c for c in order if c in out]]
out.to_csv(META / "master_sample_table.tsv", sep="\t", index=False)

# ---- checks + dictionary -----------------------------------------------------
checks = []
def chk(name, val):
    checks.append((name, val)); print(f"{name}: {val}")
chk("fam rows", len(fam))
chk("table rows", len(out))
chk("fam rows missing from table", int((~fam["fam_row"].isin(out["fam_row"].dropna())).sum()))
chk("duplicate (FID,IID) keys among fam rows", int(out[out.in_fam == "Y"].duplicated(["fam_FID", "sample_id"]).sum()))
chk("duplicate sample_id among non-fam rows", int(out[out.in_fam == "N"].duplicated("sample_id").sum()))
chk("ancient in fam without raw metadata", int(((out.sample_type == "ancient") & (out.in_fam == "Y") & (out.study == "UNMATCHED")).sum()))
for k, v in out.groupby(["sample_type", "inclusion_status"]).size().items():
    chk(f"count {k}", int(v))
for k, v in out[out.sample_type == "ancient"].study.value_counts().items():
    chk(f"ancient study {k}", int(v))

dd = pd.DataFrame([
    ("sample_id", "IID in the fam file (instance ID; '_d' = damage-restricted). Key together with fam_FID."),
    ("fam_FID / fam_row", "FID and 1-based line number in ancient.modern.iran.fam. Blank for ancient samples not in the fam file."),
    ("sample_type", "ancient or modern."),
    ("study", "Study whose metadata record is used: Lazaridis2022 / Narasimhan2019 / Mathieson2015. Where a sample is in several raw tables, the one whose coverage equals the working-workbook coverage is chosen; ties go Laz > Nar > Mathieson. Moderns: not recorded in repo."),
    ("all_matching_studies", "Every raw table that lists the sample."),
    ("original_publication", "Study in which the DNA was first published, as stated by the raw table."),
    ("empirical_cluster", "The group label used in analyses. Ancient: workbook analysis label (Analysis_Label / Analysis.Label / Archaeological.culture), falling back to the published analysis label. Modern Iranians: Iran_cluster_0-4 from the numeric 'Population' column. Other moderns: the reference population label."),
    ("empirical_cluster_basis", "How empirical_cluster was derived."),
    ("self_reported_label", "Modern: self-reported ethnicity (Sheet1) or population label. Ancient: archaeological culture / site label / period as reported by the study (there is no self-report for ancients)."),
    ("self_reported_label_type", "What kind of label self_reported_label is."),
    ("analysis_label_published", "Analysis label exactly as in the raw study table, kept to expose drift against the workbook."),
    ("sampling_location", "Site/locality (ancient). Not available for moderns."),
    ("country / broad_region", "Country; broad region (modern: Continent column; ancient: Nar broad region or workbook region where given)."),
    ("latitude / longitude", "Decimal degrees from the raw study."),
    ("sex", "Genetic sex M/F; U = unknown. Moderns: U (fam sex is -9 for all rows)."),
    ("coverage / snps_hit", "Mean coverage and SNPs hit on the autosomal 1240k targets from the chosen study. Moderns: blank."),
    ("date_mean_calBP", "Mean calibrated date BP (1950). Mathieson = mean of Min/Max date."),
    ("source_qc_status", "PASS/QUESTIONABLE/FAIL from the study's own assessment (Laz ASSESSMENT, Nar 'Passed all analysis filters' + assessment, Mathieson Popgen/Selection)."),
    ("in_high_coverage_list", "Sample is on the workbook's 'High Coverage Samples' sheet. NB that sheet lists all ancient fam samples, including coverage <0.1x, so it is not a coverage filter."),
    ("inclusion_status", f"Included | Included_flagged (in fam, but see flags) | Excluded. Excluded = not in the fam file, or failed the source study's QC. Low-coverage threshold for flags: {LOW_COVERAGE_FLAG}x."),
    ("exclusion_reason / flags / notes", "Why excluded; flags = issues that make a sample Included_flagged (low coverage, conflicting labels, no cluster, duplicate IID, no raw metadata); notes = informational only (in several studies, workbook label differs from the published label, damage-restricted instance)."),
], columns=["column", "definition"])

with pd.ExcelWriter(META / "master_sample_table.xlsx") as xw:
    out.to_excel(xw, sheet_name="sample_table", index=False)
    dd.to_excel(xw, sheet_name="data_dictionary", index=False)
    pd.DataFrame(checks, columns=["check", "value"]).to_excel(xw, sheet_name="build_checks", index=False)
    out.groupby(["sample_type", "study", "inclusion_status"]).size().rename("n").reset_index().to_excel(xw, sheet_name="summary", index=False)
print("wrote", META / "master_sample_table.tsv")
