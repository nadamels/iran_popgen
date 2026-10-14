"""Shared helpers for turning the three raw ancient-DNA metadata tables (Lazaridis 2022,
Narasimhan 2019, Mathieson 2015) into one common schema. Lifted from build_sample_table.py so the
mega-table builder uses exactly the same harmonisation rules."""
import numpy as np
import pandas as pd

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


def load_sources(meta_dir):
    """Read the three raw tables and return {study_key: harmonised DataFrame indexed by instance_id}."""
    import warnings
    warnings.filterwarnings("ignore")
    laz = pd.read_excel(meta_dir / "laz-raw-metadata.xlsx", dtype=str)
    nar = pd.read_excel(meta_dir / "nar-raw-metadata.xlsx", sheet_name=0, header=2, dtype=str)
    mat = pd.read_excel(meta_dir / "matheison-raw-metadata.xlsx", dtype=str)
    H = {"L": from_laz(laz), "N": from_nar(nar), "M": from_mat(mat)}
    for k, df in H.items():
        for col in COMMON:
            if col not in df:
                df[col] = np.nan
        H[k] = df[COMMON].dropna(subset=["instance_id"]).drop_duplicates("instance_id").set_index(
            "instance_id", drop=False)
    return H


def lookup(H, sid):
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
