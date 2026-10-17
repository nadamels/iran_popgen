#!/usr/bin/env python3
"""Build the mega sample table.

Rows    : exactly the IIDs of ONE fam (default all_anc_mod.run3.fam in the repo root), in the order of
          metadata/sample_order.txt (new IIDs placed after their fam neighbour; fam order if no such file).
          Each sample is assigned to a group by rule: ancient (listed in a raw study table), Iranian cohort
          (IID numeric or IGP_<n>) or modern global reference.
Columns : IID | Source Publication(s) | Empirical cluster (Study Pop only) | Population |
          Super-population | Sex | Archeological Period | Archeological Date |
          Original Genome Coverage (ancient samples only) | Undergone Imputation | Latitude | Longitude
Metadata: ancient  -> metadata/{laz,nar,matheison}-raw-metadata.xlsx (+ working labels in
                      full_sample_metadata.xlsx)
          modern   -> full_sample_metadata.xlsx (Modern Samples, Orignal Iranian Labels),
                      1kGP_metadata.txt, Hgdp_pops.tsv, SGDP_metadata.xlsx, genome_asia.xlsx

Outputs (metadata/): mega_sample_table.csv / .xlsx, samples_missing_sex.csv,
                     cross_fam_duplicate_iids.csv

Run: python3 scripts/build_mega_table.py [--root <repo root>] [--fam <fam file>]
"""
import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sample_metadata_lib as lib  # noqa: E402

warnings.filterwarnings("ignore")

# ---- policy -----------------------------------------------------------------------------------
# Groups (assigned per sample, see classify()); the fam itself is a single file.
FAM_FILE = "all_anc_mod.run3.fam"
ANC_GROUP = "Ancient"
IRAN_FAM = "Modern Iranians"
REF_GROUP = "Modern global references"
IRAN_ID = re.compile(r"(\d+|IGP_\d+)$")   # Iranian cohort IIDs
COLS = ["IID", "Source Publication(s)", "Empirical cluster (Study Pop only)", "Population", "Super-population",
        "Sex", "Archeological Period", "Archeological Date", "Original Genome Coverage (ancient samples only)",
        "Undergone Imputation", "Latitude", "Longitude"]
CITE = {"L": "Lazaridis et al. 2022 (Science)", "N": "Narasimhan et al. 2019 (Science)",
        "M": "Mathieson et al. 2015 (Nature)"}

# Country -> super-population for ANCIENT samples, using the vocabulary of the modern 'Continent' column.
# This is a judgement call (no source provides it); it is written to the xlsx for review.
SUPERPOP = {
    "Iran": "Iran", "Turkey": "Anatolia", "Armenia": "Caucuses", "Georgia": "Caucuses", "Azerbaijan": "Caucuses",
    "Israel": "Levant", "Jordan": "Levant", "Palestine": "Levant", "Lebanon": "Levant", "Syria": "Levant",
    "Cyprus": "Levant", "Iraq": "Levant", "Egypt": "North-Africa",
    "Spain": "Southern-Europe", "Portugal": "Southern-Europe", "Italy": "Southern-Europe",
    "Greece": "Southern-Europe", "Croatia": "Southern-Europe", "Albania": "Southern-Europe",
    "North Macedonia": "Southern-Europe", "Montenegro": "Southern-Europe", "Serbia": "Southern-Europe",
    "Bosnia and Herzegovina": "Southern-Europe", "Slovenia": "Southern-Europe", "Malta": "Southern-Europe",
    "Kosovo": "Southern-Europe",
    "Bulgaria": "Eastern-Europe", "Romania": "Eastern-Europe", "Moldova": "Eastern-Europe",
    "Ukraine": "Eastern-Europe", "Hungary": "Eastern-Europe", "Czech Republic": "Eastern-Europe",
    "Czechia": "Eastern-Europe", "Poland": "Eastern-Europe", "Slovakia": "Eastern-Europe",
    "Belarus": "Eastern-Europe", "Russia": "Eastern-Europe",
    "Germany": "Western-Europe", "France": "Western-Europe", "Netherlands": "Western-Europe",
    "Belgium": "Western-Europe", "United Kingdom": "Western-Europe", "Ireland": "Western-Europe",
    "Switzerland": "Western-Europe", "Austria": "Western-Europe", "Luxembourg": "Western-Europe",
    "Sweden": "Northern-Europe", "Denmark": "Northern-Europe", "Norway": "Northern-Europe",
    "Finland": "Northern-Europe", "Estonia": "Northern-Europe", "Latvia": "Northern-Europe",
    "Lithuania": "Northern-Europe", "Iceland": "Northern-Europe",
    "Kazakhstan": "Central-Asia", "Turkmenistan": "Central-Asia", "Tajikistan": "Central-Asia",
    "Uzbekistan": "Central-Asia", "Kyrgyzstan": "Central-Asia", "Mongolia": "Central-Asia",
    "Afghanistan": "Central-Asia",
    "Bosnia-Herzegovina": "Southern-Europe",
    "Pakistan": "South-Asia", "India": "South-Asia", "Bangladesh": "South-Asia", "Nepal": "South-Asia",
    "China": "East-Asia", "Japan": "East-Asia", "Brazil": "Americas", "Congo": "Africa", "Ethiopia": "Africa",
}
RUSSIA_ASIAN_LON = 60  # Russian sites east of this longitude (Altai, Yenisei...) -> Central-Asia


def read_fam(path):
    f = pd.read_csv(path, sep=r"\s+", header=None, dtype=str,
                    names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--fam", default=None, help=f"fam file defining the rows (default <root>/{FAM_FILE})")
    args = ap.parse_args()
    root = Path(args.root)
    meta = root / "metadata"
    fam_path = Path(args.fam) if args.fam else root / FAM_FILE

    f = read_fam(fam_path)
    assert f["IID"].is_unique, f"duplicate IIDs inside {fam_path.name}"
    iids = list(f["IID"])
    # Row order: metadata/sample_order.txt lists the established order (ancient studies, Iranian cohort,
    # modern global references last). IIDs missing from it go right after their preceding neighbour in the fam.
    order_file = meta / "sample_order.txt"
    if order_file.exists():
        known = [x.strip() for x in order_file.read_text().splitlines() if x.strip()]
        kset, inset = set(known), set(iids)
        after, prev = {}, None
        for x in iids:
            if x in kset:
                prev = x
            else:
                after.setdefault(prev, []).append(x)
        ordered = after.get(None, [])
        for x in known:
            if x in inset:
                ordered.append(x)
                ordered += after.get(x, [])
        assert sorted(ordered) == sorted(iids)
        iids = ordered
    print(f"fam {fam_path.name}: {len(iids)} samples")

    # ---- reference data ---------------------------------------------------------------------------
    H = lib.load_sources(meta)
    wb = pd.read_excel(meta / "full_sample_metadata.xlsx", sheet_name=None, dtype=str)
    WB_ANC = {
        "European Samples": dict(label="Archaeological.culture", cov="Coverage", country="Country", period="Archaeological Period"),
        "Near Eastern Samples": dict(label="Analysis.Label", cov="Coverage", country="Country"),
        "Central Asian and Iranian Sampl": dict(period="Broad chronological / archaeological period", country="Country",
                                                cov="Coverage"),
        "Southern Arc Samples": dict(label="Analysis_Label", cov="Coverage.on.autosomal.targets", country="Country",
                                     period="Archeological_Period"),
    }
    wbrec = {}
    for sheet, m in WB_ANC.items():
        for _, r in wb[sheet].iterrows():
            rec = wbrec.setdefault(r.iloc[0], {"label": [], "period": [], "cov": []})
            for key in ("label", "period", "cov"):
                if key in m and not pd.isna(lib.clean(r.get(m[key]))):
                    rec[key].append(lib.clean(r.get(m[key])))
    modern = wb["Modern Samples"].set_index("IID")
    ethn = wb["Orignal Iranian Labels"].set_index("IID")["Ethnicity"]

    k1 = pd.read_csv(meta / "1kGP_metadata.txt", sep=None, engine="python", dtype=str)
    hg = pd.read_csv(meta / "Hgdp_pops.tsv", sep="\t", dtype=str)
    sg = pd.read_excel(meta / "SGDP_metadata.xlsx", dtype=str)
    gaf = []
    for sh in ("S1b. Samples", "S1c-SouthaisanSamplesDetails"):
        d = pd.read_excel(meta / "genome_asia.xlsx", sheet_name=sh, dtype=str, header=1)
        d.columns = d.columns.str.strip()
        gaf.append(d[["GA sample ID", "Gender"]])
    ga = pd.concat(gaf).drop_duplicates("GA sample ID")
    SGSEX = {"XY": "M", "XX": "F"}
    sex_src = {  # dataset -> {id: M/F/U}   (priority = order)
        "1000 Genomes Phase 3": dict(zip(k1["sample"], k1["gender"].map({"male": "M", "female": "F"}).fillna("U"))),
        "HGDP": dict(zip(hg["Sample name"], hg["Sex"].map({"M": "M", "F": "F"}).fillna("U"))),
        "SGDP": {**dict(zip(sg["Sample ID (Collaborator)"], sg["Genetic sex assignment"].map(SGSEX).fillna("U"))),
                 **dict(zip(sg["Sample ID (Illumina)"], sg["Genetic sex assignment"].map(SGSEX).fillna("U")))},
        "GenomeAsia 100K": dict(zip(ga["GA sample ID"], ga["Gender"].map({"M": "M", "F": "F"}).fillna("U"))),
    }

    # ---- group per sample: ancient = listed in a raw study table; Iranian cohort = numeric / IGP_n IID ----
    def classify(iid):
        if lib.lookup(H, iid):
            return ANC_GROUP
        return IRAN_FAM if IRAN_ID.match(iid) else REF_GROUP
    rows = [(iid, classify(iid)) for iid in iids]
    membership = {iid: [g] for iid, g in rows}
    dups = {}
    n_fam_rows = len(iids)

    # ---- build rows ---------------------------------------------------------------------------------
    out, missing, conflicts, unmapped = [], [], [], set()
    for iid, group in rows:
        rec = dict.fromkeys(COLS)
        rec["IID"] = iid
        in_anc = group == ANC_GROUP
        why_missing = None
        if in_anc:
            hits = lib.lookup(H, iid)
            w = wbrec.get(iid, {"label": [], "period": [], "cov": []})
            wcov = [lib.num(c) for c in w["cov"]]
            match = [k for k in lib.PRECEDENCE if k in hits and any(
                abs(hits[k]["coverage"] - c) < 2e-3 for c in wcov if not np.isnan(c) and not np.isnan(hits[k]["coverage"]))]
            order = match + [k for k in lib.PRECEDENCE if k in hits and k not in match]
            if not order:
                sys.exit(f"ERROR: ancient IID {iid} ({group}) not found in any raw study table")
            h = hits[order[0]]
            rec["Source Publication(s)"] = "; ".join(CITE[k] for k in order)
            wl = sorted(set(w["label"]))
            pub = lib.clean(h["analysis_label_published"])
            rec["Empirical cluster (Study Pop only)"] = wl[0] if wl else pub
            pop = next((x for x in (pub, h["split_label"], h["archaeological_culture"]) if not pd.isna(x)), np.nan)
            rec["Population"] = pop if not pd.isna(pop) else rec["Empirical cluster (Study Pop only)"]
            country, lon = h["country"], h["longitude"]
            sp = SUPERPOP.get(country)
            if country == "Russia" and not pd.isna(lon) and lon >= RUSSIA_ASIAN_LON:
                sp = "Central-Asia"
            if sp is None:
                unmapped.add(country)
            rec["Super-population"] = sp if sp else "Unassigned"
            rec["Sex"] = h["sex"]
            if rec["Sex"] == "U":
                why_missing = "study reports sex as unknown / not determined"
            per = h["archaeological_period"]
            rec["Archeological Period"] = per if not pd.isna(per) else (w["period"][0] if w["period"] else np.nan)
            d, mu = h["date_description"], h["date_mean_calBP"]
            rec["Archeological Date"] = d if not pd.isna(d) else (f"{mu:.0f} calBP (mean)" if not pd.isna(mu) else np.nan)
            cv = h["coverage"] if not pd.isna(h["coverage"]) else next((c for c in wcov if not np.isnan(c)), np.nan)
            rec["Original Genome Coverage (ancient samples only)"] = round(cv, 6) if not pd.isna(cv) else np.nan
            rec["Undergone Imputation"] = "Yes"
            rec["Latitude"], rec["Longitude"] = h["latitude"], h["longitude"]
        else:
            raw = iid.split(":", 1)[-1]
            found = {n: m[raw] for n, m in sex_src.items() if raw in m}
            rec["Source Publication(s)"] = "; ".join(found) if found else (
                "not recorded in repo metadata (Iranian cohort)" if group == IRAN_FAM else "not recorded in repo metadata")
            in_wb = iid in modern.index
            cont = modern.at[iid, "Continent"] if in_wb else ("Iran" if group == IRAN_FAM else np.nan)
            pop_wb = modern.at[iid, "Population"] if in_wb else np.nan
            if group == IRAN_FAM:
                e = ethn.get(iid, np.nan)
                rec["Population"] = e if isinstance(e, str) else "Iranian (ethnicity not reported)"
                if isinstance(pop_wb, str) and re.fullmatch(r"[0-4]", pop_wb):
                    rec["Empirical cluster (Study Pop only)"] = f"Iran_cluster_{pop_wb}"
            else:
                rec["Population"] = pop_wb
            rec["Super-population"] = cont
            vals = {v for v in found.values() if v in ("M", "F")}
            if len(vals) > 1:
                conflicts.append((iid, found))
            sx = next((v for v in found.values() if v in ("M", "F")), "U")
            rec["Sex"] = sx
            if sx == "U":
                why_missing = ("no sex information in metadata files" if not found
                               else "matched metadata file(s) but sex not reported")
            rec["Undergone Imputation"] = "No"
        out.append(rec)
        if why_missing:
            missing.append({"IID": iid, "Population": rec["Population"], "Super-population": rec["Super-population"],
                            "Source Publication(s)": rec["Source Publication(s)"], "Group": group,
                            "Reason": why_missing})

    df = pd.DataFrame(out)[COLS]
    df["Sex"] = df["Sex"].fillna("U")
    miss = pd.DataFrame(missing)
    dupdf = pd.DataFrame([{"IID": k, "In fam files": "; ".join(v), "Row placed in": v[0]} for k, v in dups.items()])

    # ---- checks -----------------------------------------------------------------------------------
    assert df["IID"].is_unique and set(df["IID"]) == set(membership)
    checks = [("fam rows read", n_fam_rows), ("unique IIDs (= table rows)", len(df)),
              ("IIDs present in more than one fam", len(dups)),
              ("ancient rows (imputed Yes)", int((df["Undergone Imputation"] == "Yes").sum())),
              ("modern rows (imputed No)", int((df["Undergone Imputation"] == "No").sum())),
              ("sex M / F / U", " / ".join(str(int((df.Sex == s).sum())) for s in "MFU")),
              ("sex conflicts between metadata files (first listed used)", len(conflicts)),
              ("ancient rows with Super-population Unassigned", int((df["Super-population"] == "Unassigned").sum())),
              ("ancient rows with latitude+longitude", int((df.Latitude.notna() & df.Longitude.notna()).sum())),
              ("ancient rows missing coverage", int(((df["Undergone Imputation"] == "Yes") &
                                                     df["Original Genome Coverage (ancient samples only)"].isna()).sum()))]
    for c in checks:
        print(f"{c[0]}: {c[1]}")
    if unmapped:
        print("WARNING countries without a super-population mapping:", sorted(map(str, unmapped)))

    # ---- write ------------------------------------------------------------------------------------
    df.to_csv(meta / "mega_sample_table.csv", index=False)
    miss.to_csv(meta / "samples_missing_sex.csv", index=False)
    dupdf.to_csv(meta / "cross_fam_duplicate_iids.csv", index=False)
    gser = pd.Series([r[1] for r in rows])
    groups = pd.DataFrame([(g, fam_path.name, int((gser == g).sum()), "") for g in (ANC_GROUP, IRAN_FAM, REF_GROUP)],
                          columns=["Group", "Fam file", "Rows in table", "Note"])
    dd = pd.DataFrame([
        ("IID", "Sample ID exactly as in the fam."),
        ("Source Publication(s)", "Ancient: every raw study table listing the sample (Lazaridis 2022 / Narasimhan 2019 / Mathieson 2015), the one whose coverage matches the working workbook first. Modern: every reference dataset whose metadata file lists the ID (1000 Genomes, HGDP, SGDP, GenomeAsia 100K); 'not recorded' when none does."),
        ("Empirical cluster (Study Pop only)", "Only for study populations. Ancient: working-workbook analysis label, else the study's published analysis label. Modern Iranians: Iran_cluster_0-4 from the numeric Population column. Blank for global references and for Iranians without a cluster."),
        ("Population", "Ancient: the study's published analysis label (else site/split label, else archaeological culture, else the empirical cluster). Modern Iranians: self-reported ethnicity. Other moderns: population label used in the fam/workbook."),
        ("Super-population", "Modern: 'Continent' of the working workbook. Ancient: derived from the sampling country with the same vocabulary (see ancient_superpop_map sheet; Russian sites east of 60E -> Central-Asia). A judgement, please review."),
        ("Sex", "M/F/U. Ancient: genetic sex from the study. Modern: from 1kGP > HGDP > SGDP > GenomeAsia metadata. U = unknown; those samples are listed in samples_missing_sex.csv."),
        ("Archeological Period", "Ancient only. Study/workbook period code (e.g. EN, BA); blank when none is given."),
        ("Archeological Date", "Ancient only. Date string as given by the study (calBCE/BCE range), else mean calBP."),
        ("Original Genome Coverage (ancient samples only)", "Mean coverage on the 1240k autosomal targets from the chosen study (blank for moderns)."),
        ("Undergone Imputation", "Yes for every ancient sample, No for every modern sample (by rule)."),
        ("Latitude", "Decimal degrees, from the same raw study table as the rest of the ancient row (Lazaridis 2022 'Lat.', Narasimhan 2019 / Mathieson 2015 'Latitude'). Blank when the study gives none and for all modern samples (no coordinates in any modern metadata file)."),
        ("Longitude", "Decimal degrees, same source and rules as Latitude."),
        ("Row rules", "Rows = the IIDs of the single fam named in the 'fam_groups' sheet, in the order of metadata/sample_order.txt. Group (ancient / Iranian cohort / global reference) is assigned per sample from the metadata; see 'fam_groups'."),
    ], columns=["Column", "Definition"])
    smap = pd.DataFrame(sorted(SUPERPOP.items()), columns=["Country", "Super-population"])
    with pd.ExcelWriter(meta / "mega_sample_table.xlsx") as xw:
        df.to_excel(xw, sheet_name="mega_table", index=False)
        dd.to_excel(xw, sheet_name="data_dictionary", index=False)
        groups.to_excel(xw, sheet_name="fam_groups", index=False)
        dupdf.to_excel(xw, sheet_name="cross_fam_duplicates", index=False)
        miss.to_excel(xw, sheet_name="missing_sex", index=False)
        smap.to_excel(xw, sheet_name="ancient_superpop_map", index=False)
        pd.DataFrame(checks, columns=["check", "value"]).to_excel(xw, sheet_name="build_checks", index=False)
        xw.sheets["mega_table"].freeze_panes = "B2"
    print("wrote", meta / "mega_sample_table.csv")


if __name__ == "__main__":
    main()
