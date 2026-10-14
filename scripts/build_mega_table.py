#!/usr/bin/env python3
"""Build the mega sample table.

Rows    : exactly the IIDs in cleaned_fams/*.fam (nothing else), one row per unique IID, grouped by
          fam file in FAM_ORDER (modern global references last), original fam row order kept inside
          each group.
Columns : IID | Source Publication(s) | Empirical cluster (Study Pop only) | Population |
          Super-population | Sex | Archeological Period | Archeological Date |
          Original Genome Coverage (ancient samples only) | Undergone Imputation
Metadata: ancient  -> metadata/{laz,nar,matheison}-raw-metadata.xlsx (+ working labels in
                      full_sample_metadata.xlsx)
          modern   -> full_sample_metadata.xlsx (Modern Samples, Orignal Iranian Labels),
                      1kGP_metadata.txt, Hgdp_pops.tsv, SGDP_metadata.xlsx, genome_asia.xlsx

Outputs (metadata/): mega_sample_table.csv / .xlsx, samples_missing_sex.csv,
                     cross_fam_duplicate_iids.csv

Run: python3 scripts/build_mega_table.py [--root <repo root>]
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
# Output order of the datasets. Samples from one fam stay together; global references go last.
FAM_ORDER = [
    ("Mathieson2015 / European ancient", "euro_samples_merged.filtered.renamediids.dupSampVarem.flipped.amb.rel2ndrem.fam"),
    ("Narasimhan2019 + 54 ancient", "Narasimhan2019_plus54Full.fam"),
    ("Lazaridis2022 / Southern Arc ancient", "southern_arc_samples_merged.filtered.renamediids.dupSampVarem.flipped.amb.rel2ndrem.fam"),
    ("Modern Iranians", "all_avr0.2_gq20.SNPsOnly.vcf.filtered.dupSampVarem.flipped.amb.rel1strem.fam"),
    ("Modern global references", "GA_Paps_refs.related_removed.filtered.flipped.amb.filtered.dupSampVarem.rel2ndrem.fam"),
]
IRAN_FAM = "Modern Iranians"
COLS = ["IID", "Source Publication(s)", "Empirical cluster (Study Pop only)", "Population", "Super-population",
        "Sex", "Archeological Period", "Archeological Date", "Original Genome Coverage (ancient samples only)",
        "Undergone Imputation"]
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
    root = Path(ap.parse_args().root)
    meta, famdir = root / "metadata", root / "cleaned_fams"

    # ---- rows: IIDs from the cleaned fams -------------------------------------------------------
    rows, membership = [], {}
    for group, fn in FAM_ORDER:
        f = read_fam(famdir / fn)
        assert f["IID"].is_unique, f"duplicate IIDs inside {fn}"
        for iid in f["IID"]:
            membership.setdefault(iid, []).append(group)
            if len(membership[iid]) == 1:
                rows.append((iid, group))
    n_fam_rows = sum(len(v) for v in membership.values())
    dups = {k: v for k, v in membership.items() if len(v) > 1}
    print(f"fam rows {n_fam_rows}; unique IIDs {len(membership)}; IIDs in >1 fam {len(dups)}")

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

    # ---- build rows ---------------------------------------------------------------------------------
    out, missing, conflicts, unmapped = [], [], [], set()
    for iid, group in rows:
        rec = dict.fromkeys(COLS)
        rec["IID"] = iid
        in_anc = group not in ("Modern Iranians", "Modern global references")
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
    groups = pd.DataFrame([(g, fn, int((pd.Series([r[1] for r in rows]) == g).sum()), "") for g, fn in FAM_ORDER],
                          columns=["Group (output order)", "Fam file", "Rows in table", "Note"])
    groups["Fam rows in file"] = [len(read_fam(famdir / fn)) for _, fn in FAM_ORDER]
    dd = pd.DataFrame([
        ("IID", "Sample ID exactly as in the cleaned fam."),
        ("Source Publication(s)", "Ancient: every raw study table listing the sample (Lazaridis 2022 / Narasimhan 2019 / Mathieson 2015), the one whose coverage matches the working workbook first. Modern: every reference dataset whose metadata file lists the ID (1000 Genomes, HGDP, SGDP, GenomeAsia 100K); 'not recorded' when none does."),
        ("Empirical cluster (Study Pop only)", "Only for study populations. Ancient: working-workbook analysis label, else the study's published analysis label. Modern Iranians: Iran_cluster_0-4 from the numeric Population column. Blank for global references and for Iranians without a cluster."),
        ("Population", "Ancient: the study's published analysis label (else site/split label, else archaeological culture, else the empirical cluster). Modern Iranians: self-reported ethnicity. Other moderns: population label used in the fam/workbook."),
        ("Super-population", "Modern: 'Continent' of the working workbook. Ancient: derived from the sampling country with the same vocabulary (see ancient_superpop_map sheet; Russian sites east of 60E -> Central-Asia). A judgement, please review."),
        ("Sex", "M/F/U. Ancient: genetic sex from the study. Modern: from 1kGP > HGDP > SGDP > GenomeAsia metadata. U = unknown; those samples are listed in samples_missing_sex.csv."),
        ("Archeological Period", "Ancient only. Study/workbook period code (e.g. EN, BA); blank when none is given."),
        ("Archeological Date", "Ancient only. Date string as given by the study (calBCE/BCE range), else mean calBP."),
        ("Original Genome Coverage (ancient samples only)", "Mean coverage on the 1240k autosomal targets from the chosen study (blank for moderns)."),
        ("Undergone Imputation", "Yes for every ancient sample, No for every modern sample (by rule)."),
        ("Row rules", "Rows = unique IIDs of cleaned_fams/*.fam. Grouped by fam in the order of the 'fam_groups' sheet, original fam order kept inside each group. A sample in several fams appears once, in the first group; see cross_fam_duplicate_iids.csv."),
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
