#!/usr/bin/env python3
"""Define sub-populations of a sample set from PCA + ADMIXTURE results (step 3 of define_subpops.sh).

Three independent views of the structure are computed and then reconciled:
  1. GMM   Gaussian mixture on the top PCs; the number of components is chosen by BIC (k = 1..KMAX_GMM).
  2. HDBSCAN  density-based clustering on the same PCs (no k to choose; sparse samples -> noise).
  3. ADMIXTURE  hard assignment at the best K (lowest cross-validation error): the component with the
     largest ancestry fraction, if it is >= ADMIX_MIN, else "admixed".
The final sub-populations are the GMM clusters (they have a model-based posterior), after
  * clusters smaller than MIN_SIZE are set to "unassigned",
  * samples whose posterior is < MIN_POST are flagged ("low_confidence"),
and every sample gets the label each method gave it, so disagreement is visible.
Agreement between methods is summarised with the adjusted Rand index.

Usage: subpop_consensus.py <eigenvec> <eigenval|-> <cv_table.tsv|-> <Q_dir> <Q_prefix> <fam> <out_prefix> [table]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture


def read_eigenvec(path):
    df = pd.read_csv(path, sep=r"\s+", header=None, comment=None)
    first = str(df.iloc[0, 0])
    if first in ("#FID", "FID", "#IID", "IID"):  # plink2 header
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
        if "FID" in df.columns:
            df = df.drop(columns="FID")
        ids = df.iloc[:, 0].astype(str).values
        return ids, df.iloc[:, 1:].values.astype(float)
    return df.iloc[:, 1].astype(str).values, df.iloc[:, 2:].values.astype(float)


def pick_gmm(X, kmax, seed):
    rows, best = [], None
    for cov in ("full", "diag"):
        for k in range(1, kmax + 1):
            if k > len(X) // 3:
                break
            g = GaussianMixture(k, covariance_type=cov, n_init=10, random_state=seed, reg_covar=1e-6).fit(X)
            rows.append({"covariance": cov, "k": k, "BIC": g.bic(X)})
            if best is None or rows[-1]["BIC"] < best[0]:
                best = (rows[-1]["BIC"], g)
    return best[1], pd.DataFrame(rows)


def admixture_labels(qdir, qprefix, cv_path, ids, admix_min):
    qdir = Path(qdir)
    if cv_path and Path(cv_path).exists():
        cv = pd.read_csv(cv_path, sep="\t")
        k = int(cv.loc[cv["CV"].idxmin(), "K"])
    else:
        ks = [int(m.group(1)) for p in qdir.glob(f"{qprefix}.*.Q") if (m := re.search(r"\.(\d+)\.Q$", str(p)))]
        if not ks:
            return None, None, None
        k = max(ks)
    q = np.loadtxt(qdir / f"{qprefix}.{k}.Q")
    if len(q) != len(ids):
        sys.exit(f"ERROR: {k}.Q has {len(q)} rows but {len(ids)} samples")
    top = q.argmax(1)
    lab = np.where(q.max(1) >= admix_min, [f"A{c + 1}" for c in top], "admixed")
    return lab, k, q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("eigenvec")
    ap.add_argument("eigenval")
    ap.add_argument("cv")
    ap.add_argument("qdir")
    ap.add_argument("qprefix")
    ap.add_argument("fam")
    ap.add_argument("out")
    ap.add_argument("table", nargs="?")
    ap.add_argument("--npcs", type=int, default=10)
    ap.add_argument("--kmax-gmm", type=int, default=10)
    ap.add_argument("--min-size", type=int, default=5)
    ap.add_argument("--min-post", type=float, default=0.9)
    ap.add_argument("--admix-min", type=float, default=0.6)
    ap.add_argument("--min-cluster-size", type=int, default=5)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    ids, pcs = read_eigenvec(a.eigenvec)
    X = pcs[:, : a.npcs]
    print(f"{len(ids)} samples, {X.shape[1]} PCs")

    gmm, bic = pick_gmm(X, a.kmax_gmm, a.seed)
    bic.to_csv(f"{a.out}.gmm_bic.tsv", sep="\t", index=False)
    post = gmm.predict_proba(X)
    g_lab = post.argmax(1)
    g_conf = post.max(1)
    print(f"GMM: k={gmm.n_components} ({gmm.covariance_type})")

    h_lab = HDBSCAN(min_cluster_size=a.min_cluster_size, copy=True).fit_predict(X)
    print(f"HDBSCAN: {len(set(h_lab) - {-1})} clusters, {(h_lab == -1).sum()} noise")

    adm = None
    if a.qdir != "-":
        adm, kbest, q = admixture_labels(a.qdir, a.qprefix, a.cv, ids, a.admix_min)
    if adm is not None:
        print(f"ADMIXTURE: best K={kbest}")

    # final = GMM clusters, small ones removed, ordered by size and renamed SP1, SP2, ...
    sizes = pd.Series(g_lab).value_counts()
    order = [c for c in sizes.index if sizes[c] >= a.min_size]
    name = {c: f"SP{i + 1}" for i, c in enumerate(order)}
    final = np.array([name.get(c, "unassigned") for c in g_lab])
    low = g_conf < a.min_post

    out = pd.DataFrame({"IID": ids, "subpop": final, "gmm_posterior": g_conf.round(4), "low_confidence": low,
                        "hdbscan": np.where(h_lab == -1, "noise", ["H" + str(c + 1) for c in h_lab])})
    if adm is not None:
        out["admixture"] = adm
    for i in range(X.shape[1]):
        out[f"PC{i + 1}"] = X[:, i]
    if adm is not None:
        for i in range(q.shape[1]):
            out[f"Q{i + 1}"] = q[:, i]

    if a.table and Path(a.table).exists():
        t = pd.read_excel(a.table) if a.table.endswith((".xlsx", ".xls")) else pd.read_csv(a.table, sep=None, engine="python")
        for col in ("sample_id", "master_id"):
            if col in t.columns:
                t[col] = t[col].astype(str)
        key = next((c for c in ("sample_id", "master_id") if c in t.columns and t[c].isin(ids).any()), None)
        if key:
            keep = [c for c in ("self_reported_label", "empirical_cluster", "sampling_location", "country", "latitude", "longitude") if c in t.columns]
            out = out.merge(t[[key] + keep].drop_duplicates(key).rename(columns={key: "IID"}), on="IID", how="left")

    out.to_csv(f"{a.out}.subpops.tsv", sep="\t", index=False)

    # agreement between methods
    ok = final != "unassigned"
    lines = ["method comparison (adjusted Rand index; 1 = identical partitions):"]
    pairs = [("GMM", "HDBSCAN", g_lab, h_lab)]
    if adm is not None:
        pairs += [("GMM", "ADMIXTURE", g_lab, adm), ("HDBSCAN", "ADMIXTURE", h_lab, adm)]
    for n1, n2, l1, l2 in pairs:
        lines.append(f"  {n1} vs {n2}: {adjusted_rand_score(l1, l2):.3f}")
    summ = out.groupby("subpop").agg(n=("IID", "size"), mean_posterior=("gmm_posterior", "mean"),
                                    n_low_confidence=("low_confidence", "sum")).reset_index()
    summ.to_csv(f"{a.out}.subpop_summary.tsv", sep="\t", index=False)
    Path(f"{a.out}.method_agreement.txt").write_text("\n".join(lines) + "\n\n" + summ.to_string(index=False) + "\n")
    print("\n".join(lines))
    print(summ.to_string(index=False))

    # plots
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n_pairs = min(2, X.shape[1] // 2)
    fig, axs = plt.subplots(1, n_pairs, figsize=(6 * n_pairs, 5), squeeze=False)
    cmap = plt.get_cmap("tab10")
    labs = sorted(set(final))
    for p in range(n_pairs):
        ax = axs[0][p]
        for i, l in enumerate(labs):
            m = final == l
            ax.scatter(X[m, 2 * p], X[m, 2 * p + 1], s=22, c="lightgrey" if l == "unassigned" else [cmap(i % 10)], label=f"{l} (n={m.sum()})")
        ax.scatter(X[low, 2 * p], X[low, 2 * p + 1], s=70, facecolors="none", edgecolors="k", linewidths=0.8)
        ax.set_xlabel(f"PC{2 * p + 1}")
        ax.set_ylabel(f"PC{2 * p + 2}")
    axs[0][0].legend(fontsize=8, title="circled = low confidence")
    fig.tight_layout()
    fig.savefig(f"{a.out}.subpops_PCA.png", dpi=200)
    fig.savefig(f"{a.out}.subpops_PCA.pdf")


if __name__ == "__main__":
    main()
