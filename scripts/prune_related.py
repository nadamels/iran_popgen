#!/usr/bin/env python3
"""Choose the fewest samples to remove so no pair is related up to a given degree.

Reads KING --related output (.kin and/or .kin0), builds a graph whose edges are
pairs that KING itself calls related up to DEGREE (its InfType column: Dup/MZ, PO, FS, 2nd, 3rd - this
is based on shared IBD segments and can disagree with the Kinship column) OR have
Kinship >= 2^-(degree+1.5), and finds a minimum vertex cover per
connected component: brute force up to 18 samples, an exact MILP (scipy >= 1.9) for larger clusters,
and only if scipy is missing / the MILP times out a greedy + redundancy-pruning fallback, in which case
a WARNING says the removal set is not guaranteed minimal. Ties are broken in favour of removing the
sample with the higher missing-genotype rate.

usage: prune_related.py DEGREE IMISS KIN [KIN ...] -- OUT_PREFIX
"""
import re
import sys
from collections import defaultdict
from itertools import combinations

EXACT_MAX = 18          # components up to this size are solved by brute force (no scipy needed)
MILP_SECONDS = 120      # time limit per component for the exact (scipy MILP) solver

degree = int(sys.argv[1])
imiss = sys.argv[2]
sep = sys.argv.index("--")
kin_files = sys.argv[3:sep]
out = sys.argv[sep + 1]
thr = 2 ** (-(degree + 1.5))

fmiss = {}
with open(imiss) as fh:
    h = fh.readline().split()
    fi, ii, mi = h.index("FID"), h.index("IID"), h.index("F_MISS")
    for line in fh:
        p = line.split()
        if p:
            try:
                fmiss[(p[fi], p[ii])] = float(p[mi])
            except ValueError:
                fmiss[(p[fi], p[ii])] = 0.0

# KING's InfType -> degree. Matching is tolerant of spelling/case ("Dup/MZ", "MZ twin", "MZTwin", "1st", ...)
def inf_degree(label):
    lab = re.sub(r"[^a-z0-9]", "", label.lower())
    if "dup" in lab or "mz" in lab or "twin" in lab:
        return 0
    if lab in ("po", "fs", "1st", "1stdegree"):
        return 1
    m = re.fullmatch(r"([2-9])(?:st|nd|rd|th)(?:degree)?", lab)
    if m:
        return int(m.group(1))
    return None            # "UN"/"unrelated"/unknown

unknown_labels = set()
kinship = {}
n_by_inftype = 0
for path in kin_files:
    try:
        fh = open(path)
    except FileNotFoundError:
        continue
    with fh:
        h = fh.readline().split()
        if not h:
            continue
        k = h.index("Kinship")
        i1, i2 = h.index("ID1"), h.index("ID2")
        f1 = h.index("FID1") if "FID1" in h else h.index("FID")
        f2 = h.index("FID2") if "FID2" in h else h.index("FID")
        it = h.index("InfType") if "InfType" in h else None
        for line in fh:
            p = line.split()
            if not p:
                continue
            kv = float(p[k])
            by_inf = False
            if it is not None and len(p) > it:
                d_inf = inf_degree(p[it])
                by_inf = d_inf is not None and d_inf <= degree
                if d_inf is None and re.sub(r"[^a-z]", "", p[it].lower()) not in ("un", "unrelated", "other"):
                    unknown_labels.add(p[it])
            if kv >= thr or by_inf:
                if by_inf and kv < thr:
                    n_by_inftype += 1
                a, b = (p[f1], p[i1]), (p[f2], p[i2])
                kinship[frozenset((a, b))] = kv

if unknown_labels:
    print("NOTE: unrecognised InfType label(s) ignored: " + ", ".join(sorted(unknown_labels)))

adj = defaultdict(set)
for pair in kinship:
    a, b = tuple(pair)
    adj[a].add(b)
    adj[b].add(a)

seen, comps = set(), []
for n in sorted(adj):
    if n in seen:
        continue
    stack, comp = [n], []
    seen.add(n)
    while stack:
        x = stack.pop()
        comp.append(x)
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    comps.append(sorted(comp))

def fm(n):
    return fmiss.get(n, 0.0)

def solve_exact(comp, edges):
    """Minimum vertex cover by MILP (scipy >= 1.9). Ties go to the sample with more missing data.
    Returns (set, proven_optimal) or None if scipy's MILP is unavailable."""
    try:
        import numpy as np
        from scipy.optimize import Bounds, LinearConstraint, milp
        from scipy.sparse import lil_matrix
    except Exception:
        return None
    idx = {x: i for i, x in enumerate(comp)}
    n = len(comp)
    eps = 1e-6 / max(n, 1)                          # a bonus too small to ever beat one fewer removal
    c = np.array([1.0 - eps * fm(x) for x in comp])
    A = lil_matrix((len(edges), n))
    for k, (u, v) in enumerate(edges):
        A[k, idx[u]] = 1
        A[k, idx[v]] = 1
    res = milp(c, constraints=LinearConstraint(A.tocsr(), lb=1, ub=np.inf), integrality=np.ones(n),
               bounds=Bounds(0, 1), options={"time_limit": MILP_SECONDS})
    if res.x is None:
        return None
    chosen = {comp[i] for i in range(n) if res.x[i] > 0.5}
    if not all(u in chosen or v in chosen for u, v in edges):
        return None
    return chosen, res.status == 0


def solve_greedy(comp, edges):
    rem, left = set(), set(edges)
    while left:
        deg = defaultdict(int)
        for u, v in left:
            deg[u] += 1
            deg[v] += 1
        pick = max(deg, key=lambda x: (deg[x], fm(x), x))
        rem.add(pick)
        left = {(u, v) for u, v in left if u != pick and v != pick}
    for x in sorted(rem, key=lambda x: (fm(x), x)):   # drop redundant removals
        if all(m in rem for m in adj[x]):
            rem.discard(x)
    return rem


removed = set()
not_proven = []
for comp in comps:
    edges = [(a, b) for a in comp for b in adj[a] if a < b]
    if len(comp) <= EXACT_MAX:
        for k in range(1, len(comp) + 1):
            covers = [set(s) for s in combinations(comp, k)
                      if all(u in s or v in s for u, v in edges)]
            if covers:
                removed |= max(covers, key=lambda s: (sum(fm(x) for x in s), sorted(s)))
                break
    else:
        res = solve_exact(comp, edges)
        if res is not None:
            rem, optimal = res
            if not optimal:
                not_proven.append((len(comp), "MILP time limit reached; best solution found"))
        else:
            rem = solve_greedy(comp, edges)
            not_proven.append((len(comp), "scipy>=1.9 not available; greedy result"))
        removed |= rem
for size, why in not_proven:
    print("WARNING: related cluster of %d samples: removal set NOT guaranteed minimal (%s)" % (size, why),
          file=sys.stderr)

with open(out + ".removed_related.txt", "w") as fh:
    for n in sorted(removed):
        fh.write("%s %s\n" % n)
with open(out + ".removed_related_details.tsv", "w") as fh:
    fh.write("FID\tIID\tn_related\tmax_kinship\trelated_kept_samples\n")
    for n in sorted(removed):
        kept = sorted(m for m in adj[n] if m not in removed)
        mk = max(kinship[frozenset((n, m))] for m in adj[n])
        fh.write("%s\t%s\t%d\t%.4f\t%s\n" % (n[0], n[1], len(adj[n]), mk,
                                            ",".join(m[1] for m in kept) or "-"))
print("related pairs (KING InfType up to degree %d, or kinship >= %.4f): %d (%d of them only by InfType); "
      "samples involved: %d; removed: %d"
      % (degree, thr, len(kinship), n_by_inftype, len(adj), len(removed)))
