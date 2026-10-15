#!/usr/bin/env python3
"""Choose the fewest samples to remove so no pair is related up to a given degree.

Reads KING --related output (.kin and/or .kin0), builds a graph whose edges are
pairs with Kinship >= 2^-(degree+1.5), and finds a minimum vertex cover per
connected component (exact for small components, greedy + redundancy pruning
for large ones). Ties are broken in favour of removing the sample with the
higher missing-genotype rate.

usage: prune_related.py DEGREE IMISS KIN [KIN ...] -- OUT_PREFIX
"""
import sys
from collections import defaultdict
from itertools import combinations

EXACT_MAX = 18

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

kinship = {}
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
        for line in fh:
            p = line.split()
            if not p:
                continue
            kv = float(p[k])
            if kv >= thr:
                a, b = (p[f1], p[i1]), (p[f2], p[i2])
                kinship[frozenset((a, b))] = kv

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

removed = set()
for comp in comps:
    edges = [(a, b) for a in comp for b in adj[a] if a < b]
    if len(comp) <= EXACT_MAX:
        for k in range(1, len(comp) + 1):
            covers = [set(s) for s in combinations(comp, k)
                      if all(u in s or v in s for u, v in edges)]
            if covers:
                removed |= max(covers, key=lambda s: (sum(fm(n) for n in s), sorted(s)))
                break
    else:
        rem, left = set(), set(edges)
        while left:
            deg = defaultdict(int)
            for u, v in left:
                deg[u] += 1
                deg[v] += 1
            pick = max(deg, key=lambda n: (deg[n], fm(n), n))
            rem.add(pick)
            left = {(u, v) for u, v in left if u != pick and v != pick}
        for n in sorted(rem, key=lambda n: (fm(n), n)):   # drop redundant removals
            if all(m in rem for m in adj[n]):
                rem.discard(n)
        removed |= rem

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
print("related pairs (kinship >= %.4f): %d; samples involved: %d; removed: %d"
      % (thr, len(kinship), len(adj), len(removed)))
