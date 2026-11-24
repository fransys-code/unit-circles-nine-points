"""Upper bound, step 1: every 9-point structure with >= 21 triples that passes the rhombus filter
is one of the 22 structures of data/structures_9_points_21_triples.json, and none has >= 22 triples.
Needs pynauty (canonical forms); filter: rhombus.py, or the faster rhombus.cpp if its executable
is given.

A structure (k, T): k points, T a set of triples, each pair of points in at most 2 triples (two
distinct unit circles meet in at most 2 points). Removing a point of minimum degree d <= 3t/k from
a realisable structure with k points and t triples leaves a realisable structure with k - 1 points
and t - d triples. So the level X(k) of filtered structures with >= tmin(k) triples is obtained
from X(k - 1) by adding a point of minimum degree (its triples form a graph of maximum degree 2 on
the old points), provided t - floor(3t/k) >= tmin(k - 1) for every t >= tmin(k) (asserted).
Level 5 is built by brute force.

usage: python3 generate.py [rhombus_executable] [processes]
"""
import itertools
import json
import os
import subprocess
import sys
from multiprocessing import Pool

import pynauty

TMIN = {5: 3, 6: 6, 7: 9, 8: 14, 9: 20}


def canon_key(k, triples):
    m = len(triples)
    adj = {k + j: list(t) for j, t in enumerate(triples)}
    g = pynauty.Graph(k + m, directed=False, adjacency_dict=adj,
                      vertex_coloring=[set(range(k)), set(range(k, k + m))] if m else [set(range(k))])
    return pynauty.certificate(g)


def pair_mult(triples):
    mult = {}
    for t in triples:
        for p in itertools.combinations(sorted(t), 2):
            mult[p] = mult.get(p, 0) + 1
    return mult


def base_level(k, tmin):
    alltr = list(itertools.combinations(range(k), 3))
    out = {}
    for mask in range(1 << len(alltr)):
        T = [alltr[i] for i in range(len(alltr)) if mask >> i & 1]
        if len(T) < tmin or max(pair_mult(T).values(), default=0) > 2:
            continue
        out.setdefault(canon_key(k, T), T)
    return list(out.values())


def extensions(args):
    kold, T, tmin = args
    mult = pair_mult(T)
    deg = [0] * kold
    for t in T:
        for p in t:
            deg[p] += 1
    pairs = [p for p in itertools.combinations(range(kold), 2) if mult.get(p, 0) <= 1]
    res = {}
    dg = [0] * kold
    chosen = []

    def rec(i):
        d = len(chosen)
        if len(T) + d >= tmin and all(deg[u] + dg[u] >= d for u in range(kold)):
            S = [tuple(t) for t in T] + [(a, b, kold) for a, b in chosen]
            res.setdefault(canon_key(kold + 1, S), S)
        for i2 in range(i, len(pairs)):
            a, b = pairs[i2]
            if dg[a] < 2 and dg[b] < 2:
                dg[a] += 1; dg[b] += 1; chosen.append((a, b))
                rec(i2 + 1)
                chosen.pop(); dg[a] -= 1; dg[b] -= 1

    rec(0)
    return res


def py_filter(args):
    from rhombus import passes
    k, S = args
    return passes(k, S)


def run_filter(k, structs, exe, pool):
    if exe is None:
        return pool.map(py_filter, [(k, S) for S in structs], chunksize=64)
    inp = "".join("%d %d %s\n" % (k, len(S), " ".join("%d %d %d" % tuple(t) for t in S)) for S in structs)
    out = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split()
    assert len(out) == len(structs)
    return [x == "1" for x in out]


def hist(structs):
    h = {}
    for S in structs:
        h[len(S)] = h.get(len(S), 0) + 1
    return dict(sorted(h.items()))


if __name__ == "__main__":
    exe = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != "-" else None
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    for k in range(6, 10):
        for t in range(TMIN[k], k * (k - 1) // 3 + 1):
            assert t - (3 * t) // k >= TMIN[k - 1], (k, t)
    with Pool(procs) as pool:
        cand = base_level(5, TMIN[5])
        level = [S for S, ok in zip(cand, run_filter(5, cand, exe, pool)) if ok]
        print(5, "candidates", hist(cand), "pass the filter", hist(level), flush=True)
        for k in range(6, 10):
            found = {}
            for res in pool.imap_unordered(extensions, [(k - 1, S, TMIN[k]) for S in level], chunksize=8):
                for key, S in res.items():
                    found.setdefault(key, S)
            cand = list(found.values())
            level = [S for S, ok in zip(cand, run_filter(k, cand, exe, pool)) if ok]
            print(k, "candidates", hist(cand), "pass the filter", hist(level), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    ref = json.load(open(os.path.join(here, "data", "structures_9_points_21_triples.json")))["structures"]
    ref_keys = {canon_key(9, [tuple(t) for t in s["triples"]]) for s in ref}
    mine = {canon_key(9, S) for S in level if len(S) >= 21}
    assert all(len(S) <= 21 for S in level), "a structure with >= 22 triples passes the filter"
    assert len(ref_keys) == 22 and mine == ref_keys, "survivors differ from the data file"
    print("9 points: 0 structures with >= 22 triples, the 22 with 21 triples are those of the data file")
    print("UPPER BOUND STEP 1 VERIFIED")
