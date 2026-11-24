"""The last structure S* has no realisation (Python standard library only). See PROOF.md, section 3.

Part 1 (exact linear algebra): S* of the data file is S* of the model; modulo the 27 rhombus
identities, (1) the 54 point-centre vectors of the 18 transversal circles are +-E[s][l] (9 values);
(2) E[s][l] - E[s'][l'] is +-(a difference of two points or of two centres) when s != s';
(3) the three rows E[s][0] + E[s][1] + E[s][2] have the same sum; (4) for each group circle,
with xi_l = (point of the group) - (its centre), the group relations of group j hold, the three
group circles giving j = 0, 1, 2.
Part 2 (numerical sanity check of the lemma, not part of the proof): pairs of unit triples with the
same sum, computed as roots of z^3 - H z^2 + p conj(H) z - p, all satisfy the lemma.
Part 3 (enumeration): every weak order type of the 9 vectors E compatible with the lemma on the
row relations (105,660 of them) admits, for some j, no placement of xi_0, xi_1, xi_2 compatible
with the lemma on the 9 relations of group j. The excluded group of each type is compared with
data/sstar_order_types.jsonl.gz.

usage: python3 verify_sstar.py [processes]        (part 3: about 40 CPU-minutes)
       python3 verify_sstar.py --quick            (parts 1 and 2 only, seconds)
"""
import cmath
import gzip
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------- the lemma (PROOF.md, section 3.3) ----------
# Positions on the circle are integers modulo F = 2M; the antipode of position q is q + M.

def gaps(t, F):
    s = sorted(x % F for x in t)
    return [s[1] - s[0], s[2] - s[1], F - (s[2] - s[0])], s


def short_arc(t, M):
    """For a triple inside an open half-circle: (first point, offset of middle, offset of last)."""
    F = 2 * M
    g, s = gaps(t, F)
    i = max(range(3), key=lambda n: g[n])
    start = s[(i + 1) % 3]
    o = sorted((x - start) % F for x in t)
    return start, o[1], o[2]


def mixed(A, B, F):
    """Counter-clockwise order A0 < B0 <= B1 < A1 <= A2 < B2 within one turn."""
    o = (B[0] - A[0]) % F
    return 0 < o and o + B[1] < A[1] and A[2] < o + B[2] < F


def lemma(a, b, M):
    """Necessary condition for two triples of unit vectors at positions a, b to have equal sums."""
    F = 2 * M
    a = [x % F for x in a]
    b = [x % F for x in b]
    if sorted(a) == sorted(b):
        return True
    for p in set(a) & set(b):
        ra, rb = list(a), list(b)
        ra.remove(p)
        rb.remove(p)
        if (ra[0] - ra[1]) % F == M and (rb[0] - rb[1]) % F == M:
            return True
    ga, sa = gaps(a, F)
    gb, sb = gaps(b, F)
    kind_a = "acute" if max(ga) < M else "right" if max(ga) == M else "obtuse"
    kind_b = "acute" if max(gb) < M else "right" if max(gb) == M else "obtuse"
    if kind_a != kind_b or kind_a == "right":
        return False
    if kind_a == "acute":
        if len(set(a) | set(b)) < 6:
            return False
        return sum(sa[0] < q < sa[1] for q in sb) == 1 and sum(sa[1] < q < sa[2] for q in sb) == 1
    A, B = short_arc(a, M), short_arc(b, M)
    return mixed(A, B, F) or mixed(B, A, F)


# ---------- part 1: the reduction ----------

def model():
    tr = [tuple(3 * j + t for t in range(3)) for j in range(3)]
    for t in itertools.product(range(3), repeat=3):
        if sum(t) % 3:
            tr.append((t[0], 3 + t[1], 6 + t[2]))
    return tr


def row_relations():
    out = []
    for r, r2 in ((0, 1), (0, 2), (1, 2)):
        out.append([(3 * r + m, 1) for m in range(3)] + [(3 * r2 + m, -1) for m in range(3)])
    return out


def group_relations(j):
    """xi_l - xi_{l-1} = sum_{s != r} E[s][l+js] - sum_{m != l+jr} E[r][m]   (labels: E[s][l] = 3s+l, xi_l = 9+l)."""
    out = []
    for l in range(3):
        for r in range(3):
            rel = [(9 + l, 1), (9 + (l - 1) % 3, -1)]
            rel += [(3 * r + m, 1) for m in range(3) if m != (l + j * r) % 3]
            rel += [(3 * s + (l + j * s) % 3, -1) for s in range(3) if s != r]
            out.append(rel)
    return out


class Quotient:
    def __init__(self, gens):
        self.rows = []
        for g in gens:
            self.add(g)

    def reduce(self, v):
        v = [Fraction(x) for x in v]
        for piv, row in self.rows:
            if v[piv]:
                c = v[piv]
                v = [a - c * b for a, b in zip(v, row)]
        return v

    def add(self, g):
        v = self.reduce(g)
        for i, x in enumerate(v):
            if x:
                v = [a / x for a in v]
                self.rows = [(p, [a - r[i] * b for a, b in zip(r, v)]) for p, r in self.rows]
                self.rows.append((i, v))
                return

    def zero(self, v):
        return not any(self.reduce(v))


def part1():
    data = json.load(open(os.path.join(HERE, "data", "structures_9_points_21_triples.json")))["structures"]
    star = [s for s in data if "relabelling_to_model" in s]
    assert len(star) == 1
    perm = star[0]["relabelling_to_model"]
    mod = model()
    assert sorted(perm) == list(range(9))
    assert {frozenset(perm[p] for p in t) for t in star[0]["triples"]} == {frozenset(t) for t in mod}
    print("S* of the data file = model (points 0,1,2 = X, 3,4,5 = A, 6,7,8 = B; 3 group circles and the "
          "18 transversals (x_i, a_j, b_k) with i + j + k != 0 mod 3)")

    dim = 9 + 21
    cidx = {frozenset(t): 9 + i for i, t in enumerate(mod)}

    def vec(terms):
        v = [0] * dim
        for i, c in terms:
            v[i] += c
        return v

    gens = []
    for t1, t2 in itertools.combinations(mod, 2):
        common = set(t1) & set(t2)
        assert len(common) <= 2
        if len(common) == 2:
            p, q = common
            gens.append(vec([(p, 1), (q, 1), (cidx[frozenset(t1)], -1), (cidx[frozenset(t2)], -1)]))
    Q = Quotient(gens)
    print("rhombus identities p + q = M + M':", len(gens), " rank:", len(Q.rows))

    def circle(i, j):
        return cidx[frozenset((i, 3 + j, 6 + (1 - i - j) % 3))]

    E = {(s, l): vec([(l, 1), (circle(l, (l - s) % 3), -1)]) for s in range(3) for l in range(3)}

    def same(u, v, sign=1):
        return Q.zero([x - sign * y for x, y in zip(u, v)])

    used = {key: 0 for key in E}
    for t in mod[3:]:
        for p in t:
            inc = vec([(p, 1), (cidx[frozenset(t)], -1)])
            hits = [(key, sg) for key in E for sg in (1, -1) if same(inc, E[key], sg)]
            assert len(hits) == 1, hits
            used[hits[0][0]] += 1
    assert all(n == 6 for n in used.values())
    print("(1) the 54 point - centre vectors of the transversal circles are +-E[s][l], each E six times")

    diffs = [vec([(p, 1), (q, -1)]) for p, q in itertools.combinations(range(9), 2)]
    diffs += [vec([(p, 1), (q, -1)]) for p, q in itertools.combinations(range(9, 30), 2)]
    for k1, k2 in itertools.combinations(E, 2):
        d = [x - y for x, y in zip(E[k1], E[k2])]
        proved = any(same(d, w, sg) for w in diffs for sg in (1, -1))
        assert proved == (k1[0] != k2[0]), (k1, k2)
    print("(2) E[s][l] - E[s'][l'] = +-(two points or two centres) for s != s', so E[s][l] != E[s'][l']")

    def lin(rel, xi):
        v = [0] * dim
        for lab, sg in rel:
            u = E[lab // 3, lab % 3] if lab < 9 else xi[lab - 9]
            v = [x + sg * y for x, y in zip(v, u)]
        return v

    assert all(Q.zero(lin(rel, None)) for rel in row_relations())
    print("(3) the three rows of E have the same sum")
    found = {}
    for g in range(3):
        o = cidx[frozenset(mod[g])]
        for j in range(3):
            for pi in itertools.permutations(range(3)):
                xi = [vec([(3 * g + pi[l], 1), (o, -1)]) for l in range(3)]
                if all(Q.zero(lin(rel, xi)) for rel in group_relations(j)):
                    found.setdefault(g, set()).add(j)
    assert sorted(found) == [0, 1, 2] and all(len(v) == 1 for v in found.values())
    assert sorted(next(iter(v)) for v in found.values()) == [0, 1, 2]
    print("(4) group relations: circle X gives j = %d, circle A gives j = %d, circle B gives j = %d"
          % tuple(next(iter(found[g])) for g in range(3)))
    print("PART 1 (reduction) VERIFIED")


# ---------- part 2: numerical sanity check of the lemma ----------

def cubic_roots(c2, c1, c0):
    z = [complex(0.4, 0.9) ** k for k in range(3)]
    for _ in range(500):
        new = []
        for i in range(3):
            num = z[i] ** 3 + c2 * z[i] ** 2 + c1 * z[i] + c0
            den = 1
            for k in range(3):
                if k != i:
                    den *= z[i] - z[k]
            new.append(z[i] - num / den)
        z = new
    return z


def to_grid(angles, tol=1e-7):
    items = []
    for th in angles:
        v = (th - angles[0]) % (2 * math.pi)
        if v > 2 * math.pi - tol:
            v = 0.0
        half = 0
        if v > math.pi - tol:
            half, v = 1, max(v - math.pi, 0.0)
        items.append((v, half))
    cls = []
    for h in sorted(v for v, _ in items):
        if not cls or h - cls[-1] > tol:
            cls.append(h)
    M = len(cls)
    return [min(range(M), key=lambda i: abs(cls[i] - v)) + half * M for v, half in items], M


def part2(trials=20000, seed=1):
    rng = random.Random(seed)
    unit = lambda: cmath.exp(1j * rng.uniform(0, 2 * math.pi))
    stats = {}
    for _ in range(trials):
        kind = rng.randrange(4)
        if kind == 0:
            T = [unit() for _ in range(3)]
        elif kind == 1:
            u = unit()
            T = [u, -u, unit()]
        elif kind == 2:
            u = unit()
            T = [u, u, unit()]
        else:
            a = unit()
            T = [a, a * cmath.exp(1j * rng.choice([1, 2, 3]) * math.pi / 3), unit()]
        H = sum(T)
        p = T[0] * T[1] * T[2]
        if rng.random() < 0.5:
            p2 = p * cmath.exp(1j * rng.gauss(0, rng.choice([0.05, 0.5, 2.0])))
        else:
            p2 = p
        T2 = cubic_roots(-H, p2 * H.conjugate(), -p2)
        if max(abs(abs(z) - 1) for z in T2) > 1e-9 or abs(sum(T2) - H) > 1e-9:
            continue
        zs = T + T2
        ok = True
        for i in range(6):
            for k in range(i):
                d = min(abs(zs[i] - zs[k]), abs(zs[i] + zs[k]))
                if 1e-12 < d < 1e-5:
                    ok = False
        if not ok:
            continue
        pos, M = to_grid([cmath.phase(z) for z in zs])
        assert lemma(pos[:3], pos[3:], M), (T, T2)
        h = abs(H)
        label = "|H| < 1" if h < 1 - 1e-9 else "|H| = 1" if h < 1 + 1e-9 else "1 < |H| < 3" if h < 3 - 1e-9 else "|H| = 3"
        stats[label] = stats.get(label, 0) + 1
    print("pairs of triples with equal sums, all accepted by the lemma:", dict(sorted(stats.items())))
    print("PART 2 (sanity check) PASSED")


# ---------- part 3: weak order types ----------
# A weak order type is a list of m "slots" covering a half-turn [0, pi); slot i holds the signed
# vectors (v, +1) or (v, -1) whose direction (v or -v) is at that angle. Slot i + m holds the
# opposites of slot i, so the position of a signed vector is an integer modulo 2m.

SPLITS = [A for A in itertools.combinations(range(6), 3) if 0 in A]


def positions(slots):
    m = len(slots)
    pos = {}
    for i, sl in enumerate(slots):
        for v, s in sl:
            pos[v] = i if s == 1 else i + m
    return pos


def insertions(slots, v):
    for i in range(len(slots)):
        for s in (1, -1):
            yield slots[:i] + [slots[i] + [(v, s)]] + slots[i + 1:]
            yield slots[:i + 1] + [[(v, s)]] + slots[i + 1:]


class Checker:
    def __init__(self, relations, distinct):
        self.rels = []
        for rel in relations:
            labels = frozenset(v for v, _ in rel)
            for A in SPLITS:
                T = [rel[i] for i in A]
                T2 = [(v, -s) for i, (v, s) in enumerate(rel) if i not in A]
                self.rels.append((labels, T, T2))
        self.distinct = distinct

    def ok(self, slots, new):
        m = len(slots)
        pos = positions(slots)
        for grp in self.distinct:
            if new in grp:
                ps = [pos[v] for v in grp if v in pos]
                if len(set(ps)) < len(ps):
                    return False
        F = 2 * m
        for labels, T, T2 in self.rels:
            if new in labels and all(v in pos for v in labels):
                pT = [(pos[v] + (0 if s == 1 else m)) % F for v, s in T]
                pT2 = [(pos[v] + (0 if s == 1 else m)) % F for v, s in T2]
                if not lemma(pT, pT2, m):
                    return False
        return True


def extend(slots, labels, checker, out=None, first_only=False):
    if not labels:
        if out is not None:
            out.append(slots)
        return 1
    n = 0
    for s2 in insertions(slots, labels[0]):
        if checker.ok(s2, labels[0]):
            n += extend(s2, labels[1:], checker, out, first_only)
            if first_only and n:
                return n
    return n


XI = [9, 10, 11]
DISTINCT_E = [[i, k] for i in range(9) for k in range(i) if i // 3 != k // 3]
CK_ROWS = Checker(row_relations(), DISTINCT_E)
CK_GROUP = [Checker(group_relations(j), [XI]) for j in range(3)]


def excluded_group(slots):
    for j in range(3):
        if extend(slots, XI, CK_GROUP[j], first_only=True) == 0:
            return j
    return None


def canonical(slots):
    return json.dumps([sorted([list(x) for x in sl]) for sl in slots], separators=(",", ":"))


def part3(procs):
    t0 = time.time()
    types = []
    extend([[(0, 1)]], list(range(1, 9)), CK_ROWS, out=types)
    print("weak order types of the 9 vectors E compatible with the row relations:", len(types),
          "(%.0f s)" % (time.time() - t0), flush=True)
    table = {}
    with gzip.open(os.path.join(HERE, "data", "sstar_order_types.jsonl.gz"), "rt") as f:
        for line in f:
            d = json.loads(line)
            table[canonical(d["type"])] = d["excluded_group"]
    mine = [canonical(s) for s in types]
    assert len(set(mine)) == len(mine) == len(table) and set(mine) == set(table), "type list differs from the table"
    print("same list as data/sstar_order_types.jsonl.gz", flush=True)
    counts = {0: 0, 1: 0, 2: 0, None: 0}
    agree = 0
    with Pool(procs) as pool:
        for i, (key, j) in enumerate(zip(mine, pool.imap(excluded_group, types, chunksize=8))):
            counts[j] += 1
            agree += table[key] == j
            if (i + 1) % 10000 == 0:
                print(i + 1, counts, "(%.0f s)" % (time.time() - t0), flush=True)
    print("first excluded group: j=0: %d, j=1: %d, j=2: %d, none: %d" % (counts[0], counts[1], counts[2], counts[None]))
    print("agreement with the table:", agree, "/", len(mine))
    assert counts[None] == 0 and agree == len(mine)
    print("PART 3 (enumeration) VERIFIED: S* has no realisation (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    part1()
    part2()
    if "--quick" not in sys.argv:
        args = [a for a in sys.argv[1:] if not a.startswith("--")]
        part3(int(args[0]) if args else os.cpu_count())
