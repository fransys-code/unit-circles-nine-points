"""Exact linear rhombus filter (Python standard library).

Symbols: the k points, then one centre per triple. H = graph of unit distances, initially the
point-centre incidences. Rules, iterated to a fixed point, in exact rational arithmetic on formal
linear combinations of the symbols:
  * rhombus: if x1, x2, x3, x4 is a 4-cycle of H with x1 != x3 and x2 != x4 (as points of the
    plane), then x1 + x3 = x2 + x4 (x1 and x3 are the two intersections of the unit circles
    centred at x2 and x4);
  * closure: if x - y equals +-(an edge vector of H) modulo the relations, xy is a unit distance.
Two symbols are known to be distinct if they are of the same type (two points, or two centres),
or if their difference is a non-zero rational multiple of a unit vector. Contradictions: two
symbols of the same type coincide; an edge of H is zero; two edges of H are proportional with a
ratio other than +-1. passes(k, triples) is False only if the structure is unrealisable.
"""
from fractions import Fraction
from itertools import combinations


def prim(d):
    for a in d:
        if a:
            return tuple(x / a for x in d)
    return None


class Lin:
    def __init__(self, k, triples):
        self.k = k
        self.T = [tuple(t) for t in triples]
        self.N = k + len(self.T)
        self.rows = {}
        self.U = set()
        for j, t in enumerate(self.T):
            for p in t:
                self.U.add((p, k + j))

    def reduce(self, v):
        v = [Fraction(a) for a in v]
        for p, row in self.rows.items():
            c = v[p]
            if c:
                v = [a - c * b for a, b in zip(v, row)]
        return v

    def add(self, g):
        v = self.reduce(g)
        for i, x in enumerate(v):
            if x:
                v = [a / x for a in v]
                for p, row in self.rows.items():
                    c = row[i]
                    if c:
                        self.rows[p] = [a - c * b for a, b in zip(row, v)]
                self.rows[i] = v
                return True
        return False

    def nf(self, x):
        if x in self.rows:
            return tuple(-a if i != x else Fraction(0) for i, a in enumerate(self.rows[x]))
        return tuple(Fraction(int(i == x)) for i in range(self.N))

    def same_type(self, x, y):
        return (x < self.k) == (y < self.k)

    def run(self):
        """False if a contradiction is found (the structure is unrealisable)."""
        N = self.N
        while True:
            nfs = [self.nf(x) for x in range(N)]
            self.nfs = nfs
            for x, y in combinations(range(N), 2):
                if self.same_type(x, y) and nfs[x] == nfs[y]:
                    return False
            ev = set()
            dirs = {}
            for a, b in self.U:
                d = tuple(p - q for p, q in zip(nfs[a], nfs[b]))
                pd = prim(d)
                if pd is None:
                    return False
                q = abs(next(x for x in d if x))
                if dirs.setdefault(pd, q) != q:
                    return False
                ev.add(d)
                ev.add(tuple(-x for x in d))
            changed = False
            dist = [[False] * N for _ in range(N)]
            for x, y in combinations(range(N), 2):
                d = tuple(p - q for p, q in zip(nfs[x], nfs[y]))
                if (x, y) not in self.U and d in ev:
                    self.U.add((x, y))
                    changed = True
                if self.same_type(x, y) or prim(d) in dirs:
                    dist[x][y] = dist[y][x] = True
            adj = [set() for _ in range(N)]
            for a, b in self.U:
                adj[a].add(b)
                adj[b].add(a)
            for x2, x4 in combinations(range(N), 2):
                if not dist[x2][x4]:
                    continue
                com = sorted(adj[x2] & adj[x4])
                for x1, x3 in combinations(com, 2):
                    if dist[x1][x3]:
                        g = [0] * N
                        g[x1] += 1; g[x3] += 1; g[x2] -= 1; g[x4] -= 1
                        if self.add(g):
                            changed = True
            if not changed:
                return True


def passes(k, triples):
    return Lin(k, triples).run()
