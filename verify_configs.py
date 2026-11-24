"""Lower bounds: re-check the certified configurations (needs mpmath).

For each file in data/config_*.json: every point-centre incidence vector is, as an exact rational
identity, +-w_i or +-E_j with E_j in `equations`; the interval Krawczyk test proves that the system
|E_j(theta)|^2 = 1 has an exact zero in the box; with the angles in that box, the points are
pairwise distinct and the centres are pairwise distinct (interval arithmetic, 60 digits).
Conclusion: k distinct points and that many distinct unit circles, each through 3 of the points.

usage: python3 verify_configs.py [data/config_*.json ...]
"""
import glob
import json
import os
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 60
mp.iv.dps = 60
iv = mp.iv


def to_iv(fr):
    return iv.mpf(fr.numerator) / fr.denominator


def residual_and_jacobian(eqs, th):
    values, jac = [], []
    for q in eqs:
        S = [i for i, a in enumerate(q) if a]
        val = to_iv(sum(q[i] * q[i] for i in S) - 1)
        grad = {i: iv.mpf(0) for i in S}
        for x in range(len(S)):
            for y in range(x + 1, len(S)):
                i, l = S[x], S[y]
                c = to_iv(2 * q[i] * q[l])
                val += c * iv.cos(th[i] - th[l])
                s = c * iv.sin(th[i] - th[l])
                grad[i] -= s
                grad[l] += s
        values.append(val)
        jac.append(grad)
    return values, jac


def check(c):
    k, triples, r = c["points"], [tuple(t) for t in c["triples"]], c["r"]
    A = [[Fraction(a) for a in row] for row in c["A"]]
    eqs = [[Fraction(a) for a in q] for q in c["equations"]]
    n = k + len(triples)
    assert len(A) == n and all(len(row) == r for row in A)
    assert len({frozenset(t) for t in triples}) == len(triples) == c["circles"]
    unit = {tuple(Fraction(int(i == j)) for i in range(r)) for j in range(r)} | {tuple(q) for q in eqs}
    unit |= {tuple(-a for a in u) for u in unit}
    for j, tr in enumerate(triples):
        assert len(set(tr)) == 3 and all(0 <= p < k for p in tr)
        for p in tr:
            assert tuple(a - b for a, b in zip(A[k + j], A[p])) in unit, ("incidence", j, p)
    fixed = {int(i) - 1: Fraction(v) for i, v in c["fixed_angles"].items()}
    free = sorted(int(i) - 1 for i in c["free_angles"])
    assert sorted(list(fixed) + free) == list(range(1, r)) and len(free) == len(eqs)
    rad = mp.mpf(c["box_radius"])
    mid = {i: mp.mpf(c["free_angles"][str(i + 1)]) for i in free}
    X = [iv.mpf(0)] * r
    xm = [iv.mpf(0)] * r
    for i, v in fixed.items():
        X[i] = xm[i] = to_iv(v)
    for i in free:
        X[i] = iv.mpf([mid[i] - rad, mid[i] + rad])
        xm[i] = iv.mpf(mid[i])
    m = len(eqs)
    if m:
        Fm, Jmid = residual_and_jacobian(eqs, xm)
        _, JX = residual_and_jacobian(eqs, X)
        Jn = mp.matrix(m, m)
        for a in range(m):
            for b, i in enumerate(free):
                Jn[a, b] = mp.mpf(Jmid[a].get(i, iv.mpf(0)).mid)
        Y = Jn ** -1
        Yi = iv.matrix(m, m)
        JXm = iv.matrix(m, m)
        for a in range(m):
            for b, i in enumerate(free):
                Yi[a, b] = iv.mpf(Y[a, b])
                JXm[a, b] = JX[a].get(i, iv.mpf(0))
        YF = Yi * iv.matrix(Fm)
        R = iv.eye(m) - Yi * JXm
        for a in range(m):
            acc = xm[free[a]] - YF[a]
            for b in range(m):
                acc += R[a, b] * (X[free[b]] - xm[free[b]])
            assert acc.a > X[free[a]].a and acc.b < X[free[a]].b, "Krawczyk test failed"
    W = [iv.mpc(iv.cos(t), iv.sin(t)) for t in X]
    pos = []
    for x in range(n):
        z = iv.mpc(0)
        for i, a in enumerate(A[x]):
            if a:
                z += to_iv(a) * W[i]
        pos.append(z)
    sep = []
    for lo, hi in ((0, k), (k, n)):
        best = None
        for a in range(lo, hi):
            for b in range(lo, a):
                d = pos[a] - pos[b]
                assert abs(d.real).a > 0 or abs(d.imag).a > 0, ("not distinct", a, b)
                lb = float(max(abs(d.real).a, abs(d.imag).a).mid)
                best = lb if best is None else min(best, lb)
        sep.append(best)
    return k, len(triples), sep


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    paths = sys.argv[1:] or sorted(glob.glob(os.path.join(here, "data", "config_*.json")))
    for path in paths:
        k, m, sep = check(json.load(open(path)))
        print(f"{os.path.basename(path)}: certified, {k} distinct points, {m} distinct unit circles "
              f"through 3 of them (separation about {sep[0]:.3f} between points, {sep[1]:.3f} between centres)")
    print("ALL CONFIGURATIONS CERTIFIED")
