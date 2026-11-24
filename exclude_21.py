"""Upper bound, step 2: the 21 structures of data/structures_9_points_21_triples.json other than S*
have no realisation, not even a complex one. Needs Singular (tested with 4.3.2) on the PATH.

After the rhombus filter (rhombus.py), every symbol is a rational combination of free symbols;
the unit edges of H become vectors e(z) = sum c_i z_i in complex unknowns z_i, with conjugates
e(w) = sum c_i w_i (w_i stands for conj(z_i)). One unit edge is fixed to 1 (rotation). Equations:
e(z) e(w) = 1 for every unit edge. The ideal is then saturated, one after the other, by the
linear forms L(z) = difference of two points or of two centres (non-zero in any realisation).
If 1 lies in the ideal, there is no solution with all points distinct and all centres distinct.

usage: python3 exclude_21.py [processes] [timeout_seconds]
"""
import itertools
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction
from math import lcm
from multiprocessing import Pool

from rhombus import Lin


def poly(coefs, names):
    den = 1
    for c in coefs:
        den = lcm(den, Fraction(c).denominator)
    terms = ["(%d)*%s" % (int(c * den), n) for c, n in zip(coefs, names) if c]
    return "+".join(terms) if terms else "0"


def script(k, T):
    L = Lin(k, T)
    assert L.run()
    N = L.N
    free = [x for x in range(N) if x not in L.rows]
    fr = free[1:]
    zs = ["z%d" % i for i in range(len(fr))]
    ws = ["w%d" % i for i in range(len(fr))]

    def coords(x, y):
        return [L.nfs[x][f] - L.nfs[y][f] for f in fr]

    for x, y in itertools.combinations(range(N), 2):
        assert sum(a - b for a, b in zip(L.nfs[x], L.nfs[y])) == 0
    seen, eqs, first = set(), [], None
    for a, b in sorted(L.U):
        c = coords(a, b)
        key = tuple(c) if next(x for x in c if x) > 0 else tuple(-x for x in c)
        if key in seen:
            continue
        seen.add(key)
        if first is None:
            first = key
            eqs.append("%s-1" % poly(key, zs))
            eqs.append("%s-1" % poly(key, ws))
        else:
            eqs.append("(%s)*(%s)-1" % (poly(key, zs), poly(key, ws)))
    sats, seen = [], set()
    for lo, hi in ((0, k), (k, N)):
        for a, b in itertools.combinations(range(lo, hi), 2):
            c = coords(a, b)
            p = next(x for x in c if x)
            key = tuple(x / p for x in c)
            if key not in seen:
                seen.add(key)
                sats.append(poly(key, zs))
    s = 'LIB "elim.lib";\nring R = 0,(%s),dp;\noption(redSB);\n' % ",".join(zs + ws)
    s += "ideal I = %s;\nI = std(I);\n" % ",\n".join(eqs)
    s += 'if (reduce(1,I)==0) { print("EMPTY 0"); quit; }\n'
    for n, p in enumerate(sats, 1):
        s += "def s%d = sat(I, ideal(%s));\n" % (n, p)
        s += 'if (typeof(s%d)=="list") { I = s%d[1]; } else { I = s%d; }\nI = std(I);\n' % (n, n, n)
        s += 'if (reduce(1,I)==0) { print("EMPTY %d"); quit; }\n' % n
    s += 'print("NOT CLOSED"); print(dim(I)); quit;\n'
    return s, len(fr), len(eqs), len(sats)


def run(args):
    k, T, timeout = args
    s, nv, ne, ns = script(k, T)
    with tempfile.NamedTemporaryFile("w", suffix=".sing", delete=False) as f:
        f.write(s)
    try:
        out = subprocess.run(["Singular", "-q", f.name], capture_output=True, text=True,
                             timeout=timeout).stdout.strip().replace("\n", " ")
    except subprocess.TimeoutExpired:
        out = "TIMEOUT"
    finally:
        os.unlink(f.name)
    return out, "%d complex unknowns, %d equations, %d saturations available" % (nv, ne, ns)


if __name__ == "__main__":
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    timeout = float(sys.argv[2]) if len(sys.argv) > 2 else 600
    here = os.path.dirname(os.path.abspath(__file__))
    data = json.load(open(os.path.join(here, "data", "structures_9_points_21_triples.json")))["structures"]
    todo = [(i, s) for i, s in enumerate(data) if "relabelling_to_model" not in s]
    assert len(data) == 22 and len(todo) == 21
    ok = 0
    with Pool(procs) as pool:
        for (i, _), (out, info) in zip(todo, pool.imap(run, [(9, [tuple(t) for t in s["triples"]], timeout) for _, s in todo])):
            closed = out.startswith("EMPTY")
            ok += closed
            print("structure %2d: %s -> %s" % (i, info, "1 in the ideal after %s saturation(s)" % out.split()[1] if closed else out), flush=True)
    print("%d / 21 structures excluded" % ok)
    print("UPPER BOUND STEP 2 VERIFIED" if ok == 21 else "FAILED")
