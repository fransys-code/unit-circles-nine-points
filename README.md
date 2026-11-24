# Nine points carry at most 20 unit circles

Let **a(n)** be the largest number of distinct circles of radius 1 that can each pass through at
least 3 points of a set of n points in the plane ([OEIS A003829](https://oeis.org/A003829);
[Erdős problem #104](https://www.erdosproblems.com/104) asks for its growth). The known values are
a(3..8) = 1, 4, 4, 8, 12, 16 (Harborth and Mengersen 1986; Harborth 1985 for n = 8).

| n | a(n) | status |
|---|---|---|
| 9 | **= 20** | lower bound: certified configuration; upper bound: computer-assisted proof ([PROOF.md](PROOF.md)) |
| 10 | **24 ≤ a(10) ≤ 26** | lower bound: certified configuration; upper bound: search, 16 structures with 26 circles undecided |

The proof of a(9) ≤ 20 in one paragraph: keeping 3 points on each circle gives an incidence
structure; an exact linear rule (if two unit circles centred at x2 and x4 both pass through x1 and
x3, then x1 + x3 = x2 + x4) leaves 22 structures with 9 points and 21 circles, and none with more;
exact algebra (Gröbner bases) excludes 21 of them. In the last one, S\*, nine points in three groups
of three, the relations reduce to twelve unit vectors; a lemma on two triples of unit vectors with
the same sum (acute and interleaved if the sum has modulus < 1, in a "mixed" order if it is between
1 and 3) turns them into conditions on the circular order of the vectors, and an enumeration of the
105,660 possible orders shows that one of the three group circles always fails.

> **Status.** The search, the proof and the verification programs were written with the help of an
> AI assistant. Every computational step was redone by independently written code (second
> generation path, filters in two languages, a second exclusion method, a re-implemented
> order-type enumeration) with identical results, and the hand proof of the lemma was checked by a
> second AI agent and re-proved independently. **No human referee has read the proof, and it is not
> formalised.** The two implementations share two external tools: pynauty (isomorphism rejection)
> and Singular (Gröbner bases). To our knowledge the values a(9) = 20 and the bounds for a(10) are
> new: OEIS A003829 stops at a(8) (checked on 24 November 2026), and we found no other source.
> Corrections and pointers to prior work are welcome: please open an issue.

*Résumé en français : a(n) est le nombre maximal de cercles de rayon 1 passant chacun par au moins
trois points d'un ensemble de n points du plan (OEIS A003829, problème d'Erdős n° 104 ; valeurs
connues jusqu'à n = 8). Nous démontrons a(9) = 20 et 24 ≤ a(10) ≤ 26. La configuration de 9 points
et 20 cercles est certifiée (arithmétique d'intervalles). La majoration est une preuve assistée par
ordinateur : règle linéaire du losange, 22 structures restantes, 21 exclues par bases de Gröbner, la
dernière (S\*) par un lemme sur deux triplets de nombres complexes de module 1 de même somme et
l'énumération de 105 660 ordres circulaires. Tous les calculs ont été refaits par un code écrit
indépendamment ; le lemme est démontré à la main et relu par un second agent d'IA, mais ni relu par
un mathématicien ni formalisé. Nouveauté vérifiée seulement par recherche bibliographique rapide.*

Figures and explanations (French / English): https://totocaedres.fransys.io/pistes/cercles-unite/notre-travail

## Files

- `PROOF.md` — the proof (structures, rhombus rule, reduction of S\*, the lemma with its proof,
  the enumeration).
- `data/config_9_points_20_circles.json`, `data/config_10_points_24_circles.json` — the certified
  configurations: triples of points on each circle, the exact incidence data, the certified box of
  angles, and approximate coordinates of the points (`approximate_points`).
- `data/structures_9_points_21_triples.json` — the 22 incidence structures with 9 points and 21
  triples that pass the rhombus filter, with how each is excluded; for S\*, the relabelling onto the
  model of PROOF.md.
- `data/sstar_order_types.jsonl.gz` — the 105,660 weak order types of the nine vectors E of S\*, each
  with the group circle it excludes (SHA-256
  `f19cdc18a71b537a569c1462d57bddea230c403be49375d8f6b56d0626df80a4`).
- `data/undecided_10_points_26_triples.json` — the 16 structures with 10 points and 26 triples left
  undecided.
- `verify_configs.py` — lower bounds (needs `mpmath`).
- `rhombus.py`, `rhombus.cpp` — the exact rhombus filter, written twice (Python standard library;
  C++, optional and faster).
- `generate.py` — regenerates the structures and checks that exactly the 22 of the data file have
  21 triples and none has more (needs `pynauty`).
- `exclude_21.py` — excludes the 21 structures other than S\* by Gröbner bases (needs `Singular`).
- `verify_sstar.py` — S\*: the reduction (exact linear algebra), a numerical sanity check of the
  lemma, and the enumeration of the order types (Python standard library only).

## Verification

```bash
pip install mpmath pynauty            # Singular: e.g. apt install singular
python3 verify_configs.py
g++ -O2 -o rhombus rhombus.cpp
python3 generate.py ./rhombus 4       # or: python3 generate.py - 4  (pure Python filter, slower)
python3 exclude_21.py 4
python3 verify_sstar.py 8             # python3 verify_sstar.py --quick : parts 1 and 2 only
```

Outputs obtained (Python 3.12, mpmath 1.3.0, pynauty 2.8.8.1, Singular 4.3.2):

```
$ python3 verify_configs.py
config_10_points_24_circles.json: certified, 10 distinct points, 24 distinct unit circles through 3 of them (separation about 0.122 between points, 0.041 between centres)
config_9_points_20_circles.json: certified, 9 distinct points, 20 distinct unit circles through 3 of them (separation about 0.132 between points, 0.083 between centres)
ALL CONFIGURATIONS CERTIFIED
```

```
$ python3 generate.py ./rhombus 4        # about 4 minutes
5 candidates {3: 3, 4: 4, 5: 2, 6: 1} pass the filter {3: 3, 4: 4}
6 candidates {6: 24, 7: 11, 8: 3} pass the filter {6: 20, 7: 6, 8: 2}
7 candidates {9: 542, 10: 239, 11: 74, 12: 9, 13: 1, 14: 1} pass the filter {9: 310, 10: 111, 11: 21, 12: 5}
8 candidates {14: 17814, 15: 4781, 16: 389, 17: 30, 18: 1} pass the filter {14: 3944, 15: 351, 16: 22, 17: 1}
9 candidates {20: 350655, 21: 11761, 22: 615, 23: 4, 24: 2} pass the filter {20: 1675, 21: 22}
9 points: 0 structures with >= 22 triples, the 22 with 21 triples are those of the data file
UPPER BOUND STEP 1 VERIFIED
```

```
$ python3 exclude_21.py 3                # about 20 seconds
structure  0: 3 complex unknowns, 11 equations, 77 saturations available -> 1 in the ideal after 0 saturation(s)
structure  2: 3 complex unknowns, 11 equations, 85 saturations available -> 1 in the ideal after 0 saturation(s)
structure  3: 6 complex unknowns, 13 equations, 123 saturations available -> 1 in the ideal after 22 saturation(s)
structure  4: 5 complex unknowns, 13 equations, 114 saturations available -> 1 in the ideal after 1 saturation(s)
structure  5: 4 complex unknowns, 12 equations, 99 saturations available -> 1 in the ideal after 0 saturation(s)
structure  6: 4 complex unknowns, 12 equations, 97 saturations available -> 1 in the ideal after 1 saturation(s)
structure  7: 3 complex unknowns, 11 equations, 83 saturations available -> 1 in the ideal after 1 saturation(s)
structure  8: 3 complex unknowns, 11 equations, 90 saturations available -> 1 in the ideal after 1 saturation(s)
structure  9: 4 complex unknowns, 12 equations, 102 saturations available -> 1 in the ideal after 1 saturation(s)
structure 10: 3 complex unknowns, 11 equations, 85 saturations available -> 1 in the ideal after 1 saturation(s)
structure 11: 3 complex unknowns, 10 equations, 81 saturations available -> 1 in the ideal after 0 saturation(s)
structure 12: 3 complex unknowns, 11 equations, 83 saturations available -> 1 in the ideal after 0 saturation(s)
structure 13: 3 complex unknowns, 11 equations, 88 saturations available -> 1 in the ideal after 0 saturation(s)
structure 14: 4 complex unknowns, 12 equations, 106 saturations available -> 1 in the ideal after 0 saturation(s)
structure 15: 3 complex unknowns, 11 equations, 82 saturations available -> 1 in the ideal after 0 saturation(s)
structure 16: 2 complex unknowns, 9 equations, 51 saturations available -> 1 in the ideal after 0 saturation(s)
structure 17: 4 complex unknowns, 12 equations, 116 saturations available -> 1 in the ideal after 1 saturation(s)
structure 18: 3 complex unknowns, 10 equations, 72 saturations available -> 1 in the ideal after 0 saturation(s)
structure 19: 3 complex unknowns, 11 equations, 81 saturations available -> 1 in the ideal after 0 saturation(s)
structure 20: 2 complex unknowns, 9 equations, 54 saturations available -> 1 in the ideal after 0 saturation(s)
structure 21: 6 complex unknowns, 14 equations, 123 saturations available -> 1 in the ideal after 9 saturation(s)
21 / 21 structures excluded
UPPER BOUND STEP 2 VERIFIED
```

```
$ python3 verify_sstar.py 8              # about 6 minutes on 8 cores
S* of the data file = model (points 0,1,2 = X, 3,4,5 = A, 6,7,8 = B; 3 group circles and the 18 transversals (x_i, a_j, b_k) with i + j + k != 0 mod 3)
rhombus identities p + q = M + M': 27  rank: 19
(1) the 54 point - centre vectors of the transversal circles are +-E[s][l], each E six times
(2) E[s][l] - E[s'][l'] = +-(two points or two centres) for s != s', so E[s][l] != E[s'][l']
(3) the three rows of E have the same sum
(4) group relations: circle X gives j = 0, circle A gives j = 1, circle B gives j = 2
PART 1 (reduction) VERIFIED
pairs of triples with equal sums, all accepted by the lemma: {'1 < |H| < 3': 5434, '|H| < 1': 2034, '|H| = 1': 6651}
PART 2 (sanity check) PASSED
weak order types of the 9 vectors E compatible with the row relations: 105660 (70 s)
same list as data/sstar_order_types.jsonl.gz
10000 {0: 7815, 1: 1614, 2: 571, None: 0} (103 s)
20000 {0: 16518, 1: 2830, 2: 652, None: 0} (129 s)
30000 {0: 25589, 1: 3721, 2: 690, None: 0} (152 s)
40000 {0: 34185, 1: 5089, 2: 726, None: 0} (176 s)
50000 {0: 42411, 1: 6450, 2: 1139, None: 0} (205 s)
60000 {0: 49644, 1: 8381, 2: 1975, None: 0} (239 s)
70000 {0: 58316, 1: 9511, 2: 2173, None: 0} (263 s)
80000 {0: 66879, 1: 10910, 2: 2211, None: 0} (285 s)
90000 {0: 75601, 1: 12129, 2: 2270, None: 0} (310 s)
100000 {0: 84625, 1: 13042, 2: 2333, None: 0} (332 s)
first excluded group: j=0: 88860, j=1: 14064, j=2: 2736, none: 0
agreement with the table: 105660 / 105660
PART 3 (enumeration) VERIFIED: S* has no realisation (349 s)
```

## How it was found

The search programs are not published here; nothing above depends on them. The structures were
first generated and filtered as described; the 20-circle configurations were found by numerical
solution of the structures that survive, then certified. For S\* we first observed numerically that
on the realisable family "S\* minus one group circle" the third group circle always has radius > 1,
and that S\* has complex solutions (so that no Gröbner certificate can exist); the order-type argument
came after.

## References

- H. Harborth, L. Mengersen, Discrete Mathematics 60 (1986): a(n) for n ≤ 7.
- H. Harborth, 3rd Colloquium on Discrete Geometry, Salzburg (May 1985), pp. 163–168: a(8) = 16.
- P. Erdős, G. Purdy, Extremal problems in combinatorial geometry, in Handbook of Combinatorics (1995).
- OEIS A003829, https://oeis.org/A003829 ; Erdős problem #104, https://www.erdosproblems.com/104

## Acknowledgements

The search, the proofs and the verification code were written with the help of Claude (Anthropic).

## Licence

Code: MIT (`LICENSE`). Data in `data/`: CC BY 4.0 (`LICENSE-DATA.md`).
