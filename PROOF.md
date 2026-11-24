# Proof that a(9) = 20

**a(n)** is the largest number of distinct circles of radius 1 such that each passes through at
least 3 points of some set of n points in the plane (OEIS A003829, Erdős problem #104). Known:
a(3..8) = 1, 4, 4, 8, 12, 16 (Harborth and Mengersen 1986 for n ≤ 7, Harborth 1985 for n = 8).

**Theorem.** a(9) = 20. Moreover 24 ≤ a(10) ≤ 26.

The lower bounds are explicit configurations (section 5). The upper bound a(9) ≤ 20 has four
steps:

1. a realisation gives an *incidence structure* (section 1);
2. an exact linear rule (the *rhombus rule*) leaves 22 structures with 9 points and 21 triples, and
   none with 22 or more (section 2);
3. 21 of them have no realisation, by exact algebra (section 2.3);
4. the last one, **S\***, has no realisation (section 3): it reduces to 12 unit vectors subject to
   linear relations, a lemma on two triples of unit vectors with the same sum (section 3.3) turns
   the relations into conditions on the circular order of the vectors, and an enumeration of all
   105,660 possible orders of the first 9 vectors shows that the conditions cannot be met.

Steps 2, 3 and 4 use computer enumerations; section 6 lists what each program checks.

Notation: points of the plane are complex numbers; a *unit vector* is a complex number of modulus 1.

## 1. Incidence structures

Take n points and m distinct unit circles, each through at least 3 of the points. Keep 3 of the
points on each circle. This gives a *structure* (n, T): T is a list of m triples of points, and:

* the triples are distinct (three points lie on at most one circle);
* each pair of points lies in at most 2 triples (two distinct unit circles meet in at most two
  points, and through two points there are at most two unit circles);
* the m centres are pairwise distinct (a unit circle is determined by its centre).

A structure is *realisable* if there are pairwise distinct points and pairwise distinct centres
such that every point of every triple is at distance 1 from the centre of that triple. Removing a
point (and its triples), or a triple, from a realisable structure leaves a realisable structure.

So a(9) ≥ 21 would give a realisable structure with 9 points and 21 triples.

## 2. The rhombus filter and the 22 structures

### 2.1 The rhombus rule

Consider the symbols "points" and "centres", and the graph H whose edges are pairs of symbols known
to be at distance 1 (at first, the incidences point–centre). If x1 x2 x3 x4 is a 4-cycle of H with
x1 ≠ x3 and x2 ≠ x4 (as points of the plane), then x1 and x3 are the two intersection points of
the unit circles centred at x2 and x4, and

    x1 + x3 = x2 + x4.

This is a linear relation with integer coefficients. The filter (`rhombus.py`, and independently
`rhombus.cpp`) works in exact rational arithmetic in the space of formal combinations of the
symbols modulo all relations found so far, and repeats until nothing changes:

* *coincidence:* if two points, or two centres, become equal modulo the relations, contradiction;
  if an edge of H becomes 0, contradiction; if two edges of H become proportional with a ratio
  other than ±1, contradiction (both have length 1);
* *closure:* if x − y equals ± an edge vector of H, then |x − y| = 1 and xy becomes an edge;
* *distinctness* (used to apply the rhombus rule): two symbols are distinct if they have the same
  type, or if their difference is a non-zero rational multiple of an edge vector.

Every rule is a consequence of the geometry, so a structure that leads to a contradiction is not
realisable.

### 2.2 Generation (`generate.py`)

In a structure with k points and t triples, a point of minimum degree d has d ≤ ⌊3t/k⌋; removing it
leaves a realisable structure with k − 1 points and t − d triples. With the thresholds
tmin(5..9) = 3, 6, 9, 14, 20 one checks that t − ⌊3t/k⌋ ≥ tmin(k − 1) whenever t ≥ tmin(k). Hence
every realisable structure with k points and ≥ tmin(k) triples arises from a filtered structure
with k − 1 points and ≥ tmin(k − 1) triples by adding a point of minimum degree, whose triples form
a graph of maximum degree 2 on the old points (pair multiplicity ≤ 2). Starting from all
structures with 5 points, adding points this way, removing isomorphic copies (nauty canonical
forms) and applying the filter at each level gives:

| k | candidates (by number of triples) | pass the filter |
|---|---|---|
| 5 | 3: 3, 4: 4, 5: 2, 6: 1 | 3: 3, 4: 4 |
| 6 | 6: 24, 7: 11, 8: 3 | 6: 20, 7: 6, 8: 2 |
| 7 | 9: 542, 10: 239, 11: 74, 12: 9, 13: 1, 14: 1 | 9: 310, 10: 111, 11: 21, 12: 5 |
| 8 | 14: 17,814, 15: 4,781, 16: 389, 17: 30, 18: 1 | 14: 3,944, 15: 351, 16: 22, 17: 1 |
| 9 | 20: 350,655, 21: 11,761, 22: 615, 23: 4, 24: 2 | 20: 1,675, **21: 22** |

So **no structure with 9 points and ≥ 22 triples passes**, and 22 structures with 21 triples do
(`data/structures_9_points_21_triples.json`). (The filter alone already gives back the known values
1, 4, 4, 8, 12 for n = 3..7; for n = 8 one structure with 17 triples passes the filter, and it is
excluded by the algebra of section 2.3 — it is the last case of Harborth's proof.)

### 2.3 Twenty-one of the 22 structures (`exclude_21.py`)

Modulo the relations of the filter, each symbol is a rational combination of a few free symbols;
fix one of them at the origin (translation). The edge vectors of H become linear forms e(z) in
complex unknowns z_i, rotate so that one of them equals 1, and introduce unknowns w_i standing for
the conjugates. Each unit edge gives the polynomial equation e(z)·e(w) = 1. Each difference L(z)
of two points or two centres must be non-zero; saturating the ideal by these linear forms one at
a time (Singular), 1 enters the ideal for all 21 structures. Hence they have no solution, even
over ℂ with w unrelated to z, with distinct points and distinct centres.

(In the original search these 21 structures were excluded by further exact linear rules — the
triangle inequality on minimal rational relations between 3 or 4 unit vectors, and the fact that
four unit vectors with zero sum cancel in pairs, with branching — and no branch remained. The
Gröbner computation is an independent second proof.)

The 22nd structure, S\*, has complex solutions forming a curve, so no such certificate exists; a
real argument is needed.

## 3. The structure S\*

### 3.1 Description

Nine points in three groups X = {x0, x1, x2}, A = {a0, a1, a2}, B = {b0, b1, b2}, indices in ℤ/3.
The 21 circles: the three *group circles* (through X, through A, through B) and the 18
*transversal circles* (x_i, a_j, b_k) with i + j + k ≢ 0 (mod 3). Every point has degree 7. Every
pair (x_i, a_j), (x_i, b_k), (a_j, b_k) lies on exactly two transversal circles (the two values of
the third index that avoid the forbidden sum), which gives 27 rhombus identities
p + q = M + M′. Removing one group circle gives a structure with 20 triples that *is* realisable:
it is the structure of our 20-circle configuration family.

### 3.2 Reduction to twelve unit vectors

Write M(i, j) for the centre of the transversal circle through x_i, a_j and b_k with
i + j + k ≡ 1, and put

    E[s][l] = x_l − M(l, l − s)        (s, l ∈ ℤ/3),

nine unit vectors. For a group circle with centre O and points q0, q1, q2, put ξ_l = q_l − O
(three distinct unit vectors). Using only the 27 rhombus identities (exact linear algebra over ℚ,
checked by part 1 of `verify_sstar.py`):

* **(R1)** each of the 54 vectors point − centre of the transversal circles equals ±E[s][l], each
  E[s][l] six times;
* **(R2)** if s ≠ s′, E[s][l] − E[s′][l′] equals ± the difference of two points or of two centres,
  so **E[s][l] ≠ E[s′][l′]** (no conclusion within a row);
* **(R3)** the *row relations*: E[s][0] + E[s][1] + E[s][2] is the same for s = 0, 1, 2 — the three
  triangles E[s][·] inscribed in the unit circle have the same centroid c;
* **(R4)** the *group relations*: for the group circle of X (j = 0), of A (j = 1), of B (j = 2),
  with a suitable numbering of its points, for all l, r ∈ ℤ/3,

        ξ_l − ξ_{l−1} = Σ_{s≠r} E[s][l + js] − Σ_{m≠l+jr} E[r][m].

  By (R3) the right-hand side equals Σ_s (E[s][l + js] − c): the sides of the group triangle are
  determined by the E's.

Each row relation and each group relation says that a signed sum of **six unit vectors** vanishes;
splitting the six terms into two groups of three (10 ways), two triples of unit vectors have the
same sum. The lemma below turns this into a condition on their circular order.

### 3.3 The lemma: two unit triples with the same sum

**Lemma.** Let T = {a, b, c} and T′ = {a′, b′, c′} be multisets of unit complex numbers with
a + b + c = a′ + b′ + c′ = H and T ≠ T′. Then |H| < 3 and:

1. if |H| = 1: H ∈ T and H ∈ T′, and the two other elements of T are opposite, as are those of T′;
2. if |H| < 1: the six numbers are distinct, each triangle T, T′ has its three arcs < π (it is
   *acute*: the centre is inside), and the points of T and T′ alternate around the circle;
3. if 1 < |H| < 3: each triple lies in an open half-circle, and one can name the triples A, B so
   that, going counter-clockwise once around the circle,

        α1 < β1 ≤ β2 < α2 ≤ α3 < β3,

   where α1, α2, α3 are the points of A and β1, β2, β3 those of B, each triple listed from the start
   of its short arc (the arc of length < π that contains it); equality is possible only as β1 = β2
   or α2 = α3 (a double point). We call this the *mixed order* (pattern A B B A A B).

The program uses the following consequence, which only involves the circular order, equalities,
antipodes, and whether arcs are shorter than, equal to or longer than π: *if T and T′ have the same
sum, then T = T′, or (1), or (2), or (3) holds* (with (3) checked in either orientation, which is
harmless: it is a weaker condition).

**Step A: an identity.** Let x, y, z ≥ 0 be the three arcs between consecutive points of T
(x + y + z = 2π; a double point gives an arc 0). Then

    |a + b + c|² = 1 − 8 cos(x/2) cos(y/2) cos(z/2).

*Proof.* Rotate so that a = 1, b = e^{ix}, c = e^{i(x+y)}. Then
|a + b + c|² = 3 + 2(cos x + cos y + cos z). Since (x + y)/2 = π − z/2,
cos x + cos y = 2 cos((x+y)/2) cos((x−y)/2) = −2 cos(z/2) cos((x−y)/2), and
cos z = 2cos²(z/2) − 1 = −2 cos(z/2) cos((x+y)/2) − 1. Adding,
cos x + cos y + cos z = −1 − 2cos(z/2)[cos((x−y)/2) + cos((x+y)/2)]
= −1 − 4 cos(x/2) cos(y/2) cos(z/2). ∎

Since x/2 + y/2 + z/2 = π, at most one half-arc exceeds π/2. Hence: |H| < 1 iff all three arcs
are < π (acute); |H| = 1 iff one arc equals π (two opposite points); |H| > 1 iff one arc is > π,
i.e. the triple lies in an open half-circle. Also |H| ≤ 3 with equality iff a = b = c, so T ≠ T′
forces |H| < 3.

**Step B: the cubic.** Let p = abc. Since 1/a = ā etc., ab + bc + ca = abc(ā + b̄ + c̄) = p·H̄.
So T is the multiset of roots of

    f_p(z) = z³ − H z² + p H̄ z − p,

and T is determined by (H, p). As T ≠ T′, their products p ≠ p′.

**Step C: |H| = 1.** Then a + b + c + (−H) = 0 with four unit vectors. Four unit vectors
u1 + u2 + u3 + u4 = 0 cancel in pairs: put s = u1 + u2 = −(u3 + u4); if s = 0 then u2 = −u1 and
u4 = −u3; if s ≠ 0, the two unit vectors with sum s are the two intersection points of the unit
circle with the circle |z − s| = 1, so {u1, u2} = {−u3, −u4}. If −H cancels with a, then a = H and
b = −c; if a = −b, then c = H; and so on. So T = {H, u, −u}, and likewise T′ = {H, u′, −u′}. This is
case 1.

**Step D: the map g (|H| ≠ 1).** For |z| = 1, |1 − H̄z| = |z̄ − H̄| = |z − H| ≠ 0. Put

    g(z) = z² (z − H) / (1 − H̄ z).

Then f_p(z) = (1 − H̄z)(g(z) − p), and |g(z)| = 1 on the unit circle. So the triples with sum H
are exactly the fibres g⁻¹(p) on the unit circle (with multiplicity) that have three points
counted with multiplicity; a double root of f_p is a zero of g′. Write g(e^{iθ}) = e^{iG(θ)} with
G real and smooth. For |α| < 1 and the Blaschke factor B_α(z) = (z − α)/(1 − ᾱz), the derivative of
arg B_α(e^{iθ}) is the Poisson kernel P_α(θ) = (1 − |α|²)/|e^{iθ} − α|² > 0, and ∫ P_α = 2π.

*Case |H| < 1.* g = z²·B_H, so G′ = 2 + P_H > 0: G is strictly increasing and gains 6π per turn.
Each value is taken exactly three times per turn, at simple points: T has three distinct points.
Between two consecutive points of the fibre over p, G increases by exactly 2π, so it takes the
value arg p′ (≠ arg p mod 2π) exactly once in between. The fibres alternate: with step A, this is
case 2.

*Case 1 < |H| < 3.* Put b = 1/H̄ (|b| < 1). A direct computation gives
(z − H)/(1 − H̄z) = (b/b̄)/B_b(z), so

    G′(θ) = 2 − P_b(θ).

P_b is largest in the direction of H, where it equals (|H| + 1)/(|H| − 1) > 2 (because |H| < 3),
smallest in the opposite direction, where it equals (|H| − 1)/(|H| + 1) < 2, and strictly
monotone on each half-circle in between. So G′ < 0 on an open arc D = (c1, c2) around the
direction of H, G′ > 0 on the complementary open arc I, G′ = 0 only at c1 and c2, where G″ ≠ 0
(the zeros of g − p there are exactly double). Over one turn G gains 2π (= 4π − 2π).

Lift angles to ℝ with c1 < c2 < c1 + 2π, and normalise G(c1) = 0. Let δ = G(c1) − G(c2) > 0;
then G decreases from 0 to −δ on [c1, c2] and increases from −δ to 2π on [c2, c1 + 2π].
We have δ < 2π: otherwise a value v ∈ (2π − δ, 0) would be taken twice on D (as v and v − 2π)
and twice on I, four points in one fibre of a cubic. Hence:

* for v ∈ (−δ, 0), the fibre over e^{iv} consists of three distinct simple points: d_v ∈ D with
  G(d_v) = v, and x_v, y_v ∈ I with G(x_v) = v, G(y_v) = v + 2π;
* for v = 0 it is {c1 (double), x_0}; for v = −δ it is {c2 (double), y_{−δ}};
* any other value is taken once (a simple root): f_p then has a single root on the circle, and
  there is no triple.

So the triples with sum H are T_v = {y_v, d_v, x_v}, v ∈ [−δ, 0], with lifts
y_v ≤ c1 ≤ d_v ≤ c2 ≤ x_v < y_v + 2π (the last because G(x_v) = v < v + 2π = G(y_v + 2π) and G
increases on I). As v decreases, d_v moves forward and x_v, y_v move backward. So for
−δ ≤ v′ < v ≤ 0:

    y_{v′} < y_v ≤ d_v < d_{v′} ≤ x_{v′} < x_v < y_{v′} + 2π,

the last inequality because G(x_v) = v < v′ + 2π (as v − v′ ≤ δ < 2π). Equalities: y_v = d_v only
for v = 0 and d_{v′} = x_{v′} only for v′ = −δ.

It remains to see that the short arc of T_v runs from y_v to x_v, i.e. that the arc from x_v to
y_v + 2π (the gap g3) is the one longer than π. By step A each T_v has exactly one gap > π. The
three gaps g1 = d_v − y_v, g2 = x_v − d_v, g3 = y_v + 2π − x_v are continuous in v (inverse
functions of strictly monotone ones), so the sets {v : g_k(v) > π} are open, disjoint, and cover
the interval [−δ, 0]; by connectedness one of them is everything. At v = 0, g1 = 0; at v = −δ,
g2 = 0. So g3 > π for every v. With A = T_{v′} and B = T_v this is exactly the mixed order of
case 3. (This is the place where one must not assume that the two "outer" points of a triple are
the ends of its short arc without proof: the argument above shows it, including when a double root
sits at c1 or c2.) ∎

### 3.4 Weak order types and the enumeration (part 3 of `verify_sstar.py`)

Suppose S\* is realised. Take the 9 vectors E[s][l] and the three ξ_l of one group circle, and
their opposites: 24 points of the unit circle. Their *weak order type* is the list of the distinct
directions met in a half-turn starting at E[0][0], each with the signed vectors (v or −v) found
there; a signed vector gets a position in ℤ/2m (m = number of slots), the opposite of position q
being q + m. Two positions q < q′ span an arc shorter than π iff q′ − q < m, equal to π iff
q′ − q = m. So the order type determines everything the lemma uses.

The order type must satisfy:

* E[s][l] ≠ E[s′][l′] for s ≠ s′ (R2), and ξ_0, ξ_1, ξ_2 pairwise distinct;
* for each of the 3 row relations (R3) and each of the 10 splits of its six terms into two triples,
  the lemma;
* for each of the 9 group relations (R4) of that group and each split, the lemma.

The program builds all weak order types by inserting the vectors one at a time — into an existing
slot, or into a new slot after any existing slot, with either sign; this produces every weak order
type (it was also compared with an exhaustive enumeration of positions on a grid for up to 5
vectors) — and discards a partial type as soon as a relation whose vectors are all placed violates
the lemma (the restriction of a valid type to fewer vectors is valid).

**Result.** There are exactly **105,660** weak order types of the 9 vectors E compatible with
(R2) and the row relations. For **every** one of them, some group j admits no placement of
ξ_0, ξ_1, ξ_2: the first such j is j = 0 for 88,860 types, j = 1 for 14,064 and j = 2 for 2,736
(`data/sstar_order_types.jsonl.gz`). Since a realisation of S\* would give a valid placement for
all three groups at once, **S\* is not realisable**, and a(9) ≤ 20.

Two remarks. Removing one group circle leaves a realisable structure; the order types of its
realisations are found in the list, with both present groups accepted and the missing group
excluded, as it must be. And the lemma conditions are only necessary: the enumeration never uses
them to *construct* anything.

## 4. Ten points

The same chain (from the 9-point structures with ≥ 19 triples, 51,567 of them) gives 1,174
10-point structures with ≥ 26 triples passing the filter (1,172 with 26, 2 with 27, none with 28).
The two with 27 triples were proved unrealisable through an unrealisable substructure; of the
1,172 with 26 triples, 68 contain an already excluded substructure, 1,086 were proved unrealisable
directly by the exact rules and 2 through a substructure; 16 remain undecided
(`data/undecided_10_points_26_triples.json`), none of them realised numerically. Hence
a(10) ≤ 26. This part was **not** re-checked independently and is not covered by the verifiers
of this repository.

## 5. Lower bounds

`data/config_9_points_20_circles.json` and `data/config_10_points_24_circles.json` give each point
and centre as an integer combination Σ A[x][i]·w_i of unit complex numbers w_i = e^{iθ_i}, so that
every incidence vector is ±w_i or ±E_j for a few rational combinations E_j. The equations
|E_j|² = 1 are certified to have an exact solution in a tiny box of angles by the interval Krawczyk
operator, and in that box the points are pairwise distinct and the centres pairwise distinct
(`verify_configs.py`). The 20-circle configuration belongs to a two-parameter family (S\* minus one
group circle); the 24-circle one adds to it the orthocentre of one group.

## 6. What the computer checks

| claim | program | dependencies |
|---|---|---|
| configurations with 20 (9 points) and 24 (10 points) circles exist | `verify_configs.py` | mpmath |
| only 22 structures with 9 points and ≥ 21 triples pass the rhombus filter | `generate.py` with `rhombus.py` or `rhombus.cpp` | pynauty |
| 21 of them have no realisation | `exclude_21.py` | Singular |
| reduction of S\* to (R1)–(R4) | `verify_sstar.py`, part 1 | standard library |
| the lemma | hand proof above (part 2 of `verify_sstar.py` is only a numerical sanity check) | — |
| no weak order type of the 12 vectors satisfies the lemma | `verify_sstar.py`, part 3 | standard library |

In our search, each of these steps was also done a second time by independently written code: the
lower bounds with our own Krawczyk test at 120 digits on three different realisations; the 22
structures by two generation paths (removing a point of minimum degree, and of maximum degree) with
filters written in Python and in C++; the exclusion of the 21 by both the linear branching rules
and Gröbner bases; the reduction of S\* re-derived symbolically; the order-type enumeration
re-implemented with a different representation (positions on a grid) — same 105,660 types, same
excluded group for each.
