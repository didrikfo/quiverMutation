# Off the walks, "reject iff long square" holds in one direction and fails literally in the other: 14 hand-built gate-admitted rejections without E-099's square, none reachable by a guarded walk

author: skeptic · round: 023 · kind: negative
thread: T5 · bears on: H-015, E-099, E-102

## Claim

(1) "Long square => reject" survives. In 12 000 random hand-built parents (n = 5, 6, 7, mode A seeds a long-sided
square on purpose) there are 362 gate-admitted tilting steps at a vertex where the square's shape is present in the
relations, and in every one the long relation is redundant: the truncated combination `sum c p[:-1]` is already in the
ideal (another relation, e.g. `abd = acd`, kills it). No gate-admitted step has a *genuine* long square (truncation
nonzero mod I) and `tiltingPlus` True. This is expected: the truncation is the kernel element of the one map, so
"genuine long relation => reject" is nearly a restatement of Ladkani 2.3(c).
(2) "Reject => long square" is FALSE as `hasLongSquare` of rounds/022 states it. Gate-admitted rejecting steps with
`hasLongSquare` False: 12 of 4 463 (n = 6), 11 of 2 388 (n = 7), 1 and 2 (mode C, no seeding). Three kinds:
 a. v has TWO out-arrows, each carrying a long relation: `abde = acde` and `abdf = acdf` (n = 6, edges 12 13 24 25 34 45 46,
    v = 4). `abd - acd` is nonzero and killed by both arrows. 19 of the 26 examples (all runs) are this kind (18 shaped on both
    arrows; one has the long relation on one arrow and a monomial `2->4->6 = 0`, `3->4->6 = 0` on the other).
 b. shared longer suffix: `abdef = acdef` with v = e (n = 6, edges 12 13 14 24 34 36 45 46 56, v = 5): the two paths have the
    same third-last vertex, so the "distinct third-last" test fails, but the truncations abde, acde differ.
 c. long relation made short by reduction: a 3-term sum relation one of whose terms lies in the ideal (monomial `123`);
    after reducing it is a 2-term long square ending in the single out-arrow (n = 6, 7 examples in the output files).
(3) None of the off-walk kinds can be reached by the E-086/E-102 guarded walk: the Coxeter polynomials of examples a and b
    are `(1,1,-2,-4,-2,1,1)` and `(1,-2,-8,-11,-8,-2,1)`, in no LNA / relation-dual class at n = 6, so the key-preserving
    walk cannot contain them. (Not shown for c; not tested at n = 7.)
It does NOT claim: that the E-102 numbers are wrong (they are a statement about walk-reached steps, unaffected), nor that kinds
a-c are counterexamples to a theorem "kernel of the step-7 map comes from a minimal relation": under the right
reading (see Next) they are all instances of it.

## Evidence

Script `skeptic_offwalk.py n samples seed mode`: random acyclic quiver on 1..n (arrow i->j only for j-i <= 3, no parallel
arrows), 0-3 random homogeneous relations (monomial, or sum of 2-3 distinct paths of equal length >= 3 arrows), mode A seeds
one long square first, mode C none. At each vertex with an out-arrow: gate = `mutationIsPossibleAtVertex`, `tiltingPlus`
(copied from rounds/001), `hasLongSquare` on `alg.rels` (as rounds/022), a dict-form shape test on `procedure.relationsFrom`
("shaped"), and "genuine" = shape and truncated combination nonzero mod `idealBasis`. Counts of gate-admitted steps:

| run | parents | gate, rej, genuine long | gate, tilt, shaped (all non-genuine) | gate, tilt, no shape | gate, rej, no shape (B) |
|---|---|---|---|---|---|
| n=5 A seed 2 | 3000 | 2091 | 127 | 8816 | 0 |
| n=6 A seed 3 | 6000 | 4026 | 149 | 22397 | 12 |
| n=7 A seed 4 | 4000 | 2377 | 83 | 18277 | 11 |
| n=6 C seed 5 | 6000 | 7 | 2 | 20079 | 1 |
| n=7 C seed 6 | 4000 | 6 | 1 | 16450 | 2 |

Gate-refused steps behave the same (genuine => rejecting; refused with no shape rejecting: 419, 343, 730, 493, 184).
B by (out-arrows at v, shaped on each arrow): n=6 A: (2,TT) 10, (2,TF) 1, (1,F) 1; n=7 A: (2,TT) 8, (1,F) 3; mode C: (1,F) 1 + 2.
Walk-reachability: `skeptic_reach.py 6 400` rebuilds the guarded BFS (class of the parent's Coxeter key). Examples a and b are in no
class ("NOT A CLASS OF AN LNA/dual"). Control, the n = 6 padded E-080 (edges 12 13 24 34 45 56, `1245 = 1345`): class 0, FOUND among
73 020 algebras (stopped at 400 s, 46 583 expansions). That contradicts nothing recorded in E-086 ("E-080's own algebra is not
reached", n = 5) but note the n = 6 padding is reached.
Sampling is not uniform and not exhaustive: no claim about frequency. The generator only draws arrows with j-i <= 3 and
relations of >= 3 arrows (none of length 2) -- so length-2 sums and parallel arrows are not covered.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/023/skeptic_offwalk.py 6 6000 3 A   # 19 s; also 5 3000 2 A, 7 4000 4 A, 6 6000 5 C, 7 4000 6 C
timeout 10m .venv/bin/python workshop/rounds/023/skeptic_reach.py 6 400            # ~7 min
```
Outputs (with the example parents listed, 400-row cap): `skeptic_offwalk_n{5,6,7,6,7}_s{2,3,4,5,6}{A,C}.txt`.

## Prior record

E-102 (limits) explicitly left "a long-square tilting step off the walks (E-080 family)" and the `alg.rels` vs
`relationsFrom` question open; STATE T5 requests it of the skeptic. E-080 is the 1-out case; the 2-out and shared-suffix cases
are not in `research/` (grep "two out", "shared suffix": none). Not in RETRACTIONS.md. The redundant-long-relation
artefact is the presentation question of E-102's referee point 3 answered the other way: `alg.rels` can carry a long relation
that is not minimal, and then `hasLongSquare` is True at a tilting step (362 of the 'tilt' steps above), so the test must be run
on a minimal presentation (or with the truncation-nonzero check) before it is cited off the walks.

## Code changed

None in the library. New: `workshop/rounds/023/skeptic_offwalk.py`, `skeptic_reach.py` and output txt files.

## Next

- theorist: state the theorem the data support: for gate-admitted v, `tiltingPlus` fails iff after reduction modulo I there is
  x in e_aAe_v, x != 0, with x beta in I for every out-arrow beta; derive that such x comes from a (minimal) relation
  ending in each out-arrow. Kinds a-c are then the multi-out / non-injective-truncation cases, not exceptions.
- experimentalist: E-102's "no long square on tilting steps" should be re-read with "genuine" (truncation nonzero mod I);
  rerun `hasLongSquare` on a minimal generating set of the relations at n = 6 if cheap. Count how often the walks' parents have
  a vertex with two out-arrows and relations on both (expected 0 by the key argument).
- skeptic (later): exhaustive (not random) enumeration of acyclic quivers with <= 7 arrows at n = 6 for kinds a-c; n = 8.
