# The exact criterion (Ladkani 2.3(c)) rejects steps on non-monomial parents too, but only where the gate already refuses; the ALARM step is still the only gate-admitted rejection

author: scholar · round: 002 · kind: result (revision of round 001)
thread: T5 · bears on: H-015, F-016, F-038, E-032, E-055

## Response to referee

1. Negative control (required 1). Done, Evidence table A. At every vertex with an arrow out, on every LNA of
   length n and its opposite: gate refuses -> `tiltingPlus` False in 42/42 (n=5), 168/168 (n=6), 660/660 (n=7);
   gate admits -> True in 70/70, 252/252, 924/924. (The referee's n=4,5,6 counts, 10/42/168 and 20/70/252,
   agree at n=5,6; n=4 not rerun.) So the test is not vacuous and equals the gate on the starts.
2. "Guard" (required 2). Accepted. In `scholar_h015.py`, **guard = the Coxeter key of the reduced child equals
   the parent's** (`search._coxeterKeyOrNone`). It is not `procedure.isMutable` and not a search-time object.
   **Gate = `mutation.mutationIsPossibleAtVertex`** (delegates to `procedure.isMutable`). Round 001's
   "the guard never fires" therefore meant: every gate-admitted child kept its key. Skipped parents were those
   with a directed cycle in the quiver or `baseKey None` (`continue`, uncounted). Now counted
   (`scholar_nonmono.py`): **0 skipped** in every walk below; in round 001's own n=6/7 walks I did not count,
   but the same generator (distinct algebras from the LNAs) reaches none here, so I expect 0, unverified.
3. Second non-monomial negative (required 3). Found, with a caveat that changes what it shows: see Claim.
   There is NO second non-monomial step that the gate admits and 2.3(c) rejects. I say so plainly.
4. Not cited as support for H-015: agreed; stated again below.
5. n=7 depth-4 counts: not re-verified by anyone; the script is the same, the command is in Reproduction.

## Claim

On non-monomial parents (some relation has >= 2 terms, in practice a commutativity square) reached from
the LNAs of length 5, 6, 7 by unguarded BFS over distinct algebras, 2.3(c) rejects exactly the vertices
the gate refuses: n=6 depth 3: 214 rejected, all gate-refused, 714 accepted, all gate-admitted and all
Cartan-congruent (Prop 3.6); n=7 depth 2: 432 and 1008 likewise. The rejected parents look like
`b1 b2 = b3 b4` (commutative square) with the vertex at the source of the square's short side. So
`tiltingPlus` is now shown to discriminate on non-monomial algebras (a real second negative family), but it
adds nothing beyond the gate there: the only known step where the gate admits and 2.3(c) rejects remains
E-032 ALARM step 7. It does NOT show the guard is sufficient, nor that gate = 2.3(c) in general (n >= 10 is
where they differ, F-038), nor that the coding is right on cases outside these small ranges.

## Evidence

Table A, negative control (starts; every vertex with an arrow out; gate x tilt):

| n | gate refuses, tilt False | gate admits, tilt True | other |
|---|---|---|---|
| 5 | 42 | 70 | 0 |
| 6 | 168 | 252 | 0 |
| 7 | 660 | 924 | 0 |

Table B, non-monomial parents (BFS, all vertices with an arrow out, not only gate-admitted):

| n, depth | algebras expanded | non-monomial parents | gate, tilt, cong | refused, NOTtilt | gate, NOTtilt | skipped |
|---|---|---|---|---|---|---|
| 5, 2 | 98 | 12 | 40 | 8 | 0 | 0 |
| 6, 3 | 910 | 194 | 714 | 214 | 0 | 0 |
| 7, 2 | 1,188 | 240 | 1,008 | 432 | 0 | 0 |

Every gate-admitted step on a non-monomial parent had `child Cartan = r C r^T` (the `cong` column), so the
two independent computations again agree. Illegal-relation children: none.

Why the ALARM step is still the one that counts: it is the only recorded gate-admits/2.3(c)-rejects step
(E-055), and it lives at n=12 (steps of `03033030` relation dual), far from this range. The negatives here are
"gate correctly refuses", the weakest kind of negative for a would-be replacement gate: they show 2.3(c) is
not trivially True, not that it is stronger than the gate. To get a second gate-admitted rejection one must
go where F-038 says corruption is possible (n >= 10, seven clean steps deep): that is the Menu 4
overnight audit, not a round.

## Reproduction

```
N=5 DEPTH=2 timeout 10m .venv/bin/python workshop/rounds/002/scholar_nonmono.py   # 1 s
N=6 DEPTH=3 timeout 10m .venv/bin/python workshop/rounds/002/scholar_nonmono.py   # 22 s
N=7 DEPTH=2 timeout 10m .venv/bin/python workshop/rounds/002/scholar_nonmono.py   # 46 s
```
Output: Part A counts, per-depth progress, skipped parents, Table B rows, and the first non-monomial
NOTtilt parents (relations printed as arrow paths; up to 200 kept). It reuses `tiltingPlus`, `cartan`,
`rplus` from `workshop/rounds/001/scholar_h015.py` (which must stay in the repository).

## Prior record

As in round 001: 2.3(c) unused in `research/` before E-055; F-016, F-038, F-047 as cited there. Grepped again
(`tiltingPlus`, `isTilting`, `2.3(c)`): nothing new since E-055. Nothing here touches RETRACTIONS.

## Code changed

New `workshop/rounds/002/scholar_nonmono.py` only. No library code, no tests needed or run.

## Next

- Decision for the chair on STEERING q3: the requirement (a second non-monomial *negative* example) is met in
  the weak sense (refused parents, 200+ recorded) and not met in the strong sense (gate-admitted). I suggest
  promoting `isTilting` anyway as a *cross-check* beside the gate (it is exact and cheap: 22 s for 194
  non-monomial parents), with a test on one commutative-square parent (from the first example line of the
  n=5 output, vertex 2) pinning False, and on ALARM step 7. The chair should decide; I do not treat q3 as
  satisfied.
- Overnight (Menu 4 already lists it): add `--nonmono-only` output there; the decisive row is gate=True,
  tilt=False with the guard passing.
- toolsmith: nothing beyond the promotion.
