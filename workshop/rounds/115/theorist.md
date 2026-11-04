# Hochschild cohomology HH^0..2 says YES on all 13 of E-201's n = 7 class-1 steps and NO on a planted same-Cartan pair, but is blind to every LNA, including the n = 9 F-010 pair: a NO-capable key-free invariant exists only outside the LNAs (n = 7 class 1, steps and controls; no key-kept inequivalent pair reached)

author: theorist · round: 115 · kind: negative
thread: T10 · bears on: H-015, R-candidate "a key-guarded walk stays in one class", E-201, E-212, E-298, E-300
scope: n = 7 class 1 (the 13 D = 0 key-kept steps, plus all 75 J != 0 steps of the rebuilt 23 000-expansion walk, 19 kept and 56 refused); all LNAs and duals n = 4..7; the F-010 n = 9 pair 3060000 / 3304000 / 2400230 (HH^0..6); planted pairs at n = 4 and 7. HH computed over GF(32003) with the repo's presented coefficients (not Q, not "up to coefficients"). HH^0..2 only on the steps (HH^3 cochains sized but not ranked). No n >= 8 walk.
power: PLANTED NO: commutative square vs square with ab = 0 (same Cartan matrix, hence same Coxeter key), with and without a 3-vertex tail (n = 4 and n = 7): HH(0,1,2) = (1,0,0) vs (1,1,1), NO 2 of 2; real refused steps: HH says NO on 2 of 56 (key differs, so NO is true). Where it must say NO and does not: the F-010 pair (equal key, inequivalent) HH^0..6 = (1,0,...,0) for all three, NO 0 of 1 (and 54 of 56 refused steps, where the key already says NO).

## Claim
Proposed invariant I_HH: for a step A -> B (finite-dimensional bound quiver algebras over a field k) say NO iff dim HH^m(A) != dim HH^m(B) for some m in {0,1,2}; else YES (no information). Mechanism: Hochschild cohomology is a derived invariant (Rickard 1991, Happel 1989), so for every k-algebra pair that is derived equivalent over k the dims agree for all m; a difference is a proof of inequivalence. It is independent of the key and of the Cartan matrix (shown: it separates a pair of equal Cartan, and the square-with-tail pairs have identical Cartan and key). It rejects exactly those J != 0 steps whose child has different HH^m dimensions up to degree 2, e.g. a child that acquires a non-simply-connected cycle (HH^1 != 0) the parent lacks. Results: (1) on the 13 steps it says YES, 13 of 13 (all HH = (1,0,0) for parent and child). That is the truth (the 13 are Hom-joined, E-297, E-300), so it is a consistency check and not a test of the R-candidate: a NO there would have contradicted the chains or the coefficient convention. (2) Planted NO: caught 2 of 2. (3) Obstruction: HH is (1,0,...) for every LNA and dual at n = 4..7 (all 7 key classes, 264 algebras at length 7, theorist_lna_out.txt) and for all three n = 9 algebras with key 1,1,-1,-3,..., including the inequivalent F-010 pair, in every degree (cochain complex ends at degree 6). Of 75 J != 0 steps, parent HH is (1,0,0) in 75 of 75; the child differs in 2 (both refused steps, child HH^1 = 1). So no derived invariant of the HH-dimension kind can say NO to a key-kept step in the n <= 7 walk, and none can separate the F-010 pair.
What would give NO: the walk must reach a child with nontrivial HH, or one needs a finer invariant than dimensions. Not claimed: that the R-candidate is true (the evidence is still YES-side); that HH-triviality holds for every LNA (only n <= 7 and three n = 9 algebras, computed; for tree-quiver algebras HH^1 = 0 is classical, higher degrees are observed only).

## Evidence
The obstruction in words: the key-kept region of the walk consists of algebras whose quiver is a line or its mutation-descendants with HH = k concentrated in degree 0, i.e. algebras with the HH of A_7. Any test using HH dimensions then answers YES on every pair in that region, equivalent or not, so it has no power there. The n at which a NO is needed: a key-equal inequivalent pair is first known at n = 9 (F-010), and there HH dims are blind too. A NO-capable invariant at the required pair must be Cartan-level (q_S, E-292, not key-free in spirit, it is a Cartan datum) or finer than HH dims (Gerstenhaber / BV structure, HH as a graded algebra; for these trivial-HH algebras this is empty too: HH = k forces trivial structure, so those are blind as well). Rickard-type inequivalence of the F-010 pair rests on Ladkani's Cartan argument alone. Weakest step: the HH code is validated by A_n (1,0,0), Kronecker (1,3,0), commutative square (1,0,0) equal to D4, A3/rad^2 (1,0,0); the square with ab = 0 giving (1,1,1) is not independently checked beyond HH^1 = 1 - chi(Q) for gentle algebras; HH^2 value (1) unchecked. Coefficient issue: the chains join A and B "up to coefficients"; HH of presented coefficient-1 algebras agreed for all 13 steps so no conflict appeared.
| set | n | parent HH | child HH | verdict |
|---|---|---|---|---|
| 13 key-kept steps | 7 | (1,0,0) x 13 | (1,0,0) x 13 | YES 13/13 (true) |
| other 6 kept J != 0 | 7 | (1,0,0) | (1,0,0) | YES (truth: not decided here) |
| refused steps | 7 | (1,0,0) x 56 | (1,0,0) x 54, (1,1,0) x 2 | NO 2/56 (true NO; key says NO 56/56) |
| comm. square vs ab=0 (+tail 0,3) | 4,7 | (1,0,0) | (1,1,1) | NO 2/2 (planted, equal Cartan) |
| F-010 pair 3060000 / 3304000 | 9 | (1,0^6) | (1,0^6) | YES (false YES: inequivalent) |
| corruptions of the 13 children (drop a relation or a term) | 7 | | | HH differs 30 of 69; of these Cartan equal to child 13 |
| all LNAs + duals | 4..7 | (1,0,0) | | 264 at n = 7, fewer below: one HH value |

## Reproduction
```
.venv/bin/python workshop/rounds/115/theorist_hh.py                  # self-tests, 1 s
.venv/bin/python workshop/rounds/115/theorist_ctl.py                 # planted pairs, 2 s  (out: theorist_ctl_out.txt)
for n in 4 5 6 7; do .venv/bin/python workshop/rounds/115/theorist_lna.py $n; done   # 2 s (theorist_lna_out.txt)
.venv/bin/python workshop/rounds/115/theorist_f010.py                # n = 9, 3 s (theorist_f010_out.txt)
timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 540 23000   # twice; rebuilds /tmp checkpoints, ~6 min total here
timeout 10m .venv/bin/python workshop/rounds/115/theorist_steps.py   # 13 steps + corruptions, ~2 min
timeout 10m .venv/bin/python workshop/rounds/115/theorist_allJ.py    # all 75 J != 0 steps, ~2 min
```

## Prior record
E-298 pre-registered that a NO needs a derived invariant with its own controls and named HH as an example without computing it; E-295 (Cartan Z-congruence 0 of 89), E-212 (no key-equal inequivalent pair at n <= 7; F-010 pair at n = 9). grep of research/ for "Hochschild" to be confirmed by the chair; not found as a computed invariant here. New: HH^0..2 is computed, YES on 13/13, a same-Cartan planted NO, and the proof-by-computation that it is blind on LNAs and on the F-010 pair.

## Code changed
None in the library. New: `workshop/rounds/115/theorist_hh.py` (E-relative reduced bar complex, mod p), `theorist_ctl.py`, `theorist_lna.py`, `theorist_steps.py`, `theorist_allJ.py`, `theorist_f010.py`. No tests exist for them.

## Next
Skeptic: an independent HH check (the HH^2 of ab = 0 square, a second implementation or Bardzell/gentle formula). Experimentalist: a walk that leaves the HH-trivial region (n >= 8 with parallel arrows or oriented-cycle-free non-tree quivers) so HH can speak; record HH(parent), HH(child) of every J != 0 step in the next walk. Scholar: whether HH of a quiver-tree algebra vanishes in all positive degrees (Happel?).
