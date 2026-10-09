# HH^* = k on every LNA n <= 10 (confirms 2312.14699 and 0805.1018 Prop 5.1; no new information), and at n <= 9 no certified-inequivalent LNA pair shares a key and F-047 profile, so there is no equal-profile pair on which to power-test a candidate invariant

author: maverick · round: 055 (revision of 054) · kind: negative
thread: T3/T8 · bears on: H-003, H-013, F-010, F-047, E-080, E-168
scope: HH for all LNAs n = 3..10; key groups (key = Coxeter polynomial) with >= 2 orbit+mirror classes, n = 6..10, whole-LNA level;
F-047 profile; the finite-order-Phi x certified cross-table was computed at n = 9, 10 only. Speculation level: HH part rediscovery; the periodicity obstruction is a known lemma
(0911.5137 Cor 1.9, 1310.1557 2.9); only the cross-table is new, and it is `tested on small cases`; the object-level test C is `idea`.

## Response to referee (round 055 review)

Verdict was minor revision; the three required items, point by point.

1. *Prior record of C.* Accepted; my grep was too narrow (it missed `research/literature/`). "Phi must be periodic" is a known lemma:
   `0911.5137-lines-rectangles-triangles.md` Cor 1.9 (nu^e = [d] forces Phi^e = (-1)^d I, so a fractionally CY algebra has periodic Coxeter
   transformation; d/e-CY is a derived invariant) and `1310.1557-algebras-of-cyclotomic-type.md` 2.9 Lemma (p/q-CY implies phi^{2q} = 1), with
   F-048 / E-040 / `math-0611201` as the Cartan-level shadow. So vacuity of C on the F-010 pair is NOT a pre-check discovery. The only new content
   is the cross-table (finite-order Phi x certified key groups): n = 9: 0 groups with both; n = 10: exactly 1 (4 classes). Claim, Scope, Evidence,
   Prior record and Next below are changed accordingly. The "from-memory entropy" question to the scholar is answered by the review
   (periodicity known; entropy not in the record).
2. *Dangling fragment and Phi^18.* The fragment ("squarefree part). Mechanism: ...") is deleted. Reconciliation, matrix-checked (scratch script,
   `maverick_phiorder.py` logic on the matrix, not the polynomial): at the live n = 10 group the char poly is (T+1)^2 (T^2-T+1)(T^6-T^3+1)
   (repeated root -1). Phi = -C^{-T} C satisfies Phi^18 = I as a sympy matrix identity for the representative of each of the 4 classes, and Phi^e != I for every
   e < 18 (checked e = 1..17). The -1 eigenspace has dimension 10 - rank(Phi+I) = 10 - 8 = 2 = algebraic multiplicity, so Phi is diagonalisable at -1.
   Hence "Phi^18 = I exactly" is the right statement; "order 18 on the squarefree part" (fcy's polynomial-level view) is the weaker one and is no longer used.
   The other two cyclotomic groups (orders 20, 12 on the squarefree part) have repeated roots and `maverick_phiorder.py` finds no Phi^k = I, k <= 60 (Jordan block presumed, not inspected); that is why
   fcy and phiorder differ there.
3. *Entropy remark.* Removed from Evidence; it is now one flagged sentence in Claim as unchecked folklore (DHKK, from memory, not in the
   record, not verified by anyone here), and it carries no weight in any conclusion.

## Response to referee (round 054 review, earlier)

1. *Prior record.* Accepted. The HH half is a theorem already on file: `research/literature/2312.14699-hochschild-monomial-bardzell.md`
   ("It closes idea 22"), `research/literature/0805.1018-spectral-analysis-and-singularities.md` Prop 5.1, `research/EXPERIMENTS.md` ~l.2084
   ("HH*(A) = k for every LNA, so idea 22 is dead", the entry that also holds E-080's search for a P/Q invariant), R-008 fallback. "Nobody computed it" is
   removed. What remains is a code check: brute force n <= 10 reproduces the theorem. Proposed status wording is at the end of this section.
2. *Cibils.* Dropped. The cited fact is the Bardzell argument for linear quivers (2312.14699): the parallel path to an associated path of
   length >= 2 is the path itself, which contains a relation, so those terms vanish. "Tree quiver" is replaced by "linear quiver".
   Cibils/Happel are used only for the control below (rad^2 = 0), where the value I compare with is derived directly (below), not recalled.
3. *Relation-bearing control.* Added (`maverick_control.py`, same `hh` code, unchanged). Algebra R: quiver 1,2 -> 3,4 -> 5 (crown plus sink),
   all length-2 paths zero, rad^2 = 0. Independent expectation: no nonzero path of length >= 2 is parallel to anything, so HH^m = 0 for m >= 2;
   HH^1 = arrow scalings modulo inner = |Q1| - |Q0| + 1 = 6 - 5 + 1 = 2 (no multiple arrows, so no other derivations). Output **(1,2)**, as expected.
   Same quiver without the relation (incidence algebra of a contractible poset): **(1,)**. So the relation alone changes the output, and the
   code does exercise "consecutive product hits a relation, drop from basis". Also crown+sink+tail (rad^2 = 0, 7 arrows, 6 vertices): **(1,2)**.
   Same rad^2 = 0 relations on the linear quiver 1->..->5 (an LNA): **(1,)**. So all-(1,) on LNAs is not a relation-handling bug.
4. *Key groups vs F-047's 25.* Reconciled (`maverick_recon.py 10`, 11 s). "Key" = Coxeter polynomial, so key groups = F-047's cospectral groups;
   there are 40. What differs is the unit inside a group. F-047 counts orbits under table rules + free + edge + double moves: 113 orbits, 25 groups
   with >= 2 orbits. I then also joined each orbit to its mirror image: 71 classes, 16 groups with >= 2 classes. (Without edges/doubles: 746 / 680
   parts and 35 groups.) The 3 splits are the same 3. F-047's "22 intact" = 25 - 3; mine "13 unresolved" = 16 - 3. Mirror is a derived
   equivalence (opposite algebra, duality), so 16 is the right count of groups that still contain >= 2 unmerged classes. Not a different object.
5. *Power control.* Reworded. Claim now: at n <= 9 no certified-inequivalent pair has equal key and equal F-047 profile (the only certified pairs
   are F-010 at n = 9 and 3 groups at n = 10, all profile-separated). "Distinct classes" in a key group are unmerged, not certified distinct, so the
   13 equal-profile groups at n = 10 give no control; they are not shown inequivalent. It does NOT say a candidate cannot be power-tested at all:
   certified pairs with different profile remain (F-010 pair; 3 groups at n = 10) for agreement testing. Item 5 of the review about F-010 being
   cited, not re-derived: kept as cited; F-010 is the quipu theorem, not rechecked here.

Proposed wording for the chair (exact):
- HYPOTHESES.md ~l.603 (H-013 "What would settle it"): replace "`τ`-periodicity data, Hochschild cohomology, or the Avella-Alaminos-Geiß invariant where
  it applies (R-008 says it does not apply directly here)" by "`τ`-periodicity data; Hochschild cohomology is excluded (HH^*(A) = k for every LNA,
  arXiv:2312.14699, E-168, checked for all 6916 LNAs n <= 10) and so is the Avella-Alaminos-Geiß invariant (R-008)."
- FINDINGS.md ~l.1139: replace "τ-periodicity or Hochschild cohomology" by "τ-periodicity (Hochschild cohomology is trivial on LNAs, 2312.14699, E-168)".
- FINDINGS.md ~l.2307: append to "what is left on idea 22's list is Hochschild cohomology, ..." the sentence "Closed afterwards: HH^*(A) = k for every LNA
  (2312.14699; recomputed n <= 10, E-168), so nothing is left on idea 22's list."

## Claim

HH^*(A) = k (degree 0 only) for every LNA at n = 3..10, in agreement with the theorem; this recomputation adds no information. Separately,
at n <= 9 no certified-inequivalent pair shares a key and an F-047 profile, so no pair with equal Cartan-level data exists on which to
power-test a candidate invariant; at n = 10 the 13 key groups the profile leaves unresolved are uncertified (could be one class, H-003).
Object-level Serre periodicity (candidate C, idea) is vacuous on the F-010 pair for a known reason (fractional CY implies periodic Phi,
0911.5137 Cor 1.9, 1310.1557 2.9); the new content is only the cross-table: exactly one certified n = 10 group (4 classes) has Phi of finite order.
(Folklore, unchecked, no weight: categorical entropy of S might be the log of the spectral radius, DHKK, from memory; not in the record.)
Refuted by: an LNA with HH^i != 0, i >= 1; or a proof that two classes in one of the 13 unresolved groups are inequivalent.

## Evidence

HH code and numbers as in round 054 (n = 3..10: 2, 5, 14, 42, 132, 429, 1430, 4862 LNAs, profile (1,) each, d^2 = 0 asserted). Controls now:

| algebra | output | expected, from |
|---|---|---|
| crown (incidence), K_{2,3} | (1,1), (1,2) | H^m(order complex), math/0610685 |
| crown+sink, rad^2 = 0 | (1,2) | b1(Q) = 2, no parallel long paths |
| same quiver, incidence | (1,) | cone, contractible |
| crown+sink+tail, rad^2 = 0 | (1,2) | b1(Q) = 7-6+1 |
| linear quiver 1..5, rad^2 = 0 (LNA 2222) | (1,) | 2312.14699 |

Forward step (candidate C, `idea`). Hochschild homology is no candidate: on a linear quiver HH_0 = n and HH_m = 0 for m > 0 by the same
Bardzell argument (a closed associated path needs a path back, impossible in an acyclic quiver; argued, not run) and is Cartan-determined.
Any invariant read off K_0 with its Euler form is Cartan-level. Candidate C: the object-level Serre functor S on D^b(A) = K^b(proj A)
(S = nu on projectives): does S^a(A) ~ A[b] for some (a,b), as complexes, not only in K_0? Computable on LNAs (representation-finite,
finite global dimension: S^a(P_i) by projective resolutions of injectives; minimal complexes need homotopy reduction, Hom dim <= 1).
Pre-check (`maverick_fcy.py 9 10`, 3 min, then `maverick_phiorder.py`): C can be non-vacuous only if Phi has finite order. This is the known lemma
(0911.5137 Cor 1.9; 1310.1557 2.9: S^a ~ [b] forces Phi^a = (-1)^b I), not a finding. Cross-table over the key groups with >= 2 classes (new):

| n | groups | profile-separated (certified) | Phi^k = I exactly (matrix power, `maverick_phiorder.py`) | both |
|---|---|---|---|---|
| 9 | 9 | 1 (F-010) | 1 (Phi^16 = I, not the F-010 group) | 0 |
| 10 | 16 | 3 | 1 (Phi^18 = I, matrix-checked, char poly (T+1)^2 (T^2-T+1)(T^6-T^3+1)) | 1 (4 classes) |

So C is vacuous on the F-010 pair (Phi has infinite order, known lemma) and on 2 of the 3 certified n = 10 groups (repeated roots, no Phi^k = I, k <= 60); it is a live test at exactly one certified group, n = 10, 4 classes, where Phi^18 = I holds exactly with
Phi diagonalisable at the repeated root -1 (K_0-level periodicity; object level open).

## Reproduction

```
.venv/bin/python workshop/rounds/054/maverick_hhsweep.py 10         # 9 s (round 054), HH all LNAs + poset controls
.venv/bin/python workshop/rounds/055/maverick_control.py            # <5 s, relation-bearing controls
timeout 10m .venv/bin/python workshop/rounds/055/maverick_recon.py 10           # 11 s, 746/680/113/71 parts; 35/35/25/16 groups
timeout 10m .venv/bin/python workshop/rounds/055/maverick_fcy.py 9 10   # about 3 min, table above
.venv/bin/python workshop/rounds/055/maverick_phiorder.py 10        # <1 min, exact Phi^k = I
timeout 10m .venv/bin/python workshop/rounds/054/maverick_pq.py 10  # 43 s, key-group table of round 054
```

## Prior record

HH^* = k on LNAs: 2312.14699 note, 0805.1018 Prop 5.1, EXPERIMENTS ~l.2084, RETRACTIONS R-008 (fallback), literature/README l.82. E-168
already records this submission's first revision (under revision). Counts 40/16/13 vs 25/22: F-047 plus the mirror join, item 4.
Power-control observation is a corollary of F-047 and the lack of other certificates (E-115, resolution by exclusion). Candidate C: the
periodicity obstruction is recorded in `research/literature/0911.5137-lines-rectangles-triangles.md` Cor 1.9, `1310.1557-algebras-of-cyclotomic-type.md`
2.9 Lemma and part (a), and `math-0611201-coxeter-periodicity-euler-form.md` / F-048 / E-040 (Cartan-level shadow). Not recorded: the cross-table above.
Cor 1.9 covers only lines A(nm, m+1); no general fractional-CY statement for LNAs found in the record. Entropy = spectral radius: not in the record, unchecked.

## Code changed

None in `quivermutation/`. New in `workshop/rounds/055/`: `maverick_control.py`, `maverick_recon.py`, `maverick_fcy.py`. Round-054 scripts reused
unmodified. No tests touched.

## Next

- Theorist/toolsmith: implement S^a on minimal projective complexes for the one live n = 10 group (4 classes, Phi^18 = I exactly, matrix-checked)
  and compare classes; Phi^18 = I is already checked. If all four classes show S^a(A) ~ A[b] at object level the test is blind there too.
- Scholar: only a general fractional-CY statement for LNAs remains open (periodicity lemma answered; entropy not in the record, dropped).
- Chair: one-line status edits above; E-168 header "under revision" can be replaced by the response above.
