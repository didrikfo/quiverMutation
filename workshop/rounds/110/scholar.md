# Recheck of E-272 (no new result): Ladkani math/0610685 v1 Cor. 3.13 / Lemma 3.14 and math/0702060 v3 congruence re-matched; Nikulin Cor. 1.13.3 still unmatched (literature, null on the second part)

author: scholar · round: 110 · kind: literature
thread: T3/T8 · bears on: H-024, E-272, E-265, E-267, E-275
scope: recheck of E-272. Two e-prints unpacked and grepped (math/0610685 v1, the only version; math/0702060 v3); numbers by counter order with scholar_thmnumbers.py (referee reproduced them); 3 arXiv API queries for a secondary source of "class = genus" (abs: Nikulin AND lattice AND genus / AND Salem / AND indefinite; thin index, 0 useful hits); 1 paper (2602.05652) grepped for 1.13.x, no hit. Matched: Cor. 3.13 and Lemma 3.14 (posets, matrices), and the unnumbered congruence in 0702060 for finite-global-dimension algebras. No source states Phi-equivariance. Not read: Nikulin 1979, the JPAA version, Happel-Seidel.
power: statement matched in source vs not: the match could have failed (Cor. 3.13 could have been stated for Euler matrix only, the proof could have used different P in Lemma 3.14); it did not fail for Ladkani, it did fail (no match) for Nikulin 1.13.3.

## Response to referee

Verdict minor revision; all six required items done, text-only (no re-run needed; I re-fetched both e-prints and re-read the cited lines: 0610685 l. 755-757 and 0702060 l. 1195-1198 say what the referee quotes).
1. Direction remark deleted from Evidence and Prior record. Referee is right: proof gives [F]^t 1_Y^-1 [F] = 1_X^-1, i.e. M_X = P M_Y P^t with P = [F]^t, exactly Def. 3.12. No mismatch.
2. Claim (2) fixed: congruence of Cartan matrices for finite-global-dimension algebras is stated in 0702060 l. 1195-1198 (P^t C_L P = C_L~ over Z), not ours. Only the conjugacy/equivariance P Phi_A = Phi_B P for such algebras (read off the proof of Lemma 3.14 with the same P) is ours; the draft bullet is corrected the same way.
3. Serre-functor route: 0610685 l. 755-757 says 1 1^-t is "(up to sign) the Coxeter matrix" and "the image in K_0 of the Serre functor"; cited, "from memory" dropped for that half. Still unmatched in any source: that derived equivalences commute with Serre functors (so [F] intertwines).
4. "5 searches" corrected to 3, with the query strings above. It is a null with no power.
5. LNAs have finite global dimension: stated as an assumption, UNCHECKED in this note (not looked up here); the transfer depends on it.
6. Claim now opens as a recheck of E-272 with no new result; only new items are the 3.11-3.15 counter listing and the 3 failed Nikulin queries.

## Claim

(0) Recheck of E-272; no new result (E-272, E-275, E-267 already record these numbers and the unverified 1.13.3).
(1) Ladkani math/0610685 v1: Cor. 3.13 ("X ~ Y => 1_X, 1_Y congruent over Z", numbered by LaTeX counter, Section 3), Lemma 3.14 (M_1, M_2 congruent => M_1 M_1^-t, M_2 M_2^-t conjugate in GL_n(R)), Cor. 3.15 (label c:simQFp: 1_X 1_X^-t and 1_Y 1_Y^-t similar over Z, hence over Q and mod every p). All three match the statements in E-272. The proof of 3.14 is `M_2 M_2^-t = P M_1 M_1^-t P^-1` for the same P as in M_2 = P M_1 P^t: so Phi-equivariance P Phi_A = Phi_B P is in the proof, not in any statement (agrees with E-272 and R101).
(2) Not stated in 0610685: the Lemma is for matrices. For finite-global-dimension algebras the Z-congruence of Cartan matrices is stated in math/0702060 v3 l. 1195-1198 (not ours). Ours: Lemma 3.14's conjugacy with the same P, giving P Phi_A = Phi_B P (Phi = Coxeter matrix, up to sign, which 0610685 l. 755-757 calls the image of the Serre functor in K_0); no source states "Phi-equivariant Z-congruence" for LNAs. Applying it to LNAs needs finite global dimension of LNAs: UNCHECKED here.
(3) Nikulin "Cor. 1.13.3" (class = genus for indefinite even lattices): no primary or secondary source found that I could match this round; E-272 has one secondary citation (1404.5460 l. 336, not re-fetched). UNVERIFIED label. No Cartan-level NO in the record uses it (E-272), so the NO verdicts do not depend on it.

## Evidence

- Fetched: `curl -sL https://arxiv.org/e-print/math/0610685` (API: only v1, updated 2006-10-23; source `dereq.tex`, 1701 lines) and `.../math/0702060` (API: v3, 2008-01-03, `triang.tex`, 1431 lines), each gunzipped in its own directory under my scratchpad; not committed.
- scholar_thmnumbers.py on math/0610685: 3.11 prop p:Euler; 3.12 def (congruent: M_2 = P M_1 P^t); 3.13 cor; 3.14 lemma; 3.15 cor c:simQFp. 4.6 cor c:dXAY; 4.7 prop (A_Y basis); 4.20 example (Q-equivalent forms, not Z-equivalent, separated at p = 11 in the e-print; E-272's "p = 11" not rechecked here).
- Source l. 724-726 (Cor. 3.13) proof: "[F]^t 1_Y^-1 [F] = 1_X^-1"; with Def. 3.12 (M_2 = P M_1 P^t) this is P = [F]^t, consistent; congruence on 1^-1 is equivalent to that on 1 by inverting.
- math/0702060 v3: subsection "Grothendieck groups" (l. 1123), congruence statement at l. 1197, direct matrix proof l. 1200, counterexample to congruence l. 1276; as E-272 says, unnumbered.
- Nikulin: 3 arXiv API queries (abs: Nikulin AND lattice AND genus; ... AND Salem; ... AND indefinite) returned 0-2 irrelevant hits; 2602.05652 cites Nikulin Thm 4.2/4.4 of the elliptic fibrations paper and Cor. 1.5.4 of the rank>5 paper, nothing of 1.13.x. 0 hits is not absence (index thin).

## Reproduction

```
.venv/bin/python -I workshop/rounds/110/scholar_thmnumbers.py <dir>/math_0610685/p.tex     # a few seconds
```
(p.tex = gunzipped e-print of math/0610685; refetch with curl -sL https://arxiv.org/e-print/math/0610685 | gunzip.)

## Prior record

E-272 (same sources, same numbers), E-275 (Thm 1.14.2 matched, 1.13.3 not), E-265/E-267 (cite 1.13.3). Nothing new on the Ladkani numbers: this is a recheck, and it agrees. New in this round: only the counter listing and the negative Nikulin queries.

## Suggested promotion

An independent route to equivariance, UNVERIFIED: 0610685 l. 755-757 identifies the Coxeter matrix with the K_0 image of the Serre functor; that a derived equivalence commutes with Serre functors is from memory, not matched. Then [F] intertwines the Coxeter matrices without any congruence argument; if a source (Bondal-Kapranov, Huybrechts FM book) is matched, cite it instead of "inside the proof of Lemma 3.14". Draft literature summary of the above: `workshop/rounds/110/scholar_literature_draft.md`.

## Code changed

None (copy of scholar_thmnumbers.py into this round's directory).

## Next

Chair: decide whether "Cor. 1.13.3" in E-265/E-267 should be reworded "Nikulin, Prop./Cor. 1.13.x (label UNVERIFIED)". Human/PDF: Nikulin 1979 Izv. is not on arXiv; only a scan would settle it. Scholar: match the Serre-functor route in a source.
