# H-010 has no proof in reach: the interior-invariance statement quantifies over non-LNA intermediates, and what I can test is a run-of-three lemma, which holds at 3 mutations

author: theorist · round: 045 · kind: proposal (with a negative on the proof, and a small computation)
thread: T7 · bears on: H-010, F-022, F-024, E-021, E-025, E-029, E-068, E-124, H-011
scope: pair `(8:3)(9:3)` in a line with >= 6 empty vertices on each side; one bystander relation of m = 2, 3, 4 arrows; <= 3 mutations (m = 2 and m = 3 up to gap 1 only, run hit the 10 min cap), mutated vertices within 2 of pair + bystander, left and right mutations. No proof. Step 7 was not run (`--steps 7` not used).

## Claim

(1) I cannot prove H-010 from the procedure in one sitting, and I think the missing piece is stated precisely enough to close T7 as it stands (section "Why no proof"). (2) The statement the data actually support is sharper than "isolated pair": **Lemma L (run-of-three criterion).** If a pair of relations sharing >= 2 arrows is planted in the interior, a bystander relation lowers the maximum overlap within k mutations only if it shares >= 2 arrows with the pair (a run of three). Tested for k <= 3, one bystander, m <= 4: no lowering for any bystander that shares <= 1 arrow or none (22 of 22 such placements at k = 2; 10 of 10 that finished at k = 3), and lowering for `(8:3)(9:3)(10:3)` (k = 3: 4 of 6 reached LNAs are lower, min overlap 0; k = 2: 2 of 4, min overlap 1); the other run-of-three placement, `(8:3)(9:3)(10:4)`, did not lower at k = 2 (3 reached). It is false as a general statement about long runs (E-025), and I claim nothing for k > 3, two bystanders, or a bystander on the left (not run; F-025 says the two ends are not symmetric for unequal pairs, so I did not assume it).

## Why no proof (what is missing)

Statement to prove: for every sequence of admissible mutations of any length, from an LNA with an isolated run of two in the interior, every LNA reached has maximum overlap >= 2.

1. *The quantifier runs through non-LNA algebras.* Every intermediate is a quiver with relations that need not be a line. The procedure (F-015, `procedure.py`) is defined on all of them, so a proof needs an invariant defined on that whole class. I found none. The Cartan matrix and Coxeter polynomial are constant on the derived class, but that class contains LNAs with overlap 0 (F-022, six mutations with an end), so they cannot separate.
2. *Locality is available, invariance is not.* By E-124 the gate at v is `J_i = Hom(S_v, e_iA) = 0` for all i, and a nonzero element of `J_i` is a nonzero path p: i ~> v with p.b = 0 for every out-arrow b. In an LNA that means a relation ends exactly at the arrow out of v. So a right mutation is admissible only at vertices that are not "one before a relation's target" (dually for left), and step 7 only touches paths through v. Hence a depth-k sequence in a window of the line depends on the algebra only inside a window of radius about k, and is translation covariant. This gives: the interior orbit at depth k is independent of n and position (consistent with 2, 4, 4, 6 LNAs at depth 3 to 6, F-024). It does *not* give the statement for all k, because the window grows with k and the orbit on the infinite line is then the object, and a classification of that orbit is H-011's residue (F-032), not something the procedure hands over.
3. *The conserved quantity cannot be a naive overlap count.* E-029: 1416 of 7164 interior double mutations at n = 10 lower the maximum overlap when bystanders cross the relation. So "overlap invariant under interior mutation" is false without the isolation hypothesis, and with it the hypothesis must be preserved by the dynamics. Showing that it is (an isolated run stays isolated, or every way it stops being isolated is itself a recorded rewrite) is exactly the classification in point 2.
4. *Step 7 as a concrete test case does not bear on this.* The step-7 element of E-068 (commutativity `8,6,4,9 + 8,10,4,9 = 0`, kernel in `e_8 A e_4`) is the reason a mutation is *rejected*. It sits in a non-LNA algebra (a commutative square), is a statement about J, and says nothing about which LNAs are reached afterwards. It is not an instance of H-010. It could only matter if a proof of invariance had to leave the LNAs and pass through such squares; that is the claim in point 1, and I did not get further.

Recommendation: close T7 as "no proof available; H-010 stays SUPPORTED by search to depth 6 (F-024), now with Lemma L as the form to cite for bystanders". Reopen only if someone proposes an invariant on non-LNA algebras. The closest proof target is H-011's classification of the interior orbit on the infinite line, which T-threads on F-032 already cover.

## Evidence

Pair `(8:3)(9:3)` (relations at vertices 8 and 9 of three arrows, arrows 8..11), bystander `(t:m)` with t = 12 + g. g is the number of empty arrows between them (g < 0: shared arrows). Line length = end of bystander + margin + 6. `lowered` = LNAs reached with smaller maximum overlap than the start (2). Table is `workshop/rounds/045/theorist_bystander_3_2.txt` (k = 3) and the k = 2 run below.

| m | g | shared with pair | k = 2: reached / lowered | k = 3: reached / lowered |
|---|---|---|---|---|
| 2 | -1..5 | 1 or 0 | 4,5,6,6,6,6,6 / 0 | 9,11,11,12,12,12,12 / 0 |
| 3 | -2 | 2 (run of 3) | 4 / **2** | 6 / **4** (min overlap 0) |
| 3 | -1 | 1 | 1 / 0 | 1 / 0 |
| 3 | 0, 1 | 0 | 2, 2 / 0 | 2, 2 / 0 |
| 3 | 2..5 | 0 | 2 / 0 | not reached (cap) |
| 4 | -2 | 2 (run of 3) | 3 / 0 | not reached (cap) |
| 4 | -1..5 | 1 or 0 | 1,2,2,2,2,2,2 / 0 | not reached (cap) |

A run of three is therefore necessary for lowering in this data but not sufficient (`(10:4)` at k = 2 stays at 2; cf. F-022, where `(1:4)(2:4)(3:4)` stays at 2). The placement m = 2, g = -2 is inadmissible and was skipped.

What this adds to F-022/E-025: the bystander does not need to be far to be inert. A bystander sharing exactly one arrow, or adjacent with none, behaves like an absent one for three mutations even though mutations were allowed on its vertices. What it does not add: depth. Reached counts for the isolated pair agree with F-024 up to the bystander's own contribution.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/045/theorist_bystander.py 2 2   # about 1 min
timeout 10m .venv/bin/python workshop/rounds/045/theorist_bystander.py 3 2   # hits the cap after m = 3, g = 1 (about 10 min)
```

The k = 2 table is also the first 23 lines of the script's output; only the k = 3 output is committed as a file.

## Prior record

F-022 (run of three unlocks; E-025 qualification), F-024 (interior to depth 6), E-029 item 8 (bystanders crossing r lower overlap, isolated pairs only slide), H-010 caveat paragraph (the "right form of the statement is open"). Lemma L is the same observation as F-022's table made with a bystander at varying distance rather than at fixed overlapping positions; the new part is the dependence on shared arrows (>= 2 needed, 1 or 0 inert) at k <= 3. Not new: that a proof is hard (the H-010 text already says the caveat carries the weight). E-124's gate formulation (J_i) is used in point 2.

## Code changed

None to `quivermutation/`. New script `workshop/rounds/045/theorist_bystander.py`. No tests run (no library file touched).

## Next

- Skeptic: the m = 4, g = -2 row and the cap-cut rows; rerun them at k = 3 with the margin lowered to 1 to finish the table.
- Theorist (next breadth slot): state L1, the one-step statement "an interior mutation of an LNA at v, when the result is an LNA, changes relations only within distance 1 of v" and test it over all LNAs n <= 9; if true it is the locality half of any proof.
- Chair: decide whether T7 is closed. If kept, the open object is an invariant on non-LNA algebras, not a search.
