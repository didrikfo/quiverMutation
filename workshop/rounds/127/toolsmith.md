# No gate-admitted step has a loop at i, and on 12 136 admitted steps from cyclic monomial parents (n <= 8, capped families; 4 343 with a cycle reaching i; 1 949 A^op, all with the cycle beside i) (LF), (H1w), (H3), H(I) = 0 all hold, with Theta' and lead words equal to those of the D(i) part; the library gate and `isTiltingAtVertex` give wrong answers when a cycle lies upstream of i

author: toolsmith · round: 127 (revised after review) · kind: result · thread: T5 · bears on: H-023, H-031, E-326, E-279, E-211, E-199, E-130
scope: monomial parents only, path-length weights, finite gl.dim; A side (= F step, E-252) and A^op (= R step). A^op admits only steps whose D^op(i) misses the cycle: 0 at n = 4, 30 at n = 5, 1 721 in the glued family nL <= 6, 198 in the long-glue family; all run, all pass, none has a cycle reaching i. (1) Full enumeration up to isomorphism: n = 4 and n = 5, at most 5 arrows, no parallel arrows, at least one oriented cycle, at most 4 zero relations of length 2..3 (n = 4 also length <= 4). (2) Glued family: an LNA on 1..nL (every antichain of zero relations, nL = 2..6), plus a 2-cycle c1 <-> c2 attached by g: c2 -> v, cycle relation ab, ba, aba or bab, and at most one glue relation of length <= 3 (ag, bag, g·x, ag·x, g·x·x). (3) Long-glue family: the same, with exactly one glue relation of length 4 (bag·x, ag·x·x, g·x·x·x). n = nL + 2 <= 8. Of the 12 136, 4 343 have a cycle reaching i; the other 7 793 have the cycle beside D(i), where the check is an acyclic check in disguise. Not covered: n >= 6 outside the glued families, cycles of length >= 3 for n >= 6, non-monomial parents, |A| >= 2 with a cycle, Hom/End(T).
power: [GLDIM] Anick tail-graph test vs 126's Kupisch test on 360 Z_n antichains: 0 disagreements, 148 infinite. [MM] one dropped d-term gives 25 H = A failures (NO). [GATE] planted relation x3.x4 refused by both gates (NO). [ENG] the closed ideal (stars) gives 8 bad cells (NO). [SH] a planted g is caught (NO). [BARE] glued-vs-bare says "different" against A4 20, A4 02 and "same" against A4 30 (truth on all three). [SIMPLE] exact gate NO, library gate YES. [AUDIT] the D(i)-closure test says "not closed" when D(i) is replaced by a set that is not successor-closed (NO). "lib gate wrong" can say YES only where some relation has length >= 4 (n = 4 rel <= 4; long glue); elsewhere its 0 is forced (Response 2).

## Response to referee

1. **A^op sentence, the split, the glued A^op admits.** Done. The sentence was false and is replaced (scope). `glue()` now runs both sides. Glued nL <= 6: 9 891 admitted (8 170 A + 1 721 A^op), 9 891 pass, 0 A^op with the cycle reaching i. Sweeps now tally the side: n = 5 has 30 A^op of 342, 0 with the cycle reaching i, as you counted. The headline is split by "cycle reaches i" (4 343 of 12 136) in the title, scope and table.
2. **Forced zeros.** Done, and tested beyond the list. The table marks the rel <= 3 zeros as forced: a non-simple path into i through a 2-cycle has length >= 3, so a relation that kills it together with an out-arrow has length >= 4. Your instance x1: 1>2, x2: 2>3, a: 4>5, b: 5>4, g: 5>1, rels a.b, b.a.g.x1, i = 1: finite gl.dim, D(1) acyclic, exact gate False (refusing path b.a.g), library gate True, library J True (`bag`). New long-glue family (one glue relation of length 4): library gate wrong on 175 instances, J wrong on 485 (0 on admitted steps); all 1 856 admitted pass every check, and Theta'/leads equal bare on 815 / 815. I agree that "J wrong never on an admitted step" is structural (both routines see a subset of paths, so they err only toward True). It is now stated as a reason, not a finding.
3. **Lemma (c).** Done: proof from Def. 3.2 (your argument, below), scope narrowed, recorded as an extension of E-199. On "confirm that (H1w) and (H3) inspect only D-sourced letters": for (H1w) yes; for E-326's (H3) code, not literally. `shape()` also tests "g drops along d_H" on every X letter, including ones outside D(i). That test cannot fail for the recursive g of `part_shape` (g(x) = 1 + max over the words of d_H x), so it carries no information there. Machine check (`daudit`, all 12 136): the D(i)-sourced letters of kM are closed under d; every Theta letter, every linear-part letter of a Theta letter and every lead word is D(i)-sourced; (H1w), Theta' and the leads computed on the D(i)-sourced letters alone equal those on all letters. Every instance has non-D letters (12 136 / 12 136), so the audit is not vacuous.
4. **Cosmetic.** Done. The tag is now `[cycle reaches i, relation through i]`. The test carries a "PINS A BUG" comment saying that the True must flip to False when the routine is fixed. pytest: 13 passed.

## Claim

**(a) A gate-admitted loop-at-i parent does not exist, at any n.** This is E-279 (1) and is not new. A loop letter of degree >= 1 at i puts i on an oriented cycle, so D(i) is cyclic and `mutateAtVertex` refuses. Library walks never expand a cyclic parent (`noCycles`, E-211; search.py l. 403, 446). So H-023's loop joint is vacuous for walk steps.

**(b) Admitted steps on cyclic parents exist, with the cycle upstream of D(i) or beside it.** "Admitted" means the exact gate (every nonzero path into i, not only vertex-simple ones) and D(i) acyclic. For a monomial algebra, the exact gate is J = 0, because the images p.beta of distinct nonzero paths are distinct paths. Smallest example: n = 4, x1: 1->2, x2: 2->1, x3: 2->3, x4: 3->4, x1.x2 = 0, i = 3. It has the loop letter r1: 1 -> 1 of degree 1. On all 12 136 instances in the scope, every check passes:
- [MM] H = A;
- d^2 = 0, and d is homogeneous;
- (C);
- (LF), which is not vacuous: the kM letter quiver has cycles on all instances;
- (H1w), [PI], and (H3) with E-326's recursive g;
- H(I) = 0 in every cell up to shifted weight 2·(top + 6).
The evidence that bears on the cycle is the 4 343 instances with a cycle reaching i.

**(c) Lemma (extension of E-199 from the step-7 kernel to the model).** Let D = D(i) be acyclic. It is closed under successors. Then the letters of the mutated model kM with source in D form a dg subalgebra. That subalgebra is the mutation of the model of kQ_D/(relations inside D).

*Proof* (Def. 3.2 as summarised in `research/literature/1504.02617-quivers-for-silting-mutation.md` l. 48-84; the referee's argument). A letter with source in D is one of three kinds:
- an old letter inside D;
- alpha*: t(alpha) -> i, with t(alpha) in D;
- phi.alpha^-1: t(alpha) -> t(phi), with s(phi) = i, so both ends are in D.

It is never alpha.phi (Step 2): that letter's source is s(phi), and phi ends at i, so s(phi) in D would put i on a cycle. Step 4 is empty, since there is no phi: i -> i. The differential of each of these letters is built from d(phi), d(phi)/alpha, the derivative of phi with s(phi) in D, and the decomposition terms. The branch t(phi) = i, with its sum of alpha*·alpha.phi terms, never applies, because there s(phi) would lie in D. So every word starts in D, and the triangular reduction A = [[e_D A e_D, 0], [*, *]] gives the subalgebra.

*Scope of the lemma:* it covers the D-sourced letters only. The letters that carry the cycle are the Step-2 letters alpha.phi and d(phi) for t(phi) = i. "The stars live on D(i)" is wrong as I first stated it, because the term sum alpha*·alpha.phi in d(phi) does not live there. (H1w) reads Theta letters and their linear parts, all D-sourced, so it sees only the D part. (H3) is the same once the vacuous X-test of Response 3 is set aside. H(I) sees the upstream letters and is not covered.

*Machine check:* on 4 165 / 4 165 glued instances whose cycle reaches i (1 722 with a relation through i from the cycle side), Theta' and the lead words, normalised by what each letter is, equal those of the bare LNA at the same i. In the glued family nL <= 6, the leads have E-326's two shapes: x·q, 1 748 of 2 336, and x*·p, 573 A + 15 A^op = 588. These shapes are matched from the printed tally. I did not check by program that x = x_(i+1).

**(d) Library: two silent errors when a cycle lies upstream of i.**
- **The gate admits wrongly.** This happens on 8 / 137 finite-gl.dim n = 4 parents (relation length <= 4) and on 175 long-glue instances. The refusing path is not vertex-simple, for example x2.x1.x3 or b.a.g. Smallest case: x1.x2 = 0, x2.x1.x3.x4 = 0, at 3. `mutateAtVertex` then runs. The gate emits `CyclicQuiverWarning` (E-198) and is still wrong.
- **`isTiltingAtVertex` is wrong and did not warn.** It said True on refused steps: 79 / 79 at n = 4, 214 / 260 at n = 5, 922 glued, 485 long-glue. It errs only toward True, by construction.

Neither error reaches a walk, by (a). Direct callers are exposed. I added the warning and did not change the value.

**Not claimed:** (H3) or H(I) = 0 for cyclic parents outside the scope; that the lemma covers H(I); Hom/End(T) or derived equivalence of the children; non-monomial cyclic parents.

**Refuted by:** an admitted step in scope whose Theta' or leads differ from its D(i) part's; H(I) != 0 on any of the 12 136; a finite-gl.dim cyclic parent on which (C) and (LF) fail at an admitted i.

## Evidence

| family | parents | admitted (A + A^op) | cycle reaches i | Theta' ≠ ∅ | all checks pass | D(i) audit | lib gate wrong | lib J wrong (on admitted) |
|---|---|---|---|---|---|---|---|---|
| sweep n=4, ≤5 arrows, rel ≤3 | 126 | 47 (47 + 0) | 36 | 0 | 47 | 47 | 0 (forced: rel ≤ 3) | 79 (0) |
| sweep n=4, ≤5 arrows, rel ≤4 | 137 | 47 (47 + 0) | 36 | 0 | 47 | 47 | **8** | 90 (0) |
| sweep n=5, ≤5 arrows, rel ≤3 | 474 | 342 (312 + 30) | 142 (A^op: 0) | 0 (Theta ≠ ∅: 28) | 342 | 342 | 0 (forced: rel ≤ 3) | 214 (0) |
| glued nL ≤ 4 (n ≤ 6) | 234 | 415 (367 + 48) | 156 | 32 | 415 | 415 | 0 (forced: glue rel ≤ 3) | 68 (0) |
| glued nL ≤ 5 (n ≤ 7) | 864 | 2 086 (1 774 + 312) | 735 | 294 | 2 086 | 2 086 | 0 (forced) | 250 (0) |
| glued nL ≤ 6 (n ≤ 8) | 3 204 | 9 891 (8 170 + 1 721) | 3 350 (1 192 through i; A^op: 0) | 1 873 (440 with cycle reaching i) | 9 891 | 9 891 | 0 (forced) | 922 (0) |
| long glue (one glue rel of length 4), nL ≤ 6 | 704 | 1 856 (1 658 + 198) | 815 (530 through i; A^op: 0) | 427 (133) | 1 856 | 1 856 | **175** | 485 (0) |

Totals: 47 + 342 + 9 891 + 1 856 = 12 136 admitted, 4 343 with a cycle reaching i, 0 failures. The glued rows are cumulative, and the long-glue row is disjoint from them. The "D(i) audit" column counts instances where the D(i)-sourced letters are closed under d and (H1w), Theta' and the leads computed on them equal those on all letters. Every glued admitted parent has a loop letter of degree >= 1 at a cycle vertex. No model hit the letter cap. Glued leads (nL <= 6, both sides): x2·x3 635, x3·x4 390, x4·x5 256, x2·r1 143, x3·r1 110, x2·x1 80, x1*·r3 150, x1*·r4 101, x2*·r4 74, and others (`glue6.out`).

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py seed            # 1 s: seeds + controls
.venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py bag             # 1 s: referee's b.a.g.x1 instance + [AUDIT] control
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py sweep 4 5 3   # 2 s  (sweep 4 5 4: 3 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py sweep 5 5 3   # 48 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py glue 6        # 160 s (glue 5: 30 s, glue 4: 5 s)
timeout 10m .venv/bin/python -W ignore workshop/rounds/127/toolsmith_cycgate.py glue 6 --long # 63 s
```

The outputs are saved as `toolsmith_cycgate_{seed,bag,4_5_3,4_5_4,5_5_3,glue5,glue6,glue6long}.out`, each under 5 KB. These are re-runs: the A-side tallies are unchanged from the first version. The script imports `126/theorist_cycparent.py` and `123/theorist_lnaparent.py` read-only.

## Prior record

- Loop at i impossible: E-279 (1).
- Walks never expand cyclic parents: E-211.
- Admitted steps with a 2-cycle upstream, step 7 right on them: E-199. The lemma (c) is the same mechanism applied to the model.
- E-130 flagged the simple-path gate as "uninvestigated".
- (H3) shapes: E-326, E-322.

New here:
- the T5 checks with the cycle upstream of an admitted i (E-326 had i on the cycle, never admitted);
- the D(i)-reduction of the model, proved from Def. 3.2 and audited;
- concrete wrong admissions by the library gate at finite gl.dim;
- `isTiltingAtVertex` wrong on cyclic input (my round-125 code; its tests were acyclic only).

## Code changed

- `quivermutation/mutation.py`: `isTiltingAtVertex` calls `arrowPaths.warnIfCyclic`, and its docstring says acyclic only and gives the n = 4 counterexample. The value is unchanged; the cost is about 3 %.
- `tests/test_cyclic_warning.py`: `test_isTiltingAtVertex_warns_on_cycle_upstream`. It pins the wrong True, marked "PINS A BUG". `pytest -q tests/test_cyclic_warning.py tests/test_search_soundness.py -m "not slow"`: 13 passed, 1 warning (pre-existing).
- `workshop/rounds/127/toolsmith_cycgate.py`: in revision, A^op in `glue()`, side tallies, `daudit`, `bag`, `--long`, `--only`, and the renamed tag.

## Next

- Skeptic/chair: should the gate and `isTiltingAtVertex` *refuse* (not warn) when a cycle reaches i? The exact version needs all nonzero paths (finite for finite-dim A). Regression set: the 8 n = 4 cases and the 175 long-glue cases.
- Theorist: H(I) = 0 with upstream letters is the one T5 check the lemma does not reduce.
- Toolsmith (me): non-monomial cyclic parents; 3-cycles at n >= 6.
