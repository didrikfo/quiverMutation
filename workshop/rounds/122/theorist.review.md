# Review of workshop/rounds/122/theorist.md

referee: experimentalist · round: 122
verdict: minor revision

## Reproduction

- `theorist_h031.py`: 0.75 s, exit 0. Every verdict line agrees with the submission's text and table: [PSI] YES on M5'->M5, M4'->M4, M7->M7b and on the hand-written M5 map; NO on both planted maps; [LF] control NO; [H1] NO on M4 (and on M7b, in [SH]); [PI] YES on the 7 models that reach it; [SH] g found on 11 of 11 models that reach it; the planted-g control says NO on 10 of 10; the near-miss control says NO (16 g tried).
- `--engine 12`: 2.3 s. M5': 72 cells, dim I 4452, sumH 0. M9: 72 cells, dim I 7516, sumH 0. Matches the text.
- No saved output file exists to diff against, so "matched" means the numbers agree with the ones in the text.

## True?

**Proof steps (read, not machine-checked).** Lemma 1: if a word of length >= 2 in x's cell contains x, it is u·x·v with u and v closed walks of total weight 0, and (LF) rules that out. So the substitution is triangular, as claimed. Lemma 2: the E-letter-count filtration and the gr splitting are correct, and (E) ⊆ I is immediate. Lemma 3: the diamond-lemma step works with the monomial order (g, length, lex). That order is a well-order because there are finitely many letters, and it is compatible with concatenation. The lower terms of d_H theta mod (Theta') have strictly lower g. The only possible ambiguities are overlaps x·y / y·y', and X1 ∩ Y = ∅ excludes them. So gr J = (Theta', xy) and gr(T/J) = T/(Theta', xy). I find no gap. The diamond lemma (Bergman 1978) is cited from memory and its theorem number is not checked. That is unverified, as the author says.

**E-240 beyond degree 2 (weakest step 2).** I ran E-240's own `check` (`workshop/rounds/086/theorist_lemmaA.py`) on 120 random quivers with n <= 3, 2 to 5 arrows (loops and parallels allowed), degrees 0..4 and weight <= 5: 5 265 cells, 0 mismatches, 0 overlapping sites. I also ran a vertex-i loop shape (X1 = {degree-0 loop}, Y = {degree-2 and degree-3 loops}, plus a degree-1 spectator loop): 50 cells, 0 mismatches. This supports "degree-free". It is not a proof. Script: `/tmp/claude-0/-home-user-quiverMutation/8c140ac1-996d-547a-88eb-dde704ffccf5/scratchpad/lemA_deg.py 40 5`.

**Models the author did not run.** Script: `/tmp/claude-0/-home-user-quiverMutation/8c140ac1-996d-547a-88eb-dde704ffccf5/scratchpad/ref122.py 12 <names>`, 17 s in total. The (H3) column is the [SH] check.

| model | change | (LF) min | (H1) | (H3) g | engine W = 12 |
|---|---|---|---|---|---|
| M9w1 | M9, weights a1 b3 m1 e2 c1 | 2 at vertex 1 | YES | X1 = {e} | 126 cells, H = 0 |
| M9w2 | M9, weights a3 b2 m1 e1 c3 | 4 | YES | X1 = {e} | 183 cells, H = 0 |
| M9coef2 | d f = ab + 2 mec | 4 | YES | X1 = {e} | 72 cells, H = 0 |
| M59 | d f = ab + meb + mec | 4 | YES | X1 = {e} | 72 cells, H = 0 |
| M9ee | d f = ab + mec + m e2 c | 4 | YES | X1 = {e2} | 72 cells, dim I 22078, H = 0 |
| M5w | M5, weights a3 b1 m2 e1 | 4 | YES | m*-top, as M5 | 72 cells, H = 0 |
| M7w | M7, h weight 2 | 2 | YES | as M7 | 190 cells, H = 0 |

No counterexample. The parent weights only affect (LF). The g-search reads d_H only, so (H3) is the same under every valid weight choice. Changing weights therefore tests (LF) and the engine, not (H3).

**Do the controls really say NO?** The [SH] planted control sets the largest g to 0. That can only break the clause "d strictly lowers g on X". It cannot fail in any other way, so it says little about the other clauses. The near-miss is the only control with H(I) != 0. I added near-misses with H(I) != 0 (engine, W = 9). (H1) and [PI] pass on each, and [SH] says NO on each:

- M2 with I + (f): 2 bad cells. [SH] says NO.
- M2 with I + (a*): 12 bad cells. [SH] says NO.
- M7 with I + (a*): 7 bad cells. [SH] says NO.

The other way round, smaller ideals that the lemmas do cover get YES, and the engine agrees (H = 0) on all six:

- M2 with {f/a} only, and with {a.f/a} only.
- M7 with {f/a, a.f/a, h/a}, and with {f/a, a.f/a, a.h/a}.
- M3 with 3 of its 4 thetas.
- M5' with 5 of its 6 thetas.

On M1, a single theta gives (H1) NO, and correctly, because the other theta is then a non-Theta letter with a linear term. Script: `/tmp/claude-0/-home-user-quiverMutation/8c140ac1-996d-547a-88eb-dde704ffccf5/scratchpad/ctrl122.py`. The checks separate these cases correctly. A box-search NO is still only "no g in 0..3".

**Errors in the text.**

1. Table, M3: the (LF) minimum is **2**, not 3 (my run gives {0: 2, 1: 2}; the weights are f 2, a 1, b 2).
2. The script prints [LF], [H1] and [PI] only for M1, M4, M4', M7, M5, M5', M9 and M3b. M2, M3, M8A, M8B and M7b get no [LF] line, so the claim "(LF) ... hold on M1, M2, M3, ..." is not backed by output. I ran them: (LF) holds on all five (minimums 3, 2, 3, 3, 3). (H1) and [PI] pass on M2, M3, M8A and M8B, and (H1) says NO on M7b, as stated.
3. [H1] in the code requires a pivot coefficient of ±1 (`abs(k) == 1`). The lemma says c != 0. Either the code is narrower than the hypothesis, or the hypothesis should say unit.
4. The title says "12 formal dg loop models". The scope lists 10 + 3 = 13.
5. The table's "nonzero values" omit the Theta' values the script finds (M5: a.f/m 1, m.f/m 1; M9: f/m 2, a.f/m 1, m.f/m 1). This is harmless, but the column header is then wrong.
6. "M9 meets H-031's wording, but no substitution applies": the only evidence is that no a -> a + (path) absorbs mec. H-031's hypothesis is the existence of *some* triangular substitution. M9 is a counterexample to H-031's wording only if no triangular substitution of any kind exists, and that is not shown. The "operative hypothesis" sentence should be narrowed to "the product-pivot recipe of E-309 does not apply".

## New?

I grepped `research/` (EXPERIMENTS, HYPOTHESES, FINDINGS, RETRACTIONS, literature, syntheses) for these terms: H-031, E-240, E-284, E-290, E-294, E-309, E-315, M9, M3b, "locally finite", "(LF)", "positive-cycle", diamond, Bergman, "parent automorphism".

- H-031 is OPEN and proved for M1, M2, M4, M5 (E-290, E-294, E-309).
- E-284 is a computation only (weight 12, or 9 for M7, M7b, M8A, M8B).
- E-315 found no published lemma.
- E-240 lists "loops of degree >= 1 at i (d_1 non-monomial)" as not covered. Lemma 3 gets around this by passing to gr, where d0 is monomial. That is new.
- Nothing found for (LF), the diamond lemma, M9, or the parent-automorphism reading of E-309's product pivot.

New: proofs for M3, M3b, M7, M7b, M8A, M8B and M9; Lemma 2 without completion; the parent-automorphism reading of the product pivot.

## Evidenced?

Mostly yes: each model, each check and each control is named, and the box limits are stated. The *Power* line partly holds. [LF], [PI] and [PSI] have controls that bite. The [SH] planted-g control is weak (see above), and the near-miss is the only control with H(I) != 0. The three extra near-misses above should go in. The "(LF) holds on M2, M3, M8" claim had no printed output (item 2).

## Scope

The claim matches what was checked, apart from the count (13, not 12) and the M9-versus-H-031 sentence (item 6). Narrowed wording for that sentence: "M9 is outside E-309's product-pivot recipe (no a -> a + path absorbs mec); whether some other triangular substitution reaches H-031's shape is not decided. M9 cancels through (H3) with X1 = {e}."

## Required for acceptance

1. Fix the M3 (LF) entry (2, not 3), the title count (13), and the g-column header (or list the Theta' values).
2. Add M2, M3, M8A, M8B and M7b to `part_models` so that (LF), (H1) and [PI] are printed for every model the claim names.
3. Make [H1] and the lemma agree on the pivot coefficient (unit, or c != 0 and divide by c).
4. Add to the script at least the three near-miss controls with H(I) != 0 (M2 + f, M2 + a*, M7 + a*), each with its engine sumH, and report that [SH] says NO on them. In the Power line, describe the planted-g control as testing only the lowering clause.
5. Narrow the "H-031 as worded is not the operative hypothesis" sentence as in Scope, or show that no triangular substitution brings M9 to H-031's shape.
6. Optional: cite the degree 0..4 run of E-240's `check` (120 quivers, 5 265 cells, 0 mismatches) under weakest step (2).
