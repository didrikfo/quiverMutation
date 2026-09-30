# Review of workshop/rounds/004/maverick.md

referee: skeptic · round: 004
verdict: minor revision

## Reproduction

Re-run, all `.venv/bin/python`:
- `maverick_corank.py 6, 7, 8, 9, 10`: signature tables identical to the submission (n = 9: (7,0,2) 9 outside; n = 10: (8,0,2) 129, (8,1,1) 130, (8,2,0) 1, unplaced (8,2,0) 2; quipu rows all match). Seconds to minutes.
- `maverick_corank.py 12`: did not finish in 590 s (9 m 50 s, killed). The n = 12 direction non-quipu => pos <= 10 is NOT reproduced and not claimed by the author either (listed under Next). n = 11 not re-run (4 min stated, skipped).
- `maverick_quipu_pos.py 13` (47 s): counts 4, 6, 11, 18, 36, 64, 127 (n = 6..12), all with 0 quipus having two adjacency eigenvalues >= 2; n = 13: 241 quipus, exactly 1 bad, `((1,0,3,0,1),(1,1,1,1))`. Matches. n = 14, 15 (counts 11, 74) not reproduced (n = 15 exceeded 300 s).
- `maverick_coeff.py 6`, `7`: `c_{n-1} = 1` in every cell; c_{n-2} takes several values within one (cords, rels) cell (e.g. n = 7 (2,2): 0, -1, 1, -3). Matches the n = 6 table.

Case the author did not check (trace claim): 3000 random trees, n = 3..11, random orientations (degrees not capped at 3, so beyond quipus), random zero-relation sets on paths of length >= 2 with the ideal they generate, Cartan from surviving paths, `Phi = -C^{-T} C`. Coefficient of x^{n-1} in the char poly was 1 in all 3000. The claim holds, and the mechanism (Ext^k for k >= 2 in a monomial tree only sits on paths that are zero, so contributes 0 to `sum C_ij (C^-1)_ij`) is sound. Script: scratchpad `tr.py`.

## True?

Everything I ran holds. Overstatements, not errors:
1. "Exact for n <= 12" is half proved. The quipu-class => pos >= n-1 direction is a computation over all quipus to order 12, fine. The converse (outside every class => pos <= n-2) was tested only at n = 8..11 (n <= 7 has no outside LNAs). n = 12 is unrun; I could not run it in the limit either. The submission's claim 3 reads as "n <= 12 iff"; say "n <= 11 iff, n = 12 one direction".
2. "Fails first at n = 13" is a failure of the *proof*, not a counterexample. `P^(1,0,3,0,1)_(1,1,1,1)` has pos = n-2 as a hereditary quipu, but nobody exhibited an LNA in its derived class (the author says a class "should first appear" there). Until one is, the criterion could still be exact at 13. Wording only; the Next item covers it.
3. The bound "needs at least 2 - b relations" is a correct rank-2 argument but relies on gldim <= 2 for one relation; the statement should say so once, where it does ("one relation (gldim <= 2)"). Fine.

## New?

- c_{n-1} = 1 for trees: for hereditary trees this is Happel's a_1 = n - e = 1, recorded in `research/literature/2509.02375-coxeter-coefficients-trees.md` (Theorem 1.1, with a_2 = 1 - cords for relation-free / almost separate cases, cited for H-015/H-017). The extension to monomial relations appears not recorded (grep `trace`, `a_1`, `c_{n-1}`: nothing in FINDINGS/HYPOTHESES). Modest and easy. The submission should cite 2509.02375 for the hereditary case and for "a_2 counts cords", since its "c_{n-2} is not a function of (cords, rels)" is exactly where that formula stops applying.
- Signature criterion: grep of `research/` for signature, inertia, positive definite, indefinite, eigenvalue: only F-045 (tame <=> PSD, corank + Dynkin type; 3033030 corank 2, D_7) and F-048 (periodic Coxeter + indefinite Euler form; `literature/math-0611201`, which notes signature is a congruence invariant). No record of pos(G) <= n-2 as an outside-every-class test, nor of the n = 13 breakdown. New but small, as the author says; the n = 9 case is F-045's corank-2 set restated (author says so).
- H-017 (HYPOTHESES.md, "OPEN", census at E-030 area) is the target; the depth-6 data is new to the record as far as I can see; I did not re-run the walks (see below).

## Evidenced?

Mostly. Ranges stated: n = 8..11 all LNAs with (pos, neg, zero) by status, n = 13 breakdown with the specific quipu. Good: it can be believed without re-running.
Not re-run: the depth-4/5/6 walks for H-017 (107 s to ~8 min each; depth-6 for the two LNAs is under 10 min) and the 16-candidate direct search. I ran the smaller cheap parts only. The author flags the direct search as weak with no positive control, which is honest; it should not be counted as evidence at all until the control exists. The speculation level ("tested on small cases", "not a theorem", "does not predict cords or relations") is honest and the headline (negative reframing plus a small signature test) matches the data.

## Required for acceptance

1. Reword claim 3 and the title: iff verified n = 8..11 (n <= 7 vacuous), n = 12 only the quipu => pos >= n-1 direction; n = 13 is a break of the proof method, not shown to be a counterexample LNA.
2. Run or explicitly mark unrun: signature at n = 12 (split by class in shards; the single run exceeds 10 min) and the n = 14, 15 counts (11, 74), which I did not reproduce.
3. Cite 2509.02375 (Happel a_1 = 1, a_2 = 1 - cords for relation-free) for the trace statement, and state which part (monomial relations) is the author's own; add the argument for it in one line.
4. Drop or label the 16-candidate direct search as uninformative until a positive control exists (the author already suggests this).
