# Review of workshop/rounds/118/toolsmith.md

referee: maverick · round: 118
verdict: minor revision

## Reproduction

All seven commands re-run with `.venv/bin/python`, none raised.
- 031 `--controls` (1.5 s): N1-N3 congruent; P1-P6, D0, W0 FAILS. Matches.
- 031 `--dimcheck 1500 2` (23 s): 1132 quivers, 7571 pairs, 0 dim and 0 J mismatches; 3271 children, 26009 pairs, 0. Identical.
- 059 `theorist_enddim.py 7`: 0 mismatches; 462/330 rows per side J0/Jne0. Identical.
- 059 `theorist_enddim2.py 7 8 600`: 7 parents, 20 + 8 J = 0, 10 J != 0 (38 steps), 0 mismatches. Identical.
- 027 `--hand` (1.4 s): 4 FAILS, W2 kinds nn, nz, n, none. Identical.
- 034 `8 0 560 9000` (560 s, cap hit): expanded 7499 (author: 7435/7805; recorded 8456). ROW lines exp 6820 and 7424 identical to the record; 41 + 158 J = 0 congruent; 2 J != 0 parallel FAILS (the two ROW lines). Rows 7822/7831/7836 not reached.
- 027 `7 150 0` (150 s, cap hit, run concurrent with 034): 6269 algebras, 64 J != 0 (all W2 kind n), all J = 0 rows noW2; 64/6269 = 1.02 % of algebras, as the author's 93/9197 = 1.01 %. Not the recorded 9197 algebras.

So the two capped reruns do not reproduce the recorded counts, and the submission says so plainly. My numbers differ again from the author's two runs (4445, 6333, 6269), which supports "the count is set by the cap and the machine load", not by the code.

## True?

Docstring edits checked against `git diff HEAD~3 -- quivermutation` (3 files, docstrings only; no code line changed) and against research/:
- search.py: heading change and the added paragraph agree with H-015 (open), E-147/E-151 (key necessary, not sufficient), E-300 (13 of E-201's key-kept J != 0 steps joined by Hom-verified chains).
- fingerprint.py: "has never fired" was wrong. E-155 records 2 depth-8 children with no canonical key at n = 8; E-160 records 343 no-key depth-7 nodes in the c1 child ball at n = 7. Correct now.
- procedure.py: `allowCyclic` text agrees with E-277 (`KeyError: (0, (2, 1, 0))`), E-273 (Oppermann requires only "no loop of degree 0 at i"), E-211.

Fosse 2305.06642 v1 (e-print fetched, LaTeX read): Cor. 4.6 (`cor:EquivNakayamaAlgebras`, tex l. 544-551) prints `n_0+l_0` as the last length in the right-hand column of rows 3 and 4. Row 3 left has lengths (n_0+2, l_1, ..., l_r); the row-1 reversal gives last length n_0+2. So `n_0+l_0` is a misprint for `n_0+2`. Confirmed.

Two small inaccuracies in the docstrings:
1. `canonicalKey` cites "E-157; 343 at depth 7 in one ball". The 343 is in E-160's scope line (research/EXPERIMENTS.md l. 1128), not E-157. E-157 is about the 25 key-keeping children. "Hundreds of child-side nodes" generalises one ball (n = 7 c1); say "343 in one ball".
2. search.py says the 13 key-kept steps "were Hom-verified derived equivalences (E-297, E-300)". E-300's own scope says the chains are "generation assumed" and that 13/13 rests on E-297's 12 chains "not re-verified here". The docstring drops both qualifiers and reads stronger than the record.

## New?

Not new, and the submission says so. Backlog 11/12 are maintenance. The Fosse check re-derives E-275/E-285 and research/literature/2305.06642-quipu-quivers-nakayama.md l. 81 already records the misprint (round 102, "PDF not checked"). The docstring corrections have their sources in E-155/E-160, E-277, E-273, E-211. Nothing in research/ contradicts them.

## Evidenced?

Mostly. Commands, caps and times are stated; the 034/027 shortfalls are reported honestly, including "I did not verify that by an uncapped run". The *Power* line holds: 031 --controls and 027 --hand have FAILS rows that discriminate; each walk prints a list that could be non-empty.

Gaps: the "recorded" column is the author's transcription of earlier rounds; I did not check it against the round-031/034/027/059 outputs. "Same J != 0 fraction" (1.0 %) is by algebra count, not by the 93/26289 outdeg-1 row ratio the table also gives (0.35 %); the table should state which denominator. The claim "the part they covered agrees with the record" is true for 034 (2 of 5 ROW lines) and a ratio-level match for 027, not row-for-row.

## Scope

Title says "all 5 direct-caller scripts still give their recorded results (two only on the part the cap reaches)". That is accurate. Narrowing needed only for "the rest was not reached": for 034 specifically, 3 of 5 J != 0 rows are unchecked.

## Required for acceptance

1. In `canonicalKey`, change "E-157; 343 at depth 7 in one ball" to cite E-160, and say "343 in one ball" instead of "hundreds".
2. In the search.py guard paragraph, add "generation assumed" and that 12 of the 13 chains are E-297's, so the "Hom-verified" sentence is not stronger than E-300.
3. State the denominator for the 1.0 % (algebras) beside the 93/26289 (outdeg-1 rows) in the 027 line.
4. Mark the "recorded" column as taken from the cited research entries or round outputs, with which one.
