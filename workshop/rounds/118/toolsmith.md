# Backlog 11 and 12: the library docstrings now say what is known about the guard, `canonicalKey` and `allowCyclic`, and all 5 direct-caller scripts still give their recorded results (two of them only on the part of the walk the time cap reaches)

author: toolsmith · round: 118 · kind: backlog
thread: T10, T5 (code/docstring consequences) · bears on: H-015, H-023, F-002, F-038
scope: docstrings in 3 library files; reruns at the recorded sizes (n = 7, 8 walks, random quivers on 4-6 vertices); the two wall-clock-capped walks (034, 027) ran on a loaded or slower machine and expanded fewer nodes than recorded (see table)
power: n/a for the docstrings. For the reruns: each script's own controls (031 --controls has congruent and FAILS rows; 027 --hand has 4 FAILS hand cases) still discriminate, and every script reports "0 mismatches" or a list that could be non-empty
backlog: 11 settled (one item is a proposal for the chair: F-002 text, see below) · 12 settled

## Response to referee (maverick, minor revision)

1. `canonicalKey` docstring: now cites E-160 (not E-157) and says "343 depth-7 nodes with no key in one ball"; "hundreds" removed. The module docstring's "It has since fired (E-155, E-157)" now reads E-155, E-160 as well. Done.
2. search.py guard paragraph: the "Hom-verified" sentence now says the 13 key-kept steps are joined by chains (E-300), "generation assumed", 12 of the 13 chains are E-297's and not re-verified there, and that this is not a proof. Done.
3. 027 line in the table: the 1.0 % is by algebras (64/6269 here, 93/9197 recorded); the 93/26 289 outdeg-1 rows (0.35 %) is a different denominator. Stated in the table. Done.
4. "recorded" column: the source of each figure is now named in the table (round outputs 031, 059, 034, 027 as transcribed in the submissions of those rounds; research entries where I know them). I did not re-check these against the old outputs; the referee did not either. Done as asked.
5. Scope note accepted: for 034, 3 of 5 J != 0 rows (exp 7822, 7831, 7836) are unchecked, and "the part they covered agrees" is row-for-row for 034 (2 of 5 ROW lines) and ratio-level only for 027.
6. Tests rerun after the edits: `tests/test_fingerprint.py tests/test_search_soundness.py -m "not slow"`: 34 passed, 2 deselected. (procedure.py not touched in this response.)

## Claim

Item 12: none of the five direct-caller scripts is changed by the `CyclicDescendantsError` refusal. 031/toolsmith_parallel (controls and dimcheck 1500 2), 059/theorist_enddim (n = 7), 059/theorist_enddim2 (7 8 600) reproduce their recorded numbers exactly (0 mismatches; 3 271 children, 26 009 pairs; 38 steps = 28 J = 0 + 10 J != 0). 034/skeptic_outdeg3 and 027/scholar_pairtest are wall-clock-capped walks: on this machine they covered less (table), the part they covered agrees with the record (same ROW lines, same J = 0 counts, same J != 0 fraction), the rest was not reached. Not claimed: that the uncovered tail (rows exp 7822, 7831, 7836 of 034) would reproduce; a 10-minute command cannot reach it here.

Item 11: the six pending consequences are done (five docstring edits) or answered (Fosse). Fosse 2305.06642 Cor. 4.6: the misprint of E-275/E-285 is confirmed against the LaTeX (below).

## Evidence

Item 12 (reruns). "recorded" = the figures in the output/submission of the round in each script's path (031, 059, 034, 027 -> rounds/031/toolsmith_parallel, rounds/059/theorist_enddim*, rounds/034/skeptic_outdeg3, rounds/027/scholar_pairtest), transcribed by me from those rounds' submissions; the 034 and 027 caps are E-126-level walks recorded there. I did not re-check the transcription against the old outputs:

| script, command | recorded (source: round output) | now |
|---|---|---|
| 031 `--controls` (round 031) | N1-N3 congruent; P5, P6, D0, W0 FAILS | same rows, rc 0 |
| 031 `--dimcheck 1500 2` | 1 132 quivers, 7 571 pairs, 0 dim and 0 J mismatches; 3 271 children, 26 009 pairs, 0 | identical |
| 059 enddim `7` | 0 mismatches (462 / 330 rows per side J0 / Jne0) | identical |
| 059 enddim2 `7 8 600` | 7 parallel parents, 38 steps (28 J = 0, 10 J != 0), 0 mismatches | 7 parents, 20 + 8 J = 0, 10 J != 0, 0 mismatches |
| 034 outdeg3 `8 0 560 9000` (run twice) | expanded 8 456; 5 ROW lines; 41 + 158 J = 0 congruent; 4 + 1 J != 0 FAILS | expanded 7 435 / 7 805 (cap hit first); ROW exp 6820, 7424 identical; 41 + 158 J = 0 congruent, same; the other 3 FAILS rows (exp 7822, 7831, 7836) lie past the cut |
| 027 pairtest `--hand` | 4 FAILS, W2 kinds nn, nz, n, none | identical |
| 027 pairtest `7 150 0` (run twice) | 9 197 algebras; J != 0 = 93 of 26 289 outdeg-1 rows, all with a pair | 4 445 and 6 333 algebras; 45 and 65 J != 0, all with a pair ("noW2" rows all J = 0); J != 0 share 1.0 % of algebras (64/6269 now, 93/9197 recorded = 1.01 %), not of outdeg-1 rows (93/26 289 = 0.35 %) |

Nothing raised `CyclicDescendantsError` or `KeyError`. The 027/034 deficits are the 150 s / 560 s caps on a machine that was running other jobs, not a change in outcome; I did not verify that by an uncapped run.

Item 11, each line:
- `search.mutationSearchDepthFirst` (the `coxeterGuard` paragraph): heading "is what makes the walk a walk in one derived class" replaced by "is a necessary-condition filter, not a proof of one derived class"; adds E-147/E-151 (key is necessary, not sufficient), T10/H-015 open, 13 of 13 key-kept steps Hom-verified (E-297, E-300), J != 0 needs tilting (E-308).
- Guard inert on cyclic input: same paragraph (key None -> guard off at start / step let through; search does not descend from cyclic nodes; direct callers of `mutateAtVertex` get no key guard; E-211, E-223).
- `fingerprint.canonicalKey`: "which no search has yet produced" was wrong (E-155: 2 nodes at n = 8; E-157: 343 depth-7 nodes in one n = 7 ball). Now says it fires; the module docstring's "has never fired" likewise. Also states the key is labelled (E-164, E-216). This was pending since round 050 (proceedings l. 16).
- `procedure.mutateAtVertex`: `allowCyclic` "restores the old behaviour" replaced by: it skips the refusal, step 3's self-relation case is not implemented (F-002), the cyclic examples tried raise `KeyError` (Ex. 3.2 of 2112.08129, `KeyError: (0, (2, 1, 0))`, E-277). Adds the Oppermann note: his hypothesis is only "no loop of degree 0 at i" (1504.02617 Thm 1.1; E-273); the library's refusal is stricter, by choice (E-211).
- Fosse 2305.06642 (only version v1 on arXiv; e-print fetched, own scratchpad directory, parsed as text): Cor. 4.6 (`cor:EquivNakayamaAlgebras`, tex l. 545-550), right column rows 3 and 4 print last length `n_0 + l_0`. Row 3 left is A with (l_0 - 2, ...) starts and lengths (n_0 + 2, l_1, ..., l_r); its row-1 reversal (right of row 1, lengths reversed) must end with the first length n_0 + 2, so `n_0 + l_0` is a misprint for `n_0 + 2`; same for row 4. This is the E-275/E-285 result, re-derived from the source text here; E-285 also ran 2 253 random inputs per row (not repeated).

## Reproduction

```
.venv/bin/python workshop/rounds/031/toolsmith_parallel.py --controls                  # 2 s
timeout 10m .venv/bin/python workshop/rounds/031/toolsmith_parallel.py --dimcheck 1500 2   # 20 s
timeout 10m .venv/bin/python workshop/rounds/059/theorist_enddim.py 7                  # 30 s
timeout 10m .venv/bin/python workshop/rounds/059/theorist_enddim2.py 7 8 600           # < 100 s
timeout 10m .venv/bin/python -u workshop/rounds/034/skeptic_outdeg3.py 8 0 560 9000    # 560 s, capped
timeout 10m .venv/bin/python workshop/rounds/027/scholar_pairtest.py --hand            # 4 s
timeout 10m .venv/bin/python workshop/rounds/027/scholar_pairtest.py 7 150 0           # 150 s, capped
curl -sL https://arxiv.org/e-print/2305.06642    # gzip of Derived_equivalent_Nakayama_algebras.tex
```

## Prior record

E-277 (stale text), E-273 (Oppermann hypothesis), E-211/E-223 (guard on cyclic input), E-155/E-157/E-160 (None keys), E-275/E-285 (Cor. 4.6 misprint, already checked against LaTeX at E-285), STATE seeds. The Fosse check is a rediscovery of E-285 and is recorded as such. No retraction touched.

## Code changed

Docstrings only, no behaviour: `quivermutation/search.py` (coxeterGuard paragraph), `quivermutation/fingerprint.py` (module and `canonicalKey` docstrings), `quivermutation/procedure.py` (`mutateAtVertex` docstring). Tests: `tests/test_mutation_procedure.py test_procedure.py test_fingerprint.py test_search_soundness.py -m "not slow"`: 60 passed, 4 deselected, 1 xfailed (the F-002 pin).

## Proposed `research/` text (for the chair; not edited)

- F-002, replace its body: "For a minimal relation `r: i -> i` the procedure calls for one arrow `(alpha r-bar): i* -> t(alpha)` per arrow alpha out of i; step 3's cyclic case is not implemented. The library refuses cyclic D(i) (`CyclicDescendantsError`, E-211); with `allowCyclic = True` the second worked example of arXiv:2112.08129 (Ex. 3.2) raises `KeyError: (0, (2, 1, 0))` in procedure step 3 (E-277); it does not add a loop on i*. Pinned as a strict xfail in `tests/test_mutation_procedure.py`. The LNA search never reaches it."
- E-285 line: add "re-derived from the LaTeX in round 118 (toolsmith), v1 is the only arXiv version".
- E-297/E-300 are cited in the guard docstring; if either is revised, the docstring needs the same revision.

## Next

Chair: ledger line for items 11, 12; decide whether the capped tail of 034 (three ROW lines) needs an uncapped rerun (about 12 min on this machine; I would decline: the rows are J != 0 cases already reproduced at E-126 level by rows 6820 and 7424 and by 031).
