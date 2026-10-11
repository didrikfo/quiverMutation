# Review of workshop/rounds/026/toolsmith.md

referee: experimentalist · round: 026
verdict: minor revision

## Reproduction

- `.venv/bin/python -m pytest -q tests/test_longsquare.py`: 3 passed (1 s).
- Resume check, n = 7 class 0 (`--max-exp 3000` once, 61 s, exit 2; then `--max-exp 1000/2000/3000` with one `--ckpt`, exit 2 each): the result tables are identical (diff empty after dropping the first lines, which differ only in "STOPPED" vs "CHECKPOINT/RESUME" and times). Final line `level 7 pos 1696 of 1715 expanded 3000` is the same in both. Exit-2 on budget also confirmed: `--budget-sec 5 --ckpt` exits 2.
- n = 8 class 1, `--budget-sec 420 --show 0` (about 7 min, my run: 24 872 algebras, 9 490 expansions, so less load/coverage than the author's 10 106): rows `new True/old False` = 4, `new True/old True` = 12 (both match the claim exactly); out-degree >= 2 J != 0 = 29 (author 38; the author flags counts as load/cap dependent). Distinct (parent, v): 45 = 29 D-part1/2 + 8 (dim J 2) + 4 sq + 4 sq-parallel-only. The 4 parallel-only steps reproduce.

## True?

No error found. Independent checks (own scripts, scratchpad):
- Hand-built, different from the author's shapes: `2 => 5 -> 4` with `g0 e = g1 e` (no prefix arrow): old False, new True. Three-path relation (`g0`, `g1` through the doubled 2 => 5 and one through 3, all ending `5 -> 4`; coefficients 1, 1, -2): old False, new True. v with two (doubled) out-arrows: False. Relation ending in a later arrow `e f` with paths differing before `e`: False at both v = 5 and v = 4 (correct: penultimate arrows coincide).
- Searching for a positive off the author's walk: LNAs and duals, two reduced mutation steps, all vertices: n = 6, 4 452 tests; n = 7 (stopped at 400 s, partial list of LNAs), 22 344 tests: 0 positives, 0 old/new differences. So the only real positives are the n = 8 class 1 walk steps; the old-vs-new unit sweep is vacuous as the author says, and a two-step sweep does not fix that. I found no counter-case.
- Unchecked by the author, noted: the test requires `len({second-to-last arrows}) == len(paths)`, so a relation with 3 terms where two share a penultimate arrow is False even if a sub-pair is a square. Matches the old behaviour and the stated meaning, but it is a definition choice. Also `procedure.relationsFrom` with no `arrowRels` on a parallel quiver raises (stated, by design).
- "Distinct (parent, v)" counts and the claim that old/new `longsq` agree "on all rows except the 4": confirmed on my table.

## New?

Nothing found that does this: `longSquare` with arrow relations appears nowhere in `quivermutation/`. E-108 states the defect (out-degree 1 no-square rejects = parallel-arrow squares); E-105 states the minimality caveat; E-107 the n = 8 c0 out-degree 2 rejects. The 4-of-16 count and the checkpointed walk are new. Checkpointed walk predecessor: round 019 `toolsmith_walk.py` (cited by the author, different counts).

## Evidenced?

Mostly. Specific: sizes, times, expansion counts, resume diff. Gaps: (a) the "16 steps" claim and the 4 are on one walk, cap and load dependent (mine 4 and 12 as well, same split); (b) the 4 parallel-only steps are not shown to have a minimal relation (E-105 caveat), the author leaves this to the next round, so "long square" there means a relation present in the presentation, not a verified minimal one; (c) the resume check is n = 7 only and at a deterministic `--max-exp` cut, not SIGTERM or a time slice, and not a crash/restore from a periodic save; (d) the n = 9 sizing is an extrapolation (ratio 2.5, not flattening), stated as such. The overnight command was not run and the author does not claim it.

## Required for acceptance

1. State in the claim that the 4 parallel-only steps are relations present in the presentation, unverified for minimality (E-105), so "sees parallel-arrow long squares" does not imply "genuine long relation".
2. Say that the n = 6, 7 old-vs-new sweep is vacuous in the headline evidence, not only in the body, and (optionally) that a two-step sweep at n = 6 and a partial n = 7 also give 0 positives.
3. Note the "all penultimate arrows distinct" requirement (3 terms, 2 sharing, is False) as a definition choice.
4. Optional: a SIGTERM or time-slice resume check, since the overnight job relies on those, not `--max-exp`.
