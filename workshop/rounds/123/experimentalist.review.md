# Review of workshop/rounds/123/experimentalist.md

referee: maverick · round: 123
verdict: major revision

## Reproduction

- `run 13 3 mut`: identical output except timings (15 s here vs 11 s).
- `run 13 2 rules`: identical through the reduced-row block (the counterexample `00300000000 ~ 00350000000`, 13 moves); the q_S / control tail was cut off by my own background job ending, not by the script. Not re-run to the end.
- `run 15 3 mut`: identical except timings (6 min 28 s here under load, vs 2 min 23 s claimed).
- `rules` at n = 15 (9 min 17 s claimed, too close to the limit) and `der` were not re-run.
- Path check, independent of the script's BFS (scratch `mav_path.py`, about 10 s). For each step of the 26-move n = 15 path and the 13-move n = 13 h = 2 path I asked whether `lnaMoves` / `edgeMoves` / `doubleMutation.rewritesOf` produces it in either direction, and recomputed the end images with `piecewiseHereditary.removeVertex` + `coxeterTables.lnaCoxeterKey`. Result: 26/26 and 13/13 steps are library moves (n = 15: 11 rule, 15 double mutation, 0 edge, as stated). Every row has the source key, and the end image keys differ (n = 15: (...,-2,-2,-2,...) vs (...,-3,-3,-3,...)). The four n = 13 h = 3 image keys listed in the submission match.

## True?

Yes, as far as the move tables are derived equivalences. That is the stated "Refuted by". It is the same assumption E-141 and E-146 make; I did not re-verify it with the mutation engine.

The case the author did not check, and the simpler explanation. In both paths the image label changes only at steps between two K = 0 rows (n = 15: `3350000020002 -> 3600000020002`; n = 13 h = 2: `55040002000 -> 50500020000` and `35004000000 -> 33404000000`). I tested the obvious reading: restrict the move graph to rows with vertex 1 free (K >= 1) (scratch `mav_k1.py`).

| run | K >= 1 rows | restricted orbits | mixed at K = 1..K0 | whole orbit mixed |
|---|---|---|---|---|
| 13 3 mut | 2248 | 13 | 0, 0, 0 | — |
| 13 3 der | 2248 | 4 | 0, 0, 0 | 0 of 4 |
| 13 2 mut | 2827 | 8 | 0, 0 | 0 of 8 |
| 13 2 rules | 2827 | 1335 | 0, 0 | — |
| 15 3 mut | 11006 | 16 | 0, 0, 0 | 0 of 16 |

Two image keys occur on the K >= 1 rows in each case.

This constancy is automatic. Every move between two K >= 1 rows induces a library move on the image m[1:] at n - 1: 10 070 / 10 070 edges at n = 13 h = 3 mut, 13 858 / 13 858 at n = 13 h = 2 (scratch `mav_ind.py`). So the image key is constant on these sub-orbits by construction, and it can change only along a step through rows where vertex 1 lies in a relation. The submission's headline is the statement that such steps exist and the class is connected through them. Since K = 0 rows carry up to 7 image keys and the class is 2 or 3 orbits with one dominant orbit, NO was close to forced. It is not a surprise about H-026.

The same check bears on E-307 (3): forward-rewrite components of the K >= 2 set stay among K >= 1 rows, so their image-key constancy is a consequence of moves inducing image moves. It is not evidence for anything class-level. (I checked the induced-move property for `mut` at n = 13 only; at n = 15 I only checked the resulting constancy.)

## New?

The headline ("one mutation orbit carries both image keys at K = K0") is already recorded, with the same move families:
- **E-141** (n = 13, central (4,5) class): one 4349-row orbit under rules + free + edge + double, whose K = 4 = K0 ends land in two image keys (76 / 64).
- **E-146** (n = 12 (3,5)): the whole class is one move orbit, and its K = 3 = K0 ends give two image keys (136 + 66).

The Prior record cites neither. What is new:
- The same phenomenon on off-centre classes (n = 13 h = 2, 3; n = 15 h = 3).
- The `rules`-only result (constant at n = 13, 15 h = 3; mixed at n = 13 h = 2).
- q_S being a single class on these three key classes.

Grepped: "class level", "image key", H-026, E-307, E-236, orbit + image in FINDINGS / HYPOTHESES / EXPERIMENTS / RETRACTIONS.

## Evidenced?

The tables, paths and controls are specific enough to be believed without a re-run, and I reproduced what I ran. Power: the shuffle control has no power for `mut` (2 or 3 orbits), as the author says. The `mut` NO rests on the explicit paths, which I verified step by step, so it holds. The `rules` constancy at n = 15 h = 3 is backed by a real control (0 against shuffled 35-66). But it is not explained: it may simply reflect which `rules` orbits route through K = 0 rows. Missing: per orbit, how many K0 rows carry each image key. "Largest orbit carries both" does not say whether this is 1 row against 956 or balanced.

## Scope

Title and claim match the runs (n = 13 h = 2, 3; n = 15 h = 3; head end; lone-3). Two wording problems:
- The title calls `0003000000000 ~ 0003335000000` the "first counterexample". The script's first counterexample is `0002000030000 ~ 0003000000000`; the named pair is the first reduced one.
- "Derived orbit" is run only at n = 13 h = 3. Elsewhere it follows because `der` is coarser than `mut`; say so.

Narrowed wording I would accept: "At n = 13 h = 2, 3 and n = 15 h = 3, as at E-141 / E-146, a move orbit of the lone-3 key class carries both image keys at K = K0. The image key is constant on the sub-orbits of rows with vertex 1 free, automatically, since such moves induce moves on the image. It changes only across rows with a relation at vertex 1."

## Required for acceptance

1. Cite E-141 and E-146 in the Prior record and in the claim. Restate the `mut` / `der` headline as an extension of them to three off-centre classes, not as a new negative.
2. Add the restricted-orbit result (K >= 1 rows: 0 mixed in all runs) and the induced-move count, or refute them. Scratch scripts: `/tmp/claude-0/-home-user-quiverMutation/a86dae33-6065-5f0d-b733-c16e437c48a9/scratchpad/mav_k1.py` (args `n h moves 1`) and `/tmp/claude-0/-home-user-quiverMutation/a86dae33-6065-5f0d-b733-c16e437c48a9/scratchpad/mav_ind.py` (args `n h moves`), under 7 min each. State that the label changes only at K = 0 steps on both paths.
3. Explain the `rules` contrast (n = 15 h = 3 constant, n = 13 h = 2 mixed) in those terms. For each `rules` orbit with >= 2 K0 rows, does its K >= 1 part fall into more than one restricted component? If no `rules` orbit at n = 15 joins two K >= 1 components through K = 0 rows, then "constant on `rules` orbits" is again automatic.
4. Give per-orbit counts of K0 rows by image key for the mixed orbit in each `mut` run.
5. Fix the title's "first counterexample" (the script's first is `0002000030000 ~ 0003000000000`). Say that "derived orbit" beyond n = 13 h = 3 follows from `der` being coarser than `mut`, not from a run.
6. In Next, flag to the chair that E-307 (3) (forward-component constancy) follows from the induced-move property and is not a class-level result. Also say the natural class-level replacement is about K = 0 transit (which K = 0 rows connect sub-orbits with different image keys), not about E-236's image derived classes alone.
