# Review of workshop/rounds/130/toolsmith.md

referee: maverick · round: 130
verdict: minor revision

## Reproduction

- `pytest -q tests/test_search_soundness.py -m "not slow"`: 7 passed, 6.3 s. Matches.
- `toolsmith_f041refuse.py control`: off 1 meeting (child key met), on 0 meetings, counter 1; ball(P, 1) on 7 / off 8 keys, 1 removed (the child). Matches, 2 s.
- `run 10 0 0 --direct --idx 0,303`: pair 0 (22230300 / 02230300) keys 59 / 87, pair 303 (00030302 / 00030300) keys 156 / 219, 0 refusals, 1 meeting both ways, 0 removed / gained. Same keys and meetings as the TSV rows (pair 303 was `inferred` there; direct now). 20 s.
- `summary 10` and `summary 11` from the committed TSVs: 206 / 206 and 924 / 924 joined, classes 113 / 113 (204 orbits) and 288 / 288 (609 orbits), 0 refusals, 0 removed. Matches the table.
- Not re-run: the full n = 10 / 11 sweeps (25 min and 75 min wall).

## True?

I found no error. Is the detector actually reached? The counter in the script counts only `False` returns, so "0 refusals" alone cannot tell "never False" apart from "never called". I wrapped the script's own counting wrapper to count calls too (scratchpad `mav.py`, imports the round's module, so the same path runs). n = 10 pairs 0, 37, 75, 150, 225, 303: 321, 181, 421, 698, 460, 943 calls. n = 11 pairs 131, 761, 1512: 576, 732, 1454 calls. 5 786 calls in all, 0 False. So the test does run at every gate-admitted step. The monkeypatch reaches it because `search.py` imports `mutation` as a module and looks the name up at call time (l. 462).

The simpler explanation, which the submission does not give: in these balls the option *cannot* fire. If dim e_i A e_v <= 1 for every i, Ladkani's injectivity condition on that one class [p] says only that p beta is not in I for some beta. That is the gate's own "does not die against every out-arrow" condition. So a gate-admitted step at a vertex v with dim e_i A e_v <= 1 for all i is tilting. This argument assumes the gate implements its docstring. My counts agree with it. Of the 5 786 sampled steps, 765 are at vertices with >= 2 out-arrows, and 0 have any dim e_i A e_v >= 2. E-318's refused step (vertex 3 of `_e318Parent`) has out-degree 1 and dim e_i A e_v = 2. Its six J = 0 siblings all have dim <= 1. So the planted control fires on exactly the feature these LNA balls lack. The control proves the wiring works. It does not show that this population contains any step where the answer could have been NO. That is fine for a negative, but the *Power* line should say it. The general point: a J != 0 step requires a vertex v and an i with two independent paths i ~> v modulo I. On that reading, the null result says depth-3 F-only balls around LNAs at n <= 11 stay thin at every admitted step. That is a cleaner statement than "the option costs nothing", and it says where to look for a cost: E-332 at depth 6, and E-316's mixed walks.

## New?

- n = 10 was already recorded. E-321 (R123) found no key-kept J != 0 edge in the depth-3 balls of every F-041 pair at n = 8, 9, 10. E-323 part (6) (R125) found the option prunes no F-041 merge at n = 8, 9. The submission says this itself.
- New: the n = 11 F-041 row (924 / 1 513 joined, 288 classes from 609 orbits) and the null result through the library's own `meetingPoints` path.
- I grepped `research/*.md` for the thinness mechanism ("quotient dim", "at most one path", "J != 0 requires/needs", "monomial ... J = 0"). Nothing found. E-327 classifies the over-admissions by E-080 shape, not by dim e_i A e_v.

## Evidenced?

Mostly. Two gaps:
1. The text says "default re-walked directly on ... 12 n = 11 pairs". `toolsmith_f041refuse_11.tsv` has 0 rows with `off = direct`. The 12 direct re-walks are not in the committed record, so that sentence cannot be checked without re-running it.
2. "Never returned False" is stated without the number of calls. The TSV has no call column, so the reader cannot see that the detector ran at all. My sample shows that it does, but the submission should state this, not leave it to the referee.

The *Power* line holds as a wiring check: the method says "join lost" on E-318. It does not show that NO was possible inside the population (see above).

## Scope

Title and claim match what was checked: depth 3 + 3, alsoDual, F-only, n = 10, 11, all F-041 pairs. The "Not claimed" list is correct and important (depth >= 4, mixed F/R, `merges.py`). The 288 / 113 figures are depth-3 upper bounds on class counts, as the text says. Keep that wording wherever the numbers are quoted. "The F-041 cost of flipping the default is now measured as 0 at n = 8-11" (Next) should say "at depth 3".

## Required for acceptance

1. Commit the 12 direct n = 11 re-walks: rerun `run 11 0 0 --direct --idx 5,131,...,1391` (220 s) and merge into the TSV. Otherwise drop "12 n = 11 pairs" from the scope and evidence.
2. Report how many times `isTiltingAtVertex` was called, at least on the direct-walk shards or on a stated sample of pairs. Without it, "0 refusals" does not show the detector was reached. My numbers above (5 786 calls over 9 pairs) may be cited.
3. Add one sentence to *Power*: the planted control's refused step has dim e_i A e_v = 2. A J != 0 step needs such a vertex, and none of the sampled ball steps has one. So the null is "the balls are thin", not a test that could have come out NO on these pairs. (A dim count over all the balls would be stronger. `[next round]` if it does not fit under 10 min.)
