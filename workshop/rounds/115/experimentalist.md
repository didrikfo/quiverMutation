# K0 = min(h, n-4-h) holds at n = 19 h = 5, 6 and n = 15 h = 3, the falsifier runs pre-registered in H-026 (key level, head ends, K >= K1 only); key classes of a two-3 row are not of the lone-3 shape (last split K = 2, 1, none and >= 3 seen, no law)

## Response to referee

Verdict minor revision; four required items.
1. libdel at n = 19: done. `libdel` run on n = 19 h = 6 (`timeout 10m ... kloc.py run 19 6 h=6 libdel`, 47 s, output `kloc_n19h6_libdel_out.txt`): 1982 class rows, 0 mismatches against `pwh.removeVertex`; table identical to before. n = 19 h = 5 (157 s plus libdel) not re-run with libdel; the deletion code is now checked at n = 19 h = 6, n = 15 h = 3 and n = 13 (2,7) only.
2. Title: fixed. Now says the three points are H-026's pre-registered falsifiers, adds "none" to the non-lone thresholds, and says "key classes of a two-3 row". Scope line says "class" means key class, not a move orbit (unchecked), and that (2,7), (3,9) thresholds are limited by missing high-K rows.
3. H-026 cited in the Claim: the result is a confirmation of an existing prediction (7 to 10 fits); the non-lone tables are the only new, exploratory part.
4. The "max keys = 3, depth tied to min(p, ...)" reading from two points is removed.
Also accepted: "last split on K >= K1" wording (kept in scope and Claim); K < K0 not run; P4 was refined after P1-P3 and held trivially.

author: experimentalist · round: 115 · kind: result
thread: S-1 · bears on: H-026 (E-301 off-centre law), F-candidate non-lone check
scope: key level ("class" = Coxeter key class only; not checked to be one move orbit), head-end deletion of vertex 1, lone-3 classes n = 19 h = 5 (K >= 5), n = 19 h = 6 (K >= 6), n = 15 h = 3 (K >= 2); non-lone: five key classes of two-3 rows (n = 13 x4, K >= 1; n = 15 x1, K >= 2); the (2,7) and (3,9) thresholds are limited by absent high-K rows, not measured. K < K1 not run (n = 19 K1 = 3 would take ~50 min); so K0 is "last split on K >= K1"; with K1 = K0 in both n = 19 runs they show a split at K0 and none above, nothing about K < K0. Not class level, not tail.
power: the tables can say NO and do: the non-lone classes give 1 key (n = 13 p=3,6) or 3 keys (p=4,8), so the counter does not always report 2/1. Control: own deletion vs pwh.removeVertex(n,r,1), 0 mismatches (n = 15 h=3: 4597 rows; n = 13 p=2,7: 504 rows; n = 19 h=6 now also checked: 1982 rows, 0 mismatches; not run with libdel at n = 19 h=5, (4,8), (3,6), (3,9)). No planted-error test of the threshold itself.
backlog: n/a

## Claim

These three points are exactly the falsifier runs pre-registered in H-026 (research/HYPOTHESES.md: "n = 19 h = 5 or 6 (predicted 5, 6), or n = 15 h = 3 (predicted 3)"), so this is confirmation of an existing prediction, not a new claim. The three predictions held: K0 = 5 (n = 19, h = 5), 6 (n = 19, h = 6), 3 (n = 15, h = 3), each equal to min(h, n-4-h), each with 2 image keys at K = K0 and a single key for every K > K0. With the seven points of round 114 this makes ten fits of the off-centre law; still a fit, not a theorem. Refuted by a lone-3 class where the last split K differs. The non-lone part is exploratory (the one concrete prediction, below, was a weak test): it shows that classes with two separated 3s behave differently but gives no law.

## Evidence

Lone-3 classes, image-key counts per K (key labels arbitrary), head ends, K = number of leading vertices no relation starts at:

| n | h (class) | rows in class | K | counts per image key | keys for K>K0 |
|---|---|---|---|---|---|
| 19 | 5 (5,10) | 5344 | 5 | 2416 / 1025 | 1 key at K = 6..10: 1091, 486, 222, 71, 33 |
| 19 | 6 (6,9) | 1982 | 6 | 856 / 513 | 1 key at K = 7, 8, 9: 412, 136, 65 |
| 15 | 3 (3,8) | 4597 | 2 | 2183 / 522 | |
| 15 | 3 | | 3 | 956 / 257 | 1 key at K = 4..8: 410, 171, 68, 21, 9 |

Observations: n = 15 h = 3 is split at K = 2 and 3 and not at 4 (K0 = 3 = h). The counts 412, 136, 65 (n = 19 h = 6) equal the n = 17 central table of round 114 and 513 equals a count in the n = 21 table; 257, 522, 1025 recur from earlier tables. Unexplained, noted only.

Non-lone key classes of a two-3 row. Definition: row r0 of n-2 entries, zero except 3 at two 0-based positions p < q; class = all rows (relation-length rows of the LNA with n vertices, K >= K1) with the same Coxeter key as r0. Note this is a key class, which may contain rows not of the two-3 shape. Prediction (logged first, only for p=2, q=7, n = 13): split up to K0' = 2 and counts different from the lone-3 sequence. Result:

| n | (p,q) | rows | K=1 | K=2 | K=3 | K=4 | max split K |
|---|---|---|---|---|---|---|---|
| 13 | (2,7) | 504 | 235/96 | 158/15 | none | none | 2 (no rows K >= 3, so "1 key above" untestable) |
| 13 | (4,8) | 597 | 191/177/60 | 94 | 44 | 31 | 1 (3 keys at K=1) |
| 13 | (3,6) | 147 | 78 | 39 | 30 | | none (1 key everywhere) |
| 15 | (3,9) | 2064 | K=2: 648/282/679 | K=3: 413/42 | none at K>=4 | | >= 3 (3 keys at K=2, 2 at K=3; K=4.. no class rows) |

So the P4 prediction "K0' = 2" held trivially for (2,7) only because the class has no rows with K >= 3; the other three cases were not predicted. What they show: the number of image keys can be 3, the split can end at K = 1 or earlier than the first 3's position (4,8) , and there can be no split ((3,6)). The lone-3 counts 780/265/129/65 do not reappear. No formula guess is supported by four classes.

What would count as a falsifier of the law: any lone-3 class (h, n-4-h) with last split K != min(h, n-4-h). Of the non-lone: a class with two 3s at positions (p,q) whose last split exceeds p + something is not a falsifier of anything yet, since no non-lone law is stated.

## Reproduction (repo root; all `TAIL=0`)

```
export TAIL=0
timeout 10m .venv/bin/python workshop/rounds/115/experimentalist_kloc.py run 19 5 h=5     # 157 s
timeout 10m .venv/bin/python workshop/rounds/115/experimentalist_kloc.py run 19 6 h=6     # 51 s
timeout 10m .venv/bin/python workshop/rounds/115/experimentalist_kloc.py run 15 2 h=3 libdel   # 42 s
timeout 10m .venv/bin/python workshop/rounds/115/experimentalist_nonlone.py run 13 1 p=2 p=7 libdel   # 9 s
timeout 10m .venv/bin/python workshop/rounds/115/experimentalist_nonlone.py run 15 2 p=3 p=9   # 30 s
```
Outputs: `kloc_n19h5_out.txt`, `kloc_n19h6_out.txt`, `kloc_n15h3_out.txt`, `nonlone_extra_out.txt` in the same directory (n = 13 (2,7) output is in the log above, rerun to get it).

## Prior record

E-299, E-301 (law fitted on 7 points), round 114 tail submission, R-015 (not an issue here: the script holds tuples, not digit strings). The n = 19 h = 5, 6 and n = 15 h = 3 points are new; I did not grep research/ for a non-lone two-3 class table and found none cited in STATE.md.

## Code changed

New `experimentalist_kloc.py` (copy of round 114's tail script, unchanged logic) and `experimentalist_nonlone.py` (adds `p=` positions to define the target row). No library changes, no tests.

## Next

Class level for the lone-3 law; a non-lone class big enough to have K > K0' rows (p,q small, n >= 15) with predictions logged for the key count; theorist: derive min(h, n-4-h).
