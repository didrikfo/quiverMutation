# The E-180 dimension assertion holds on all 370 printed n = 7 path edges and all 35 n = 10 E-169 edges (0 failures, 0 unevaluable)

author: toolsmith · round: 063 · kind: result
thread: T10 · bears on: E-180, E-186, E-174, H-015
scope: n = 7 classes 1 and 2 (the 26 + 18 E-157/E-160/E-163 plan paths of `rounds/057/toolsmith_endt2_run.py`: 214 + 156 = 370 edges, both sides of each path); n = 10 the 5 E-169 paths (35 edges). Dimension assertion only; the equal-dimension rejection power on these edges is NOT tested (E-186 does that only at the failing steps).

## Claim

For every one of these 405 edges (parent algebra a, vertex v, reduced child c = the repo's mutation) the sum over vertex pairs of dim K Q_c/I_c, computed by linear algebra from the child's arrows and relation list (`quotient_dims`, `dimcheck` of E-180), equals the child's dimension matrix pair by pair, and that matrix equals dim End(T) pair by pair (`endquiver`). So the E-180 assertion that `compare2` makes inside its verdict is satisfied everywhere the E-186 note said it was uncovered. It does not claim the End(T) iso test passes beyond what E-174 recorded, that relation sets are correct (the assertion cannot see a wrong ideal of the right dimension, E-186), or anything about J.
A refutation would be one edge with `dimcheck` not None or with End(T) dimensions different from the child's.

## Evidence

| set | edges | dimcheck holds | End(T) dims = child dims | unevaluable |
|---|---|---|---|---|
| n = 7 class 1 (26 paths) | 214 | 214 | 214 | 0 |
| n = 7 class 2 (18 paths) | 156 | 156 | 156 | 0 |
| n = 10 E-169 (5 paths) | 35 | 35 | 35 | 0 |

Per-edge lines (dim Q/I, dim End(T), dim child, number of relations) are in `toolsmith_dimedges_c1.txt`, `_c2.txt`, `_n10.txt` (all under 200 KB). Every edge has at least one relation in the printed children where relations exist (no edge with `nrels 0` at n = 7), so the assertion was not vacuous there.
Power control on the same edges (`toolsmith_dimedges_ctrl.py`): deleting each single relation of the child on every edge gives 833 (class 1) + 395 (class 2) = 1 228 cases, rejected 1 228, none redundant (all raise dim K Q/I' by 1..5). That is the E-180 control repeated on path edges (previously only failing steps); it says nothing about equal-dimension wrong ideals.
The script reuses `compare2`'s own pieces (`TiltEnd`, `endquiver`, `child_info`, `quotient_dims`, `dimcheck`) unchanged; the n = 10 rerun of `rounds/057/toolsmith_endt2_n10.py` (which calls the updated `compare2`) also gives 35 iso verdicts, i.e. the assertion and the iso test agree there.
Pickles were rebuilt (`/tmp/tsm/c1.pkl` ~9 min, `c2.pkl` ~7 min, resumable collect runs of round 050; the counts of paths agree with E-174: 370 edges).

## Reproduction

```
DEADLINE=520 timeout 10m .venv/bin/python -u workshop/rounds/050/toolsmith_collect.py 7 CLS 20000 /tmp/tsm/cCLS.pkl 100   # CLS = 1, 2; rerun until it finishes
CLS=1 timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges.py n7 0:26    # ~1 min
CLS=2 timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges.py n7 0:18    # ~1 min
timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges.py n10             # seconds
CLS=1 timeout 10m .venv/bin/python -u workshop/rounds/063/toolsmith_dimedges_ctrl.py      # drop-one control, same for CLS=2
```

## Prior record

E-180 (assertion on the 25 failing steps, 125/125 drops rejected), E-174 (iso on 370/370 and 35/35 edges, before the assertion existed), E-186 (listed these edges as not covered). New: the assertion itself on the edges, and the single-drop control at edge level. Consistent with E-181 (dim End(T) = Cartan(rewrite) for n <= 8). No retraction touched.

## Code changed

New: `workshop/rounds/063/toolsmith_dimedges.py`, `toolsmith_dimedges_ctrl.py` (+ outputs). No library change, no tests needed.

## Next

Experimentalist/skeptic: the equal-dimension wrong-ideal control (E-186 family) on a sample of path edges with parallel arrows, if anyone wants the edge-level power; the premise J = 0 still rests on Hom(T,T[-1]) and the unread Thm 1.1 (scholar). Toolsmith housekeeping unchanged: `canonicalKey` cap 5040 and docstring rewords.
