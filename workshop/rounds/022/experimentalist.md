# The long-sided square is selective for rejection: 0 of 479 761 tilting steps on guarded walks at n = 5..7 have it, against every rejecting parent

author: experimentalist · round: 022 · kind: result
thread: T5 · bears on: E-097 (its named missing control), E-084, E-095

## Claim

On guarded walks (class 0 and class 1 where run) at n = 5, 6, 7, no step at which `tiltingPlus` holds has E-097's long-sided square at the mutated
vertex v (`hasLongSquare`: v has one out-arrow v->e and a relation of >= 2 paths from one start, all ending x, v, e with distinct x), while every
rejecting step has it. Strict A5 (arrows a>b, a>c, b>v, c>v, v>e) is NOT selective: it occurs at 0.5 to 1.2 % of tilting steps as well. So the
long square, not the strict A5 quiver shape, is what separates the rejecting parents here. It does not claim the long square implies rejection (the
control counts only steps the walks reach, and the test needs a relation, so quiver-only A5 steps never pass it), nor that the 0 holds off the
walks (E-078's hand-built family contains tilting steps with a long relation; not run here).

## Evidence

Guarded BFS (as E-097's script), every gate-admitted step, de-duplicated on (canonical key of parent, v); 480 s cap per run, so counts are lower bounds.

| n, class | expansions | tilting steps (distinct parent,v) | tilting with long square | tilting with strict A5 | rejecting (parent,v) | rejecting with long square | rejecting with strict A5 |
|---|---|---|---|---|---|---|---|
| 5, 0 | 6 240 (done) | 16 620 | 0 | 60 | 0 | - | - |
| 5, 1 | 5 460 (done) | 13 680 | 0 | 0 | 0 | - | - |
| 6, 0 | 41 825 (cap) | 150 592 | 0 | 807 (0.54 %) | 1 842 | 1 842 (100 %) | 1 299 (70.5 %) |
| 6, 1 | 58 749 (cap) | 199 988 | 0 | 2 306 (1.15 %) | 0 | - | - |
| 7, 0 | 24 704 (cap) | 98 881 | 0 | 1 316 (1.33 %) | 262 | 262 (100 %) | 262 (100 %) |

- Every one of the 807 + 2 306 + 1 316 + 60 tilting steps with strict A5 lacks the long square ("A5 but not long"), so in the tilting steps A5 occurs without any relation of that shape.
- Every rejecting step has v with exactly one out-arrow (1 842 of 1 842; 262 of 262); among tilting steps 99 631 of 150 592 (n = 6 c0) do, and none has the long square, so one-out-arrow is not the separating feature.
- Strict A5 at rejecting parents: 70.5 % at n = 6 (E-097: 767 of 1 123 = 68.3 %, different cap and run), 100 % at n = 7 (E-097: 156 of 156). Rate with A5 only: 1 299 / (1 299 + 807) = 62 % of A5-quiver steps at n = 6 c0 reject; with the long square, 100 % (1 842 of 1 842 + 0).
- The n = 6 rejecting count here (1 842 distinct (parent,v), one per parent as E-095 found) is larger than E-097's 1 123 only because this run did less work per step and so reached more of the walk in 480 s (41 825 expansions vs E-097's 46 584 algebras with a heavier script); the rates are what to read.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/022/experimentalist_shapectl.py 5 --class 0          # 24 s
timeout 10m .venv/bin/python workshop/rounds/022/experimentalist_shapectl.py 6 --class 0 --budget-sec 480   # 480 s, 4 runs in parallel
timeout 10m .venv/bin/python workshop/rounds/022/experimentalist_shapectl.py 7 --class 0 --budget-sec 480
```
Outputs: `workshop/rounds/022/experimentalist_shapectl_n{5c1,6c0,6c1,7c0}.txt` (n5c0 in the submission text above, 24 s run, not saved to file). The cap makes counts non-deterministic; use `--maxexp N` for a fixed cap.

## Prior record

E-097's Limits name this control as missing ("no control of how many tilting parents have the long square"; `hasLongSquare` loose). E-084/E-095 call rejecting parents A5-shaped; E-097 corrects this to long square. Grep of `research/` finds no tilting-side count. New: the control, 0 of 479 761.

## Code changed

None in the library. New script `workshop/rounds/022/experimentalist_shapectl.py` (shape tests copied from `rounds/021/scholar_step7_entries.py`; no kernel computation, so about 2x faster per step). No tests touched.

## Next

- Theorist: 0 of 479 761 suggests "tiltingPlus fails at v on the walk iff v has a long square" is a conjecture worth a proof or counterexample: the long square (relation a->..->x->v->e from two paths, v one out-arrow) is E-066's commutativity element into a vertex with one outgoing arrow. Skeptic: the test is loose; find a long-square tilting step off the walks (E-078 `long` family, rounds with relations of 3 paths).
- Experimentalist: converse on the walks is already this table (all rejecting have it); run n = 6 classes 2-3, n = 7 classes 1+ with `--maxexp`; n = 8 c2's 2 rejections (path (17,8,5,6,8,8,2,5)) shape check.
