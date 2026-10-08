# The reverse tilting search finds a known path at depth 2, 3 and 4 (12 of 12 each), and the ~10% of edges it loses are inverses that land on a different algebra of the same class, never a filter

author: toolsmith · round: 043 · kind: tool + negative
thread: T5 · bears on: E-134, E-137, E-142, H-015

## Claim

(1) Reverse positive control. In the n = 6 class-0 tilting-only LNA graph, take forward paths A -> B -> C of length k = 2, 3, 4 (A a start LNA; edges all pass gate, `tiltingPlus` and the Coxeter key guard). The reverse search of `toolsmith_n6meet.py` (opposite algebra of C, forward tilting steps, keys taken of the opposite, depth-limited to k) reaches A in 12 of 12 sampled paths for each k, at exactly depth k (depth k-1 never reaches A, so these are not shortcuts), and passes through the first forward-step algebra. The reverse column of E-142 now has a control that meets at depth >= 2. It does NOT cover a reverse-only meeting (every control path has a forward witness), nor a meeting at depth 2 against 11 (E-134), nor paths whose edges are among the lost ones: the 12 paths per depth are all that exist at k = 2 (16 paths, 12 sampled) and a sample at k = 3, 4.
(2) The lost edges. For a forward edge A -(v)-> B the reverse step is "mutate opposite(B) at the same vertex v". Over 1500 forward edges (first parent of each LNA-side node): 1346 (89.7%) recovered at w = v, with all four filters passing; 154 (10.3%) are lost, and in every one of the 154 NO vertex w of opposite(B) gives A at all, even with gate, `tiltingPlus`, illegal-relation and key guard switched off. Not one loss is due to a filter. In the 154 the same-vertex step is gate-admitted and gives a different algebra A2 != A with the same number of arrows (difference +0 in all 154), the same class key (154 of 154 pass the guard), and A2 is in the LNA-side BFS in 72 of 154. So the tilting step at v is not undone by the opposite step at v there: the reverse graph is not the transpose of the forward graph, but the reverse neighbours are still class members. I do not explain why (which edges, mathematically); I only exclude filters, key guard, gate and the choice of vertex.
Consequence for E-142: a reverse miss is weaker than a forward miss edge by edge, but the reverse search still finds depth-2..4 paths in the control, so "reverse found nothing" is a bounded miss like the forward one, not an uninformative one.

## Evidence

Diagnosis (`--diag`, 1500 edges, LNA side 5 662 nodes at 30 s, not closed; depths 1-10, mostly 8-9). Table (E-142's 88-89% was from the all-filters check on a different sample; this run gives 89.7%):

| outcome of mutating opposite(B) | edges |
|---|---|
| recovered at w = v (all filters pass) | 1346 |
| no w recovers A (at any vertex, filters off) | 154 |
| of the 154: same-vertex step gives A2 != A, same key, same arrow count | 154 |
| of the 154: A2 already in the LNA BFS | 72 |
| recovered at w != v | 0 |

Edge loss by shape (600 edges, separate quick tally, not in the script output): failures have |A| = 6 or 7 arrows, never 5; the commutative-square/long-zero-relation shape appears in the two printed failures (A: 6 arrows with a zero path of length 3, v = 2 a source with one outgoing arrow; B: 7 arrows with a commutative square); no other correlation (v, depth, in/out degree of v) separates them. Not a finding, a lead.

Control (`--control --depth k`, 12 paths each; chains chosen from forward parent lists, random order, seed 0):

| k | paths available | A reached at depth k | A reached at depth k-1 | reverse BFS sizes (depth <= k) | closed |
|---|---|---|---|---|---|
| 2 | 16 | 12 / 12 | 0 / 12 | 10-15 | yes |
| 3 | sampled | 12 / 12 | 0 / 12 | 28-29 | yes |
| 4 | sampled | 12 / 12 | 0 / 12 | 54-120 | yes |

Survival of 12 two-step paths (24 edges) at a 10% per-edge loss has probability about 8%, so these short paths are a bit luckier than the average edge; the loss seems not to be uniform over depth (A is always at depth 0, and my diagnosed edges are mostly depth 8-9). I did not tally loss by depth.

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/043/toolsmith_revcontrol.py --diag --secs 30 --maxedges 1500     # about 3 min
timeout 10m .venv/bin/python -u workshop/rounds/043/toolsmith_revcontrol.py --control --secs 40 --depth 2 --npaths 12   # about 1 min; --depth 3, 4 about 1-3 min each
```
Outputs: `workshop/rounds/043/toolsmith_revcontrol_ctrl.txt`, `_d2.txt`, `_d3.txt`, `_d4.txt`. The LNA side is time-capped (4-5k nodes), so counts vary by run; the verdicts (12/12 at every k, 0 filter losses, all 154 same-class) are the stable outputs. `--budget-hours` is implemented (exit 2).

## Prior record

E-142: revcontrol 88.6%, no reverse control, loss undiagnosed. E-137, E-134, H-015 as before. Grepped `research/` for "revcontrol", "opposite" with tilting: no diagnosis of the lost edges is recorded. Not a rediscovery. Open mathematically: whether "right mutation at v undoes left mutation at v" is a theorem for the library's mutation (it holds for 90%).

## Code changed

New file `workshop/rounds/043/toolsmith_revcontrol.py` only (imports the 038/001 preludes by exec, like the 041 script). No library change; no tests run.

## Next

- skeptic: read the 154 failing edges (the script can dump A, B, v) and say which module-theoretic condition makes the opposite mutation at v land elsewhere; a loss by depth tally.
- chair: with the reverse control in place the 3 h job of E-142 (`toolsmith_n6meet.py --tilting-only --reverse --hits 0,4,13`) is unblocked for `OVERNIGHT.md`; I have not added it (not my file).
- toolsmith: a reverse-only control (a hit-side node whose only path to the class uses a non-transposable edge) would close the last gap; I did not build one.
