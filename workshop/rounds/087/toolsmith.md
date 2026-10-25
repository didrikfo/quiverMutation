# The 3304000 ball to depth 6 holds 6 key-kept perI J != 0 steps among 11 355 parents (3 isomorphism classes of child), all real; five of the six children (all three classes) join 3304000 at total 6-8 by J = 0 moves

author: toolsmith · round: 087 · kind: result
thread: T10 · bears on: E-239, E-234, E-229, E-204, E-219/E-226 (J = 0 premise)
scope: (round 087 response: joins run for 4606, 8063, 8170, 8196, 11123; 8427 not run) n = 9, root 3304000, E-212 `moves` ball to depth 6 (15 401 nodes), all 11 355 depth-6 parents scanned for perI-flagged gate-admissible F/R steps (key-kept = child key equals the class key; key of 3060000 used, equal by E-239); Hom at S = -1, 0; join = J = 0 chain to depth-5 root ball with iso meet, child ball depth 3, one representative per child class; generation by T assumed. Not covered: depth <= 5 parents (E-239), 3060000 at depth 6, multi-vertex J != 0 not flagged by perI, a second J != 0 step after the first, any NO test.

## Response to referee

1. Missing RANGE lines: done. 0-4000 and 4000-7828 re-run from par6.pkl and appended to `toolsmith_d6scan_out.txt`: 0-4000 186 moved, 0 kept (187 s here, two scans in parallel); 4000-7828 921 moved, 1 kept (4606), DONE at 277 s with explicit end 7828 (first run was a budget stop at 7828; this one was given 7828 as its end). Totals 186+921+61+656 = 1 824 moved, 6 kept reproduce. The ball sizes line (1, 11, 64, 282, 1088, 4046, 15 401) remains stated from the build, not in a file.
2. Join move lists: done. The quoted lists were from an earlier meet; several meets of equal total exist. Replaced in the output file by what `skeptic_join.py meet` prints now (same as the referee's): 4606 total 8 `F2 F1 F3 F6 F5 || R1 R2 R4`; 8063 total 6 `F2 F1 F6 || R2 R4 R5`; 8170 total 7 `F2 F1 F6 F5 || R2 R2 R4`. Totals unchanged.
3. Title narrowed: it now says five of six children (the three class representatives and two more) join.
4. Optional join of other class members: done for 8196 (class of 8170; total 7, `F2 F1 F6 F5 || R4 R2 R4`) and 11123 (class of 8063; total 6, `F2 F1 F6 || R4 R7 R5`), both Hom-verified, iso meet True. The isomorphic-children step is thus checked twice. 8427 not run.
Not done: independent reproduction of the build and the J != 0 definition recheck (referee's note; unchanged).

## Claim

Sizing: the depth-6 scan fits in commands of at most 10 minutes once the build is checkpointed; no OVERNIGHT proposal is needed for this side. Build 8 min (budget stop) + 2.7 min (resume), scan 4 shards of 1-2.5 min. Result: among the 11 355 depth-6 parents of the 3304000 ball, 1 830 perI candidate steps (1 824 key-moved + 6 key-kept); exactly SIX key-kept (parent indices 4606, 8063, 8170, 8196, 8427, 11123): 6 of 11 355 parents (0.053%), versus 1 of 2 958 at depth 5 (E-239, 0.034%) and 10 of 4 770 for 3060000 at depth 5 (0.21%); these are flagged-and-kept counts and no rate comparison between the roots is claimed (events are not independent, see classes). All six are real J != 0 steps: tiltingPlus False, Hom(T,T[-1]) != 0 at S = -1 (total 2, 2, 4, 4, 2, 2), Coxeter key parent = child, Cartan det 1/1. Child isomorphism classes (quiver+relations up to vertex relabelling): 3 among the six, {4606}, {8063, 8427, 11123}, {8170, 8196}; the E-233 child (3060000) is a fourth, distinct class; classes {8170, 8196} have 9 arrows and 2 relations (a shape not seen at depth 5: two perI vertices carry multiplicity 2 at vertex 2 resp. 4), the other four have 10 arrows and 4 relations. Children of 4606, 8063, 8170, 8196, 11123 join 3304000 by J = 0, Hom-verified chains with iso meet (root side || child side as printed by `skeptic_join.py meet`; equal-total meets are not unique): 4606 total 8 (F2 F1 F3 F6 F5 || R1 R2 R4), 8063 total 6 (F2 F1 F6 || R2 R4 R5), 8170 total 7 (F2 F1 F6 F5 || R2 R2 R4), 8196 total 7, 11123 total 6. So again a key-kept J != 0 step exists and does not leave the J = 0 class of its own root; the join checks only that (Hom + Cartan have no NO power, E-228). Does not claim: the other members of each class were joined (isomorphic children have isomorphic balls, but 8427 was not run), that no J != 0 step follows another (child of a J != 0 step as parent not scanned), anything about index 2026's class (E-239's pickle was not rebuilt; its child may coincide with one of the three). Would be refuted by: a child failing to join its own root's J = 0 component given enough depth.

## Evidence

| shard (parents) | perI moved | kept | time |
|---|---|---|---|
| 0-4000 | 186 | 0 | 101 s |
| 4000-7828 | 921 | 1 (4606) | 8.4 min budget stop |
| 7828-8000 | 61 | 0 | 13 s |
| 8000-11355 | 656 | 5 | 140 s |
| total 11 355 | 1 824 | 6 | |

Ball sizes 1, 11, 64, 282, 1088, 4046, 15 401 (depth 5 matches E-239 exactly; so does the 10 564 gate / 69 J count in the root ball depth 5 re-built for the join). Build timings: depth 5 160 s (slower than E-239's 116 s on this machine), depth 6 expansion of 2 958 nodes 155 s after the resume. The kept parents have paths of length 6 (list in `toolsmith_d6scan_out.txt`); steps are R at v = 7 (four), R at 9, R at 2; parents are in the diamond family (perI on {5, 8} or {5, 9}), 8170/8196 add a second multiplicity-2 perI vertex. Positive control: the scan code is E-239's `jcands` unchanged; the depth-5 index 2026 was not re-found (depth 6 parents only); the control of E-239 (`ctl`) was not re-run. Hom outputs by `082/skeptic_kept.py test`; joins by `082/skeptic_join.py`.

## Reproduction

From the repository root; every command at most 10 min (`--budget-hours 0.14`, exit 2 = rerun the same command):
```
mkdir -p /tmp/ts87
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/087/toolsmith_d6scan.py --plan
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/087/toolsmith_d6scan.py build 3304000 /tmp/ts87/st.pkl --budget-hours 0.14   # exit 2 after 8.2 min at depth 5 idx 2007/2958; rerun: 2.7 min
.venv/bin/python -W ignore workshop/rounds/087/toolsmith_d6scan.py parents /tmp/ts87/st.pkl /tmp/ts87/par6.pkl                     # 11 355
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/087/toolsmith_d6scan.py scan /tmp/ts87/par6.pkl 0 4000 --budget-hours 0.14   # 101 s; then 4000 8000 (exit 2 at 7828: resume 7828 8000), 8000 11355 (140 s)
.venv/bin/python -W ignore workshop/rounds/085/toolsmith_d5pkl.py /tmp/ts87/par6.pkl /tmp/ts87 4606 8063 8170 8196 8427 11123
.venv/bin/python -W ignore workshop/rounds/085/toolsmith_distinct.py /tmp/ts87/k{4606,8063,8170,8196,8427,11123}.pkl workshop/rounds/082/skeptic_kept.pkl   # 2 s
.venv/bin/python -W ignore workshop/rounds/082/skeptic_kept.py test /tmp/ts87/k4606.pkl   # Hom, per index
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side R 3304000 5 560 /tmp/ts87/ckR     # 152 s
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side C pkl:/tmp/ts87/k4606.pkl:child 3 150 /tmp/ts87/ckC4606   # 15 s; then meet ... (8063, 8170, 8196, 11123 likewise)
```
Outputs kept: `toolsmith_d6scan_out.txt` (scan lines, truncated to 200 chars per line), `toolsmith_d6distinct_out.txt`.

## Prior record

E-239: 1 kept at depth 5 (this extends to depth 6; the depth-6 scan covers new parents only, depth <= 5 not rescanned). E-229: 0 of 69 at depth <= 4. E-204's "no J != 0 step enters the balls" is false for both roots (E-234, E-239); six more events here. Grep of `research/RETRACTIONS.md` not needed for these terms; nothing found touching depth 6.

## Code changed

New script only: `workshop/rounds/087/toolsmith_d6scan.py` (E-239's scanner with `--plan`, resumable build checkpointing every 300 expansions, `--budget-hours` with exit 2, scan resume index). No library change, no tests run. 3060000 depth 6 not run: same script, `KEYROOT` default is already the shared key; sizing: build about 12 min of expansion (4 770 front nodes at about 0.1 s) in two commands, parents about 17 000-19 000, scan about 7 min in 3 commands.

## Next

Skeptic: recheck index 4606 and the 8170 class (the new shape). Toolsmith: the J != 0 child as a parent (scan one level below the six children) to see whether a second J != 0 step follows; 3060000 depth 6. Experimentalist: a NO-capable certificate is still the gap (these joins are J = 0).
