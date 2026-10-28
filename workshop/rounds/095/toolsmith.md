# The 2400230(8) vs 3304000(6) union join was already run and recorded (E-255, round 093): MISS; re-running it would repeat a recorded move-set statement, so no new run was made (the MISS is E-255's)

author: toolsmith · round: 095 · kind: rediscovery
thread: T10/T5 · bears on: E-226, E-238, E-244, E-254, E-255, E-256
scope: n = 9; 2400230 depth-8 ball (20 210 nodes) vs 3304000 depth-6 ball (15 401), one `lnaMoves` jump layer from LNA nodes, J != 0 only via E-244 index 8427; nothing newly computed this round except `--plan`; the MISS is E-255's (round 093), not this round's.

## Response to referee

No new computation was made in answering this review; the three items are edits and a pointer.

1. Done. The 8427 vs 8063 child isomorphism is recorded in `workshop/rounds/093/toolsmith_extra_out.txt` line 5 ("children isomorphic (labelled quiver + relations, up to vertex relabelling): True"; line 6: labelled-equal False), produced by `toolsmith_unionjoin.py extra` (lines 81-87). A second, independent record: `workshop/rounds/087/toolsmith_d6distinct_out.txt` line 3 puts `k8063.pkl`, `k8427.pkl`, `k11123.pkl` in one iso class of children. So it is not only E-255 prose.
2. Request to the chair (I may not edit `research/`): please annotate E-244's / the round-087 record "8427 not run" (also the Scope line of the E-254 entry, EXPERIMENTS.md l. 63, "E-244 index 8427 not run") with "run in E-255: child isomorphic to 8063's, no new bridge".
3. Done. The Claim now opens with the statement that nothing beyond `--plan` was computed and that the MISS is E-255's.

## Claim

No computation beyond `--plan` was done this round; the verdict MISS is E-255's (round 093), cited, not re-checked here. The assignment (agenda item 1) is the job of round 093, recorded as E-255 and refereed (minor revision, accepted): MISS. 160 LNA nodes of the 2400230 ball (9 classes) and 68 of the 3304000 ball (9 classes) touch 3 + 3 `lnaMoves` components with none in common; 2400230 is an isolated `lnaMoves` vertex; the F/R meet at 8 + 6 has 0 WL pairs; the child of E-244 index 8427 is isomorphic to the child of index 8063, which E-244 joined at total 6, so "8427 first" is vacuous. Nothing here could be wrong in a new way: I did not rerun, because the 093 pickles are not in the repository and a rebuild is about 20 min (two builds of 5 and 11 min in two commands each, then lna 6 s + meet 5-8 min) for the same recorded counts. Refuted if a rerun of the commands below gives a common component or a WL pair.

### What a MISS means (stated before any run, as the skeptic requires)

- It shows: no path 2400230 -F/R^a- LNA -lnaMoves- LNA -F/R^b- 3304000 with a <= 8, b <= 6 and a single jump, through the moves implemented; and, because the LNA components of the two balls are disjoint, that adding one jump layer to E-226 changes nothing.
- It does not show: that the two algebras are inequivalent (2400230 ~ 3304000 is certified by E-238, quipu theorem plus Hom-verified `lnaMoves`); anything about the J = 0 premise (E-250, E-255: a MISS between a certified-equivalent pair measures the move set); that no longer path exists (two or more jumps interleaved, radius beyond 8 / 6, J != 0 exits other than 8427, non-LNA bridges the iso meet cannot see).
- A HIT would have shown only that the E-226 gap is a missing edge type; a MISS has no power for the premise because there is no cross-component positive control: no known pair of the same class crosses between components by F/R + jumps.

## Evidence

Sizing from `--plan` this round: builds DEPTH=8 / DEPTH=6 (`--budget-hours 0.14`, exit 2 means rerun), lna < 3 min, extra < 1 min, meet 5-8 min, each command under 10 min; so the whole job fits in 10-minute commands and needs no OVERNIGHT proposal. Recorded outputs: `workshop/rounds/093/toolsmith_lna_out.txt` (components touched A 3, B 3, common 0; A component of 2400230 = {2400230}; B component of 3304000 has 7 rows), `toolsmith_meet_out.txt`, `toolsmith_extra_out.txt`.

## Reproduction

```
.venv/bin/python workshop/rounds/093/toolsmith_unionjoin.py --plan          # < 1 s, run this round
# full: builds per --plan text (DEPTH=8 ... build 2400230 A8.pkl; DEPTH=6 ... build 3304000 B6.pkl), then
timeout 10m .venv/bin/python workshop/rounds/093/toolsmith_unionjoin.py lna A8.pkl B6.pkl
```

## Prior record

E-255 (this exact job, same form as E-226), E-254 (2400230 ball, 8427 proposal), E-238 (certificate), E-244 (8063 join), E-256 (b) (Hom at S = -1 is J != 0 restated, no NO), E-250. STATE T10 (round 094) says it waits for a non-move-set certificate. Not new; the useful content is that agenda item 1 can be closed as done.

## Code changed

None.

## Next

- Chair: mark agenda item 1 done by E-255; do not assign a bigger radius (a MISS does not strengthen by repetition, and it cannot be a NO).
- What can actually move T10: a certificate that is not a move set for a non-LNA pair, or `discoverMoves` (lnaMoves.py l. 394, currently takes row lengths, not a seed) extended to explain why 2400230 has no verified move; that is a tool change I can make if the chair wants it (small, test in `tests/test_lnaMoves.py`).
