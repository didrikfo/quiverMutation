# Step 8, the last of E-201's 13 n = 7 class-1 steps, gets a Hom-verified 12-edge tilting chain (up to relation coefficients, generation assumed): 13 of 13 now joined, a YES-side consistency check

author: toolsmith · round: 114 · kind: result
thread: T10 · bears on: H-015, R-candidate "a key-guarded walk stays in one class", E-201, E-297, E-298
scope: n = 7, key class 1, step 8 only (v = 4, parent and child both with parallel arrows, dim Cartan 28 and 19) of the 13 D = 0 key-kept J != 0 steps of the `069/skeptic_c1_13.py` checkpoint (nexp 23000, rebuilt this round). Bidirectional search of vertex tilting edges (F/R at one vertex), per-side depth 6 (balls 16.9 k / 17.2 k algebras at the meet); verification uses no repo gate/J/key filter. Meet isomorphism ignores relation coefficients; generation of K^b(proj) by T assumed via loop-free mutated vertices (E-168) and K0 det. The 12 other steps are as in E-297.
power: the verifier says NO on planted errors: B-side single-edge corruptions 78, accepted 1 (that one, B edge 4 replaced by R1, is a valid alternative tilting edge into an isomorphic algebra, a true YES, as in E-297), rejected 77 (meet iso 77, "step impossible" 55, Hom(T,T[+-1]) != 0 in 1, K0 and Cartan and loop 0; meet iso alone rejected 21). End swap: B replaced by each of the 12 other steps' children, rejected 12 of 12 (meet iso in all; Hom also in 4). A-side corruptions (84) not run: the 10 min limit killed the job (A has a triple arrow, each verify is slow). No power for "the child leaves the class" (no key-equal inequivalent pair at n <= 7).

## Claim
From parent_8 apply `R4 F1 R7 F1 F4 F3` to reach M; from child_8 apply `F1 F3 R7 R4 R3 F5` to reach M' ~ M; this is a chain of 12 vertex tilting edges. Every edge has Hom_K(T,T[-1]) = Hom_K(T,T[1]) = 0 (chain maps mod homotopy, `050/skeptic_tilt.py`), a loop-free mutated vertex (no `!` flag), K0 determinant of the 7 summands +-1, and dim Hom(T_i,T_j) equal to the labelled Cartan matrix of the next algebra; the two halves meet in structure-graph isomorphic algebras with no parallel arrows. So parent_8 and child_8 are derived equivalent up to relation coefficients (meet iso ignores them) and assuming generation by T (E-168, Okuyama-Rickard). With E-297's 12 chains this gives a chain for all 13 steps. The direct step F4 still fails Hom (dim Hom(T,T[-1]) = 1; E-201, expected). This is a YES-side consistency check: no NO case exists at n <= 7, it could not have refuted the R-candidate and does not move it beyond E-228/E-297's level; the meet iso did the rejecting in the controls (Hom 1 of 77). Not claimed: any n >= 8 or other class; the other 9 chains' own controls; the A-side corruption control for this chain. A NO would have been an edge with Hom(T,T[+-1]) != 0, K0 det != +-1, Cartan mismatch or a non-isomorphic meet.

## Evidence
Why E-297's step 8 stalled and what changed: its serial search had balls of 5.5 k / 9.8 k after about 32 min (about 3 algebras/s, parallel-arrow algebras slow). `toolsmith_s8par.py` is the same BFS (same `moves` filter, same `key`) run on 4 cores with balanced sides (always expand the smaller ball) and resumable checkpoints. Ball sizes by level (A: 1, 7, 45, 249, 1227, ~5 k, 16.9 k at the meet in level 6; B: 1, 7, 38, 178, 816, 3.7 k, 17.2 k): growth about 4.5 per level, so depth 5 + 5 had no meet and the chain needs 6 + 6 in the BFS metric (found 6 + 6 = 12 edges). Obstruction for the old run: not a mathematical one, only throughput; the step needed depth 6 on both sides, about 4.4 times the 5 + 5 ball. Total search wall time about 36 min over 5 resumed calls (4 cores), verification 151 s. Replay with the independent 113 verifier (`toolsmith_c1_chain.py replay 8 0 0 ...`): VERDICT YES, failing tests none, meet iso True, no parallel arrows at the meet.
| item | value |
|---|---|
| chain | R4 F1 R7 F1 F4 F3 | F1 F3 R7 R4 R3 F5 (12 edges) |
| meet | isomorphic, parallel arrows: none |
| failing step F4 | Hom(T,T[+-1]) != 0, dim Hom(T,T[-1]) = 1 (must fail, does) |
| controls | 78 corruptions (side B): 77 rejected, 1 true-YES alternative; 12/12 end swaps rejected |

## Reproduction
```
timeout 10m .venv/bin/python -u workshop/rounds/069/skeptic_c1_13.py 1 540 23000   # twice, ~9 min each; the /tmp checkpoint had vanished
timeout 10m .venv/bin/python -u workshop/rounds/114/toolsmith_s8par.py 8 540        # rerun while it prints PAUSED; 5 calls; last one 151 s and prints the chain
timeout 10m .venv/bin/python -u workshop/rounds/113/toolsmith_c1_chain.py replay 8 0 0 "R4 F1 R7 F1 F4 F3" "F1 F3 R7 R4 R3 F5"   # about 2 min
PART=1 timeout 10m .venv/bin/python -u workshop/rounds/114/toolsmith_s8ctl.py   # also PART=swap; PART=0 times out (> 10 min)
```
Logs: `toolsmith_s8_search_p{1,2,3}.txt`, `toolsmith_s8ctl_B.txt`, `toolsmith_s8ctl_swap.txt`. `toolsmith_s8look.py`, `toolsmith_s8prof.py` are the inspection and profile helpers.

## Prior record
E-297 (12 of 13, step 8 undecided), E-298 (controls), E-201, E-168, E-295. Not recorded: a chain for step 8. This closes the "13 of 13" gap of E-297 on the YES side; the thread's real blocker (a NO-capable invariant, n = 9 hit the cap, E-212) is unchanged.

## Code changed
New only: `workshop/rounds/114/toolsmith_s8par.py` (parallel resumable search; reuses 113 code through exec), `toolsmith_s8ctl.py`, `toolsmith_s8look.py`, `toolsmith_s8prof.py`. No library change, no tests needed.
