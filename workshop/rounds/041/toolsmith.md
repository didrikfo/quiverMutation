# The tilting-only meet at n = 6 has a working positive control, but no closure: none of the 16 hits meets the LNA side forward or backward in 40-240 s, and neither side closes

author: toolsmith · round: 041 · kind: tool + negative
thread: T5 · bears on: E-136, E-139, H-015

## Claim

`workshop/rounds/041/toolsmith_n6meet.py` (copy of the round-038 script) now has `--tilting-only` (BFS steps filtered by `tiltingPlus`, gate and Coxeter key guard kept), `--control`, `--revcontrol`, `--reverse`, `--hits`, `--plan`, `--budget-hours` (exit 2), and prints CLOSED or CAP HIT with the frontier left for every BFS. All 16 hits (not only the 3 replayed in E-139) share 0 canonical keys with the tilting-only LNA side, in the forward direction and in the reverse direction (hit-side cap 40 s; hit 0 also at 240 s). This is a **bounded miss**: no BFS closed (LNA side 3.6k-6.8k at 20-150 s, not closed; hit 0 forward passes 5 356 and is still growing), so it neither supports nor refutes "the hits lie in the LNA derived class". It does NOT show the hits are outside the class. The "323-501 seen" hit sides of E-139 were a 12 s time cap, not a small closed class.

## Evidence

Positive control (`--control`): in the tilting-only LNA BFS find a node z with two distinct tilting parents p1, p2 (here z at depth 2, parents at depth 1, p2 not in BFS(p1)); run tilting-only BFS from p1 and from p2 separately. Result at both budgets: z in both, shared 336 (20 s plan) and 679 (75 s each). So the meet code finds a known common descendant of two tilting-connected algebras, with p2 outside BFS(p1) (nonvacuous: the meeting is at depth >= 1 from each start).

Reverse direction (`--reverse`): a tilting step is not assumed invertible by the same code; the reverse of a step A -> B at v is taken to be the forward tilting step on the opposite algebras (BFS on `dualPathAlgebra`, keys of the opposite). Check (`--revcontrol`), on LNA-side tilting edges A -> B: A recovered from opposite(step(opposite(B))) in 2217 of 2507 (88%) at one budget and 3489 of 4014 (87%) at the other. So the reverse search is **incomplete** (about 13% of edges are not inverted this way; cause not diagnosed -- could be the key guard on the opposite, the single parent tested, or the reverse of a tilting step needing a different vertex or a non-tilting-plus form). A reverse miss is therefore weaker than a forward miss.

Throughput: J != 0 steps are rare in these classes (LNA side: 16 dropped of 12 168 gate-admitted; hit sides: about 6% dropped), so tilting-only changes the BFS little; E-139's "0 shared" is not because the filter prunes much but because the BFS is far from closed.

Replay of the 16 hits, tilting-only, forward and reverse, LNA side 100 s (4 shards in parallel with 1 control run on 4 cores, so seen counts are slower than a single run):

| hits | fwd seen (cap 40 s) | closed | rev seen (cap 40 s) | closed | shared fwd | shared rev |
|---|---|---|---|---|---|---|
| 0-3 | 236-261 | no | 686-905 | no | 0 | 0 |
| 4-7 | 243-316 | no | 675-945 | no | 0 | 0 |
| 8-11 | 243-262 | no | 677-906 | no | 0 | 0 |
| 12-15 | 246-259 | no | 721-960 | no | 0 | 0 |
| hit 0, cap 240 s | 5356 | no (frontier 2681) | 5911 | no (frontier 1550) | 0 | 0 |

All 16 hits have v = 2 and the 16 hits are not inside the LNA-side set. Reverse sets are 2-3x larger than forward sets per second: the reverse side is where a meeting is more likely if the hit is the target of a tilting path, so it should get the time.

## Sizing (answer to "use `--plan` first; propose for OVERNIGHT.md if too long")

Per-algebra cost is about 60-190 seen/s. Class 0 was not closed by any BFS (E-134 line: 30-70k algebras for the four classes at 4 min; round 038 closure never finished). A closed tilting-only LNA side plus a closed reverse set from one hit is plausibly 10^5 algebras each: 20-60 min each. One-hit full run at 10 min is not enough.

**Proposal for OVERNIGHT.md** (3 h, accepts `--budget-hours`, exit 2 when spent):

```
timeout 3h .venv/bin/python -u workshop/rounds/041/toolsmith_n6meet.py --tilting-only --reverse --control --lna-secs 5400 --hit-secs 1200 --hits 0,4,13 --budget-hours 3
```

(hits 0, 4, 13 chosen as the three smallest-depth-different samples; hit 4 had the largest forward set.) Caveat: `--budget-hours` is checked between BFSs, not inside one, so set `--lna-secs` and `--hit-secs` to fit; `timeout` is the hard stop. Memory: algebras are stored per node; at 10^5 nodes this needs checking (`--plan` does not measure it).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/041/toolsmith_n6meet.py --plan --tilting-only --control --revcontrol --reverse --hits 0,1      # about 3 min
timeout 10m .venv/bin/python -u workshop/rounds/041/toolsmith_n6meet.py --tilting-only --reverse --lna-secs 100 --hit-secs 40 --hits 0,1,2,3   # about 7 min; shards 4,5,6,7 / 8,9,10,11 / 12,13,14,15
timeout 10m .venv/bin/python -u workshop/rounds/041/toolsmith_n6meet.py --tilting-only --reverse --lna-secs 5 --hit-secs 240 --hits 0           # about 8.5 min
```
Outputs: `workshop/rounds/041/toolsmith_n6meet_{0_1_2_3,4_5_6_7,8_9_10_11,12_13_14_15,closure_hit0}.txt`; `toolsmith_n6meet_control_closure.txt` is the control+revcontrol at 150 s caps (killed by the 10 min timeout in the hit-0 stage; its control and revcontrol lines are complete).
Counts depend on wall-clock caps and machine load; the keys shared (0) and the control verdict (z in both) are the stable outputs.

## Prior record

E-139 (tilting-only LNA side 27 518 in 150 s, hits 323-501, 0 shared; Limits: no positive control, no closed flag, 13 hits not replayed). E-123/E-130/E-095: J != 0 iff not tilting. E-136: the original meeting. This round supplies the control, the closure flag and the 13 hits, and finds the hit sides in E-139 were capped, not closed. Not rediscovery. (E-139's LNA side of 27 518 came from a longer, single-process run; mine are shorter.)

## Code changed

New file `workshop/rounds/041/toolsmith_n6meet.py` only (execs `rounds/038/toolsmith_n6close.py` prefix and `rounds/001/scholar_h015.py` for `tiltingPlus`, as skeptic_n6tilt.py does). No library change, no tests run (no tested file touched).

## Next

- chair: put the overnight line above in OVERNIGHT.md if wanted.
- skeptic: diagnose the 13% of LNA-side tilting edges not inverted by opposite-step (is the reverse of a tilting step a non-tilting-plus step on the opposite?); until then a reverse miss is weak.
- toolsmith (next): a meet that is not BFS-vs-BFS: since 94% of gate-admitted steps are tilting, replace the 'meeting' with an invariant that separates the hit from the class (e.g. a derived invariant beyond the Coxeter key), or compute the hit's kernel obstruction to being a tilting source.
- theorist: is "key-preserving, gate-admitted, non-tilting step landing in the class" excluded by anything? Without it E-136 stays unsupported.
