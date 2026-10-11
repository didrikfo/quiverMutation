# Yes: guarded walks from LNAs reach gate-admitted, non-tilting parents at n = 6, 7, 8, 9, within 5-8 steps, and the Coxeter guard refuses every one of them

author: scholar · round: 014 · kind: result
thread: T5 · bears on: H-015, E-032, E-057, E-059, E-068, E-080, F-038, STEERING round 011 q1 (`isTilting`)

## Claim

Walk (gate = `mutationIsPossibleAtVertex`, guard = Coxeter key of the reduced child = parent's key) from the LNAs
and relation duals of length n, one Coxeter-key class at a time. Reached parents with a vertex k where the gate
admits, `tiltingPlus` (AI 2.32(b) = Ladkani 2.3(c)) is False, exist at **n = 6 (class 0, parents at distance 8),
and in every one of the 10 smallest classes sampled at n = 7, 8, 9 (first parents at distance 5..7)**. All of them
have the A5 shape (a commutative square into d, d with one outgoing arrow, a relation through d>e). **In every
rejection the guard also refuses the step** (key of the child moves); among about 1.29e6 guard-admitted steps
none fails `tiltingPlus`. So the answer to the round question is yes, and it is the gate alone that is unsound at
n <= 9; the guarded walk never takes such a step. n = 5 is closed (11,700 distinct algebras, 2 classes): no
rejection, and A5's own key is on no n = 5 LNA class. It does not claim: the whole n = 9 class space (only the 4 smallest
of 19 classes, to the first rejecting level), a count of rejections per n (stopped at the first rejecting
level), or that `tiltingPlus` False and "guard refuses" are the same thing (see the n = 8 line below).

## Evidence

Size first. LNAs plus duals: n = 5: 28, 6: 84, 7: 264, 8: 858, 9: 2860, in Coxeter-key classes of 2..32 (n = 6) up
to 2..754 starts (n = 9; 19 classes). BFS over distinct algebras (`fingerprint.canonicalKey`), 15-30 ms per
algebra. Closure: n = 5 closes at depth 14 / 16 in 17 s; n = 6 classes had 81k-125k algebras after 540 s and were
still open (new level about x1.5). So full closure is out of reach for n >= 6; the n = 7..9 runs are
stop-at-first-rejecting-level, 540 s budget, 4 in parallel.

| n | class (starts) | distinct algebras | first rejecting parents at distance | rejections at that level | guard-admitted steps seen (tilt True) |
|---|---|---|---|---|---|
| 5 | 0 (12), 1 (16) | 6240, 5460, closed | none | 0 | 16,620 + 13,680 |
| 6 | 0 (2) | 5,616 (to first level); 80,906 in the 540 s run | 8 | 8 (2,217 by depth 14, all guard-refused) | 9,993 |
| 6 | 1, 2, 3 (24, 26, 32) | 119k, 119k, 125k, open | none to depth 11, 12, 15 | 0 | 258k, 274k, 334k |
| 7 | 0, 1, 2 (8, 12, 14) | 3,019, 12,033, 16,784 | 5, 6, 6 | 5, 4, 5 | 4.7k, 18k, 28k |
| 8 | 0, 1, 2 (2, 8, 18) | 2,255, 13,769, 54,326 | 6, 6, 7 | 4, 5, 1 | 3.6k, 22k, 89k |
| 9 | 0, 1, 2, 3 (2, 2, 12, 12) | 1,399, 4,871, 17,241, 17,193 | 5, 6, 6, 6 | 5, 4, 5, 5 | 2.4k, 7.8k, 31k, 30k |

(Distance = the level at which the failing parent is expanded; the parents are found at their first BFS level, so
distance is the shortest in the class up to the canonical-key quotient.) Replay, independent of the BFS (same
`tiltingPlus` and gate code): `scholar_replay.py` rebuilds three recorded paths from their start step by step,
asserting every prefix step is gate-admitted with the key unchanged, then at the last vertex gate True,
`tiltingPlus` False, child key != base. n = 9 class 0: start = LNA with relations `1234, 3456, 4567, 6789`, 5
guarded steps (mutate at vertices 1,1,4,3,1), parent has 9 vertices, relations `[3,1,6,7] = [3,4,6,7]`,
`[5,3,1,6] = [5,3,4,6]`, only arrow out of 6 is 6>7. Same pattern replayed for n = 8 class 2 (7 steps) and n = 7
class 0 (5 steps). Each ends with child key != base.

Why E-057/E-059 found none: E-057 stopped at n = 6 depth 6, n = 7 depth 4; E-059 at n = 6 depth 3, n = 7 depth 2;
the first rejecting parents are at distance 8 (n = 6) and 5 (n = 7), just past those depths. Their statement
"the guard never fired, gate = guard here" holds only to those depths.

A5 itself (E-080) is not the reached algebra. Coxeter-key check (`scholar_a5key.py`): A5 and its pad-before-a
versions (pre = n - 5) have a key on no LNA class at n = 5..8; the other pads have LNA-class keys at n = 6..9, so the key
neither excludes nor shows reach for them. The parents that are reached carry an extra relation (square
commutes after a prefix as well) and are not the E-080 algebras up to padding.

One loose end, not pursued: n = 8 class 2 has 10 steps with gate True, `tiltingPlus` True and key moved
(`M` lines of `scholar_walk_n8_c2.txt`; all parents have parallel arrows or 3-term / duplicated relations), the
reverse of the rejections. Either `tiltingPlus` is not sufficient with parallel arrows or the rewrite is wrong
there; the Cartan congruence was not computed in this script.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/014/scholar_walk.py 9 --plan             # class sizes, 90 s
timeout 10m .venv/bin/python workshop/rounds/014/scholar_walk.py 6 --class 0 --stop-on-reject   # 25 s
workshop/rounds/014/scholar_sweep.sh "7:0 7:1 7:2 8:0 8:1 8:2 9:0 9:1 9:2 9:3"      # 4 at a time, about 12 min
timeout 10m .venv/bin/python workshop/rounds/014/scholar_replay.py 9 0                # 20 s; also "8 2", "7 0"
timeout 10m .venv/bin/python workshop/rounds/014/scholar_a5key.py 9                   # 40 s
.venv/bin/python -m pytest -q tests/test_gate_without_tilting.py                      # 1 s
```
Raw outputs: `scholar_walk_n{5..9}_c*.txt`, `scholar_replay_out.txt`, `scholar_a5key_out.txt`, `scholar_plan_n*.txt`.
Note: `scholar_walk_n6_c0.txt` is from the 540 s run before path recording was added (no path columns).

## Prior record

E-080 (hand-built, "not shown reachable"), E-068 (n = 10 parent reached by a guarded 6-step walk then refused at
step 7; E-032 "ALARM"), E-057/E-059 (none at n <= 7, shallow), F-038 (guard measured at n = 7: 97 quivers whose polynomial
moved). New: gate-admitted rejections occur from n = 6, from LNAs, and are all guard-refused; so E-068 is not special
to n = 10 or to a long walk. Not new: that the guard refuses a gate-admitted non-tilting step (that is what F-038/R-005 built it for).

What it bounds: at n <= 9, to the depths and classes above, the guard is sufficient for `tiltingPlus`: 0 of 1.29e6
guard-admitted steps fail it. The gate alone is not (rejections in 11 of the 14 classes sampled at n = 6..9; the other 3 are n = 6
classes 1-3, open, none to depth 11-15). So the gap between gate and Ladkani's criterion is real and walkable, and it is exactly where the guard
does its work; `isTilting` in the library would only matter for a walk without the guard.

## Code changed

New: `tests/test_gate_without_tilting.py` (A5: gate admits d, dim e_a A e_d = 2, dim e_a A e_e = 1, so the map to
the single arrow out of d is not injective; control abd = acd: 1 -> 1). Uses only library functions, no `isTilting`.
Ran `tests/test_gate_without_tilting.py` (2 passed). Scripts: `scholar_walk.py`, `scholar_sweep.sh`, `scholar_replay.py`,
`scholar_a5key.py`.

## Next

- Chair: STEERING r011 q1 condition (a gate-admitted rejection from a guarded walk from an LNA) is met at n = 6..9.
  Decide `isTilting`; it is redundant under the guard (all rejections are guard-refused, key moves).
- Theorist: the 10 gate+tilt+key-moved steps (n = 8) with parallel arrows; is `tiltingPlus` complete there?
- Toolsmith: n = 6 class 0 closure (80k at 540 s, still growing) and a per-class depth-first rejection depth for the other 15 classes at n = 9
  would be about 10 min per class; I would not run it overnight. Not needed unless someone doubts the sample.
