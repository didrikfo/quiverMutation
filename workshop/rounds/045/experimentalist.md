# The 10 F-041 n = 8 merges also hold in a J = 0 walk (implied by gate = J=0 = key kept on all 2396 depth-3 edges); guarded depth-4 LNA-started forward edges are all J = 0 (n = 6, 7 full; n = 8 first 244 of 429 LNAs)

author: experimentalist · round: 045 · kind: negative
thread: T5 (T10 guard audit b) · bears on: F-041, E-037/E-038, E-145, E-084, H-012, H-015
scope: n = 8, the 10 single two-arrow deletions left open by `lm.ALL_MOVES` with `free=False, edges=True, doubles=True` (the F-041 definition; a different move set gives a different list), meeting depth 3 per side (+ relation duals). Plus edge tally of guarded depth-4 forward BFS from LNAs only (not from non-LNA starts, not on the dual): n = 6 (42/42 LNAs), n = 7 (132/132), n = 8 (first 244 of 429 LNAs in enumeration order, 250 s cap). NOT covered: pipeline merges at depth 5-6 (merges.py/classify.py), E-094-style deep walks (distance >= 5), F-036 reflection bridges, n = 9.

## Response to referee

Verdict received: minor revision. All four required items done.

1. Library cross-check, all 10 pairs: DONE. `experimentalist_crosscheck.py` (output `experimentalist_crosscheck.txt`, ~35 s) runs `search.meetingPoints(A, B, 3, alsoDual=True)` and the author's `reach(..., 'guard')` on every pair. Each pair: library 1 meeting, split (3, 3); author 1 meeting, split (3, 3); identical meeting-key set on 10 of 10 pairs. This closes the "Next: Skeptic" item (the reimplementation agrees with the library on all 10, not 2 of 10).
2. Headline reworded: DONE. "Meets in every mode" is implied by the edge identity (gate = J=0 = key kept on all 2396 edges, so the three walks are the same graph). It is a consequence, not an independent result; the tilt-mode runs add no information beyond that identity. Title changed.
3. Dual side: DONE, moved into the Claim. J is tested on the algebra actually mutated, which on the dual side is the opposite algebra. "All J = 0" there means the opposite step passes `tiltingPlus`, not that the forward step is a tilting step. The forward-dual correspondence was not checked. "Tilting-only walk" must not be read as covering the forward dual side.
4. Sweep counts: DONE. 56 692 is summed over LNAs and over distinct nodes per BFS, so one algebra is counted once per LNA whose ball contains it; it is not independent evidence (the claim "no J != 0 step" is unaffected). The n = 8 subset is the first 244 of 429 LNAs in enumeration order, cut by a time cap, so a sample.
Also from the Scope section of the review: the "all 10" now states its dependence on the move set; "all LNAs" for n = 6, 7 is narrowed to forward-only walks from LNAs (mergeReport also walks the dual and from non-LNA starts).
Not done: the n = 7 / n = 8 sweeps were not rerun (unchanged, referee reran n = 6 only); the depth-2 control rerun was not repeated (referee did not rerun it, not asked).


## Claim

The class merges recorded in F-041 at n = 8 (10 pairs, joined by `meetingPoints` at depth 3 + 3) do not depend on a J != 0 step: with the walk restricted to gate + `tiltingPlus` (J = 0) and the Coxeter-key guard OFF, all 10 pairs still meet, at the same total 6, and the gate-admitted edge set is identical to the guarded one (2396 of 2396 edges J = 0, key kept). In the broader sample (guarded BFS depth 4 from LNAs, n = 6, 7, 8) all 56 692 gate+key-admitted edges are J = 0; no J != 0 step occurs. This does not say that no recorded merge uses one: E-145's J != 0, key-preserving steps are walk descendants first seen at distance >= 5 (n = 7), beyond every depth walked here, and the deep guarded merges (pipeline depth 5-6, E-094) were not replayed. On the dual side J is tested on the opposite algebra (the one actually mutated), so J = 0 there says the opposite step passes `tiltingPlus`; that it corresponds to a forward tilting step was not checked. The tilt-mode meeting of all 10 is implied by the edge identity above, not a separate finding. The 56 692 count is summed over LNAs (same algebra counted repeatedly) and the n = 8 part is the first 244 of 429 LNAs. It would be refuted by a recorded merge whose tilting-only walk fails to meet.

## Evidence

Pairs (derived with `freeMoves.derivedOrbits(8, ...)`, exactly the test `test_every_single_two_arrow_deletion_at_length_eight_is_a_mutation`; the 10 are reproduced, not assumed). Modes: guard = gate + key guard (recorded search); tilt = gate + `tiltingPlus` only, key guard off, with a per-step check whether the child key equals the base key; both = control.

| pair (n = 8 row -> row) | guard: meets, total | tilt: meets, total | edges gate / J=0 / key kept (tilt) |
|---|---|---|---|
| 230300 -> 030300 | yes, 6 | yes, 6 | 242 / 242 / 242 |
| 230302 -> 030302 | yes, 6 | yes, 6 | 140 / 140 / 140 |
| 230302 -> 230300 | yes, 6 | yes, 6 | 140 / 140 / 140 |
| 230400 -> 030400 | yes, 6 | yes, 6 | 269 / 269 / 269 |
| 240030 -> 040030 | yes, 6 | yes, 6 | 269 / 269 / 269 |
| 250002 -> 050002 | yes, 6 | yes, 6 | 278 / 278 / 278 |
| 250002 -> 250000 | yes, 6 | yes, 6 | 278 / 278 / 278 |
| 304002 -> 304000 | yes, 6 | yes, 6 | 269 / 269 / 269 |
| 400302 -> 400300 | yes, 6 | yes, 6 | 269 / 269 / 269 |
| 030302 -> 030300 | yes, 6 | yes, 6 | 242 / 242 / 242 |

Each pair has exactly 1 meeting (shared labelled quiver) in each mode. Totals: 2396 gate-admitted edges; 0 dropped by `tiltingPlus`; 0 dropped by the key guard; 0 tilting steps that move the key. Control: depth 2 per side gives no meeting for any pair in either mode (matches `test_meeting_in_the_middle_reaches_twice_as_far`: the depth-3 meeting is needed). Mode "both" equals "guard" (identical output).

Edge tally for the guarded walk shape of `mergeReport` (depth 4, gate + key guard, no dual): n = 6: 3 270 edges / 42 LNAs, J != 0: 0; n = 7: 16 126 edges / 132 LNAs, 0; n = 8: 37 296 edges / 244 of 429 LNAs (cap), 0. (Edge counts are over distinct nodes per BFS, summed over LNAs; not independent evidence.)

What this means for T5: at depths <= 4 the Coxeter-key guard and `tiltingPlus` accept the same steps on LNA-started walks, in line with E-084 (0 of ~1.3e6 guard-admitted failures), and E-145's gate-admitted key-preserving J != 0 steps do not appear this shallow. So F-041's n = 8 table and H-012 at n = 8 hold for the tilting-only walk too. On the dual side the J test is applied to the algebra actually mutated (the opposite); I did not check separately that the opposite step is the dual of a forward tilting step (cf. E-142 `--revcontrol`).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/045/experimentalist_crosscheck.py   # library cross-check, 10 pairs, ~35 s
.venv/bin/python workshop/rounds/045/experimentalist_t10b.py --plan                 # lists the 10 pairs, 5 s
timeout 10m .venv/bin/python -u workshop/rounds/045/experimentalist_t10b.py --depth 3   # 3 modes x 10 pairs, ~1 min (output: experimentalist_t10b_d3.txt)
timeout 10m .venv/bin/python -u workshop/rounds/045/experimentalist_t10b.py --depth 2 --modes guard,tilt   # control (..._d2.txt)
timeout 10m .venv/bin/python -u workshop/rounds/045/experimentalist_t10b_sweep.py 7 4 250   # 83 s; n 6: 8 s; n 8: 250 s cap
```

## Prior record

F-041 / E-037 / E-038 (the 10 pairs, guarded); E-084 (guard-admitted steps never fail `tiltingPlus`, ~1.3e6, first rejection at distance 5-8); E-145 (J != 0 key-preserving steps, n = 7 c1, c2, walk descendants). Not recorded before: a tilting-only replay of any recorded merge. The result is the expected one given E-084 and is a negative for "merges depend on J != 0" at this depth, not a surprise.

## Code changed

None to `quivermutation/`. New: `workshop/rounds/045/experimentalist_t10b.py`, `experimentalist_t10b_sweep.py`, `experimentalist_crosscheck.py` (+ `.txt`), outputs `experimentalist_t10b_d2.txt`, `_d3.txt`, `_sweep_n7d4.txt`, `_sweep_n8d4.txt` (all small). No tests run (no library file touched).

## Next

- Done in response: library `meetingPoints` agrees on all 10 pairs (see Response to referee).
- Experimentalist (overnight): the deep recorded merges, E-094's n = 8 c2 depth-8 walk (63 221 algebras) with per-edge `tiltingPlus` tally, and the pipeline depth-5/6 merges at n = 7, 8; there is where J != 0 steps live (distance >= 5).
- Toolsmith: store the witness path of each recorded class merge (merges.py currently does not), so audits can replay instead of re-search.
