# H1 for floating rules of width 9..11 holds at lengths w+1..14, 0 failures in 1 692 applications (plus an incomplete length 15)

author: maverick · round: 079 · kind: negative
thread: T4 · bears on: H-020 (H1 of E-159)
scope: floating `VERIFIED_MOVES` only, widths 9, 10, 11 (18 rules each, 54 total); every rule at every length L with w < L <= 14 (complete); L = 15 partial (sharded runs killed at 9 min). Anchored rules, lengths >= 16, L = 15 remainder, engine independence not covered; widths 4..8 not re-checked here. Weakest cell: w11 L14 (90 applications, 18 rules). Speculation level: supported (exhaustive within cells, but same `verifyMove` code as E-159/E-165/E-188).

## Response to referee

Verdict was minor revision (toolsmith, review `maverick.review.md`). Reproduction, arithmetic and code path were confirmed; no numerical error found.

1. Claim wording: done. "Several lengths beyond w+2 where the width allows" is replaced by the referee's narrowed wording (widths 9..11 at w+1..14; new cells named; 4..8 unchanged from E-165/E-188/E-159). The Claim no longer says anything about widths 4..8.
2. Weakest cell: done. The Claim and the Scope line now state that w11 at L14 is 90 applications over 18 rules (5 per rule on average) and is the weakest cell. The Observation under Evidence is unchanged.
3. Not independently re-run: noted. The w9 L14 (756 applications, 280 s) and w10 L14 (252 applications, 239 s) cells were not independently re-run in review; the log lines in `maverick_wide_log.txt` are the only record of them. All cells L <= 13 and w11 L14 were reproduced by the referee.
4. Referee's remarks without a required item: the L15 per-shard progress files were not kept and are only progress, not counted (unchanged); the Prior record does not claim a search for an older "w+1/w+2 only" statement beyond E-165's scope line (the referee noted the same).

Title unchanged: it already matches the evidence (w+1..14, 1 692 applications, 0 failures).

## Claim

No counterexample: each of the 54 floating rules of window width 9..11 passes `verifyMove` (predicted LNA, every step admissible, Coxeter polynomial unchanged) at every LNA length from w+1 to 14, so widths 9..11 are now checked at lengths w+1..14 (new: w9 at 12..14, w10 at 13..14, w11 at 14); widths 6..8 and 4..5 are unchanged from E-165/E-188/E-159. Weakest cell: w11 at L = 14 is 90 applications over 18 rules (5 per rule on average), so the w11 evidence beyond w+2 rests on that one thin cell. It does NOT claim the rules hold at all lengths (a theorem), nor says anything about anchored rules, nor length 15 as a whole. A single failing application at any larger length refutes the "length independence" reading of H1 for that rule.

## Evidence

Applications confirmed (every rule, no cap), failures 0 in every cell:

| L \ w | 9 | 10 | 11 |
|---|---|---|---|
| 10 | 18 | - | - |
| 11 | 36 | 18 | - |
| 12 | 90 | 36 | 18 |
| 13 | 252 | 90 | 36 |
| 14 | 756 | 252 | 90 |

(w+1 is the smallest length; cells L <= w not applicable.) Total 1 692 applications in complete cells; the rule count per width is 18 in every cell. Note E-165's `lengthsToCheck` only covered w+1..w+2 for w >= 7, i.e. w = 9: 10, 11; w = 10: 11, 12; w = 11: 12, 13. So new beyond that are: w9 at 12, 13, 14; w10 at 13, 14; w11 at 14 (1 530 applications of the total). The count of applications grows ~3x per length, so one length is ~4x the previous in time (w9: 14 s, 70 s, 280 s at 12, 13, 14).

Length 15 (partial, 4-way shards, 9 min timeout, 12 processes on 4 cores): applications before the kill, 0 failures: w9 660 (5 of 18 rules; shards 1 and 2 never finished a rule), w10 378 (9 of 18 rules), w11 224 (16 of 18 rules, last rule of two shards cut). Those are not cells and are not counted as checked; they are in `maverick_wide_log.txt` only as progress.

Observation: the application count at a fixed length is largest for the narrowest width and shrinks as w grows towards L, so widths 9..11 are checked on far fewer LNAs than widths 6..8 at the same L; a rule of width 11 at L = 14 has 90 applications, i.e. few windows. This is why the widths-9..11 check is weaker per cell than the table suggests.

## Reproduction

```
for w in 9 10 11; do for L in 10 11 12 13 14; do [ $L -gt $w ] && timeout 9m .venv/bin/python workshop/rounds/079/maverick_wide.py $L $w | tail -1; done; done
```
Cells L <= 13: under 80 s each; L = 14: 280 s (w9), 239 s (w10), 202 s (w11), three in parallel. Length 15 shard form: `timeout 9m .venv/bin/python workshop/rounds/079/maverick_wide.py 15 W S 4` for S in 0..3; needs about 15-25 min per width on one core, so a full cell is an `OVERNIGHT.md` item (see Next). Log of the DONE lines: `workshop/rounds/079/maverick_wide_log.txt`.

## Prior record

E-159 (w <= 5), E-165 (w 6..8 at one length; 54 rules of width 9..11 beyond w+2 named unchecked), E-188 (6..8 gaps), E-207 (second implementation, w 6..8; widths 9..11 named open). `maverick_wide.py` is E-188's script (`experimentalist_gaps.py`) copied with the same filter; so the same limitation: shared package matching and engine. Nothing in RETRACTIONS touched. New: widths 9..11 at lengths past w+2.

## Code changed

None (new script `workshop/rounds/079/maverick_wide.py`, a copy of the round 062 script; no tests needed).

## Next

- Overnight proposal: L = 15 for widths 9, 10, 11 in full (est. 3 x 4 shards, about 25 min core-time per width-shard set of 4 for w9; the w9 and w10 cells are the large ones) and L = 16 for w = 11 only.
- Skeptic: is there any structural reason the rule widths stop at 11 (18 rules per width, 54 total)? A proof of H1 for rule families (the 18 look parametric) would replace all of this length sampling; theorist.
- Anchored rules still unchecked above w+2; toolsmith could extend `lengthsToCheck` to w+4 for 9..11 (L <= 14 is now known cheap).
