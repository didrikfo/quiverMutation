# The K-threshold law K0 = 6 holds at n = 17 and K0 = 7 at n = 19 (key level, lone-3 key class, head ends, via a shift reduction)

author: maverick · round: 111 · kind: result
thread: S-1 · bears on: H-026, E-146, E-263, E-173
scope: n = 15 (control), 17, 19 · lone-3 Coxeter-key class ((6,7)/(7,6) at 17, (7,8)/(8,7) at 19) · key level only (lnaCoxeterKey of the n-1 image), not derived class or move orbit · head-end deletions only (tail ends by mirror closure, not computed) · n = 19 covers K >= 7 only · speculation level: tested (supported at key level)
power: the same code at n = 15 reproduces E-263 (K=5: two image keys, K=6: one) and shows the split at K=5, so it does say "two keys" where that is known true; it also shows K0-1 splits at n = 17, 19 (the NO side of "K >= K0-1 transports"). No planted-error control.

## Claim

Law "K >= K0 transports, K >= K0-1 does not" for the lone-3 key class, K0 = (n-1)/2 - 2 (H-026), at key level and head ends:
n = 17: K0 = 6 as predicted. Image keys at n = 16 for K >= 5, 6: 2, 2; for K >= 7: 1 (the class's longest free end is 7). n = 19: K = 7 gives 2 image keys, K = 8 gives 1 (K0 = 7 as predicted). This is a fourth and fifth length for the key-level law; it is not a derivation, not class level, and not about non-lone cores. It would be refuted by a K >= 7 head end at n = 17 whose image key differs.
Note the wording: here K0 is the largest K that still splits (E-263's usage), so "K >= 7 transports" at n = 17.

## Evidence

Reduction. A row with a head end of K >= K1 free vertices is 0^K1 followed by an arbitrary row of length n - K1 (checked n = 12, K1 = 3, 4: 1430 and 429 rows, equal sets). Tail trailing zeros are NOT a shift (13260 vs 1430), so only head ends are enumerated; the class is mirror-closed, so tail ends repeat head ends (E-263 n = 15 totals are exactly twice my head-only counts). This turns the 34-million-row scan of n = 17 into 58 786 key computations (K1 = 5), complete for every row with a head end >= 5.

| n | K1 | rows enumerated | in key class | K=5 | K=6 | K=7 | K=8 |
|---|----|-----|-----|-----|-----|-----|----|
| 15 (control) | 4 | 16 796 | 358 head ends | K4 136/222, K5 65/71 | K6 0/33 | | |
| 17 | 5 | 58 786 | 1 007 | 265/412 | 129/136 | 0/65 | |
| 17 | 6 | 16 796 | 330 | | 129/136 | 0/65 | |
| 19 | 7 | 58 786 | 651 | | | 265/257 | 129/0 |

(cells = rows per image key 0/1; the n = 15 line is half of E-263: 444+272=716, 142+130=272, 66.) removeVertex returned None for none of the ends. The n = 17 K1 = 6 and K1 = 5 runs agree on the shared K = 6, 7 cells; the pattern of counts (265, 129, 65 ...) repeats across n = 17, 19, a regularity I do not explain.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py check          # shift test, seconds
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 15 4       # control, 77 s
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 17 6       # 116 s
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 17 5       # 370 s (alone)
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 19 7       # 518 s (in parallel with another job)
```

## Prior record

E-146 (n = 12 full, K0 = 3, 4, 5 prediction), E-135, E-173 (n = 17 (6,7) one-row witness that head and tail images of the lone row differ), E-263 (n = 15 whole class, 2.67M rows), E-236, H-026 (falsifier "K0 = 6 at n = 17"). New: the whole-class statement at n = 17 and n = 19 and the shift reduction that makes it a 6-minute job instead of an overnight 28 CPU-hour scan (3 ms/row x 34M rows, E-173's estimate scaled). The (4,7)-type messy classes of E-263 were not examined at n = 17.

## Code changed

None in the library. New script workshop/rounds/111/maverick_k17.py. No tests run (no library file touched).

## Next

- Tail-end check (mirror closure assumed, not computed) and orbit/class level: do the 2 image keys at K0 lie in one derived class? Needs derivedOrbits(16), too large; theorist for a derivation of K0.
- Even n (12, 14, 16), n' small: does K0 = floor(n/2) - 2 (n = 12 gives 4, E-146)? Cheap with this script; experimentalist.
- No OVERNIGHT.md entry needed for S-1 item 3 at key level.
