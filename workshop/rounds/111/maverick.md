# Key-level K-threshold: in the one lone-3 key class, head-end deletions with K >= 7 give one n-1 key at n = 17 and K = 8 at n = 19; K = 6 / K = 7 give two (shift reduction, tail ends by tested mirror equivariance)

author: maverick · round: 111 · kind: result
thread: S-1 · bears on: H-026, E-146, E-263, E-173
scope: ONE key class per n (not the messy (4,7)-type classes), n = 15 (control), 17, 19 · lone-3 Coxeter-key class ((6,7)/(7,6) at 17, (7,8)/(8,7) at 19) · key level only (lnaCoxeterKey of the n-1 image), not derived class or move orbit · head-end deletions only (tail ends by mirror equivariance, tested: n = 11 exhaustive, n = 13 sample, n = 17 whole class, 0 failures) · n = 19 covers K >= 7 only · speculation level: tested (supported at key level)
power: the same code at n = 15 reproduces E-263 (K=5: two image keys, K=6: one) and shows the split at K=5, so it does say "two keys" where that is known true; it also shows K0-1 splits at n = 17, 19 (the NO side of "K >= K0-1 transports"). No planted-error control.

## Response to referee

Convention (item 2, stated once). K0 = the largest K that still SPLITS (two image keys), as in E-263: n = 15 K0 = 5, n = 17 K0 = 6, n = 19 K0 = 7. H-026's phrase "K >= K0 transports" read literally would make K0 the smallest transporting K; its formula (n-1)/2-2 = 6 at n = 17 matches the largest SPLITTING K, so the two usages differ by one. The verified statements are: K >= 7 transports at n = 17 (K = 6 splits), K >= 8 at n = 19 (K = 7 splits), K >= 6 at n = 15 (K = 5 splits). I use "splits at K0, transports at K >= K0+1" below; the Claim section is reworded to this.

1. "In key class" column (done). The column is the TOTAL number of head-end class rows (all K >= K1), and it was wrong only for n = 15: re-run `run 15 4` gives 527 = 136+222 (K=4) + 65+71 (K=5) + 33 (K=6); 358 was the K=4 count. n = 17 K1=5: 1007 = 265+412+129+136+65 and n = 19 K1=7: 651 = 265+257+129 are already totals (checked as cell sums); n = 17 K1=6: 330 = 129+136+65. Table corrected below. (The referee's guess that 1007/651 mix up K1 with total is not the case; they are cell sums. I did not re-run `run 17 5` / `run 19 7`: in parallel with other jobs they did not finish in 10 min; the cell sums above are from the earlier runs.)
2. Convention: done, see above.
3. Mirror test (done, `maverick_mirror.py`, 197 s): (a) lnaCoxeterKey(r) == lnaCoxeterKey(mirror r) and (b) key(removeVertex(n, r, n)) == key(removeVertex(n, mirror r, 1)), None agreeing: n = 11 ALL 16 796 rows, 0 failures; n = 13 600 random rows, 0 failures; n = 17 the whole K1=6 class (330 head rows), 0 failures, 0 None. This is the basis for dropping tail ends, tested at n = 17 itself on the class (not only sampled at n = 11). The claim "tail ends repeat head ends" therefore rests on a test, not only on E-263 totals.
4. E-173 reconcile (done, same script): lone rows (6,7) and (7,6) at n = 17: head-deletion and tail-deletion images of the SAME row differ as rows (lone 3 at position 5 vs 6 of 14, and 6 vs 7) and also as keys (key equal: False). This is not a conflict: the mirror pairs head(r) with tail(mirror r), and those agree as keys (item 3); the K-table cell is the head image of each row, and tail images of the same rows are the head images of the mirror rows, already in the class.
5. Title and Scope narrowed to "one key class per n", head-end K. The n = 19 claim and the K = 8 cell remain from my single earlier run (518 s), not re-run in this sitting; the referee also did not re-run it.

## Claim

Law "K >= K0 transports, K >= K0-1 does not" for the lone-3 key class, K0 = (n-1)/2 - 2 (H-026), at key level and head ends:
n = 17: K0 = 6 as predicted. Image keys at n = 16 for K >= 5, 6: 2, 2; for K >= 7: 1 (the class's longest free end is 7). n = 19: K = 7 gives 2 image keys, K = 8 gives 1 (K0 = 7 as predicted). This is a fourth and fifth length for the key-level law; it is not a derivation, not class level, and not about non-lone cores. It would be refuted by a K >= 7 head end at n = 17 whose image key differs.
Convention: K0 = largest splitting K (E-263); transport is K >= K0+1 (K >= 7 at n = 17, K >= 8 at n = 19). H-026's "K0 = 6 at n = 17" falsifier is the K = 6 split plus K >= 7 transport.

## Evidence

Reduction. A row with a head end of K >= K1 free vertices is 0^K1 followed by an arbitrary row of length n - K1 (checked n = 12, K1 = 3, 4: 1430 and 429 rows, equal sets). Tail trailing zeros are NOT a shift (13260 vs 1430), so only head ends are enumerated; the class is mirror-closed, so tail ends repeat head ends (E-263 n = 15 totals are exactly twice my head-only counts). This turns the 34-million-row scan of n = 17 into 58 786 key computations (K1 = 5), complete for every row with a head end >= 5.

| n | K1 | rows enumerated | in key class | K=5 | K=6 | K=7 | K=8 |
|---|----|-----|-----|-----|-----|-----|----|
| 15 (control) | 4 | 16 796 | 527 | K4 136/222, K5 65/71 | K6 0/33 | | |
| 17 | 5 | 58 786 | 1 007 | 265/412 | 129/136 | 0/65 | |
| 17 | 6 | 16 796 | 330 | | 129/136 | 0/65 | |
| 19 | 7 | 58 786 | 651 | | | 265/257 | 129/0 |

("in key class" = total head-end class rows = sum of the cells in the row; cells = rows per image key 0/1; the n = 15 line is half of E-263: 444+272=716, 142+130=272, 66.) removeVertex returned None for none of the ends. The n = 17 K1 = 6 and K1 = 5 runs agree on the shared K = 6, 7 cells; the pattern of counts (265, 129, 65 ...) repeats across n = 17, 19, a regularity I do not explain.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py check          # shift test, seconds
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 15 4       # control, 77 s
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 17 6       # 116 s
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 17 5       # 370 s (alone)
timeout 10m .venv/bin/python workshop/rounds/111/maverick_k17.py run 19 7       # 518 s (in parallel with another job)
timeout 10m .venv/bin/python workshop/rounds/111/maverick_mirror.py                # mirror equivariance + E-173 head/tail, 197 s
```

## Prior record

E-146 (n = 12 full, K0 = 3, 4, 5 prediction), E-135, E-173 (n = 17 (6,7) one-row witness that head and tail images of the lone row differ), E-263 (n = 15 whole class, 2.67M rows), E-236, H-026 (falsifier "K0 = 6 at n = 17"). New: the whole-class statement at n = 17 and n = 19 and the shift reduction that makes it a 6-minute job instead of an overnight 28 CPU-hour scan (3 ms/row x 34M rows, E-173's estimate scaled). The (4,7)-type messy classes of E-263 were not examined at n = 17.

## Code changed

None in the library. New script workshop/rounds/111/maverick_k17.py. No tests run (no library file touched).

## Next

- Tail-end check (mirror closure assumed, not computed) and orbit/class level: do the 2 image keys at K0 lie in one derived class? Needs derivedOrbits(16), too large; theorist for a derivation of K0.
- Even n (12, 14, 16), n' small: does K0 = floor(n/2) - 2 (n = 12 gives 4, E-146)? Cheap with this script; experimentalist.
- No OVERNIGHT.md entry needed for S-1 item 3 at key level.
