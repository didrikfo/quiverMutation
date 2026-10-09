# At n = 15 the lone-3 key class (5,6)/(6,5) already shows the K0 = 5 failure by a one-row witness; the full class scan (2 674 440 rows, ~33 min on 4 cores) is sized but not run

author: toolsmith · round: 055 · kind: result (sizing + witness)
thread: S-1 · bears on: E-144 (K-threshold law), E-139
scope: n = 15 (and 11, 13, 16, 17 as controls) · single-3 cores only, comparison by Coxeter key (not derived class) · the class-level counts (analogue of 136 + 66) and move orbits at n = 15 NOT computed

## Claim

Option (a), cheap half. (1) The full n = 15 key scan of E-144's style does not fit one 10-minute command but shards: 2 674 440 LNA rows, key cost ~3.0 ms/row, so ~8 000 CPU-s, ~33 min wall on 4 cores as 24 shards of ~6 min each (each shard also pays a ~23 s enumeration). (2) It is not needed to decide the K0 = 5 prediction: the lone 3 at (h,K) = (5,6) and (6,5) share a key, and for both, head and tail deletion (via `removeVertex`) give images with different Coxeter keys at n = 14, with the deleted end having K >= 5 in each case. That is the same kind of witness as E-139 (n = 13, K0 = 4) and E-144 (n = 11, K0 = 3). Hence "K >= 5 free ends transport" fails at n = 15 within the lone-3 key class, provided the two lone 3s are in the class, which they are by construction (same key). Does NOT claim: that the class is one move orbit, how the ends split in counts, that K >= 6 transports (no n = 15 row has both h, K >= 6 for a single 3; the control at n = 16 (6,6) is self-mirror and shows equal keys, which is vacuous as in E-144), or anything about derived classes (key only).
Refuted by: a different `removeVertex` convention giving equal keys for those two images (checked by script, they are the library's), or a key collision/cap in `lnaCoxeterKey` (key of lone rows is computed, not capped).

## Evidence

| n | pair (h,K) | same key | head vs tail image keys equal |
|---|---|---|---|
| 11 | (3,4),(4,3) | yes | no, no |
| 13 | (4,5),(5,4) | yes | no, no |
| 15 | (5,6),(6,5) | yes | no, no |
| 15 | (4,7),(7,4) | yes | no, no |
| 17 | (6,7),(7,6) | yes | no, no |
| 16 | (6,6) | -- | yes (vacuous, h = K) |

n = 11, 13 reproduce E-144 / E-139; n = 15 is the predicted K0 = 5; n = 17 would be K0 = 6 (also failing at (6,7): predicted by the pattern K0 = (n-5)/2 for odd n, observed on single-3 rows for n = 11, 13, 15, 17 only). All n = 15 single-3 rows fall in 6 key classes, each a mirror pair {(h,K),(K,h)} (printed by the script).

Sizing (counts: `allRelationLengths`): n = 13 208 012, n = 14 742 900, n = 15 2 674 440 rows (enumeration 1.3 s, 6.2 s, 23.1 s). Key cost 1.68, 2.13, 2.98 ms/row at n = 12, 13, 15 (first 5000 rows of a 20 000-row prefix; one timing each, other personas share the box). n = 12 scan was 118 CPU-s for 58 786 rows, consistent. Class stage (union-find with moves and ends) for a class of maybe 10^4-10^5 rows is untimed; at n = 12 (2746) it took ~5 s, at n = 13 (5023) about 5 s, so probably minutes, but this is an extrapolation.

## Reproduction

```
.venv/bin/python workshop/rounds/055/toolsmith_s1witness.py      # about 3 s
.venv/bin/python workshop/rounds/055/toolsmith_s1plan.py 15      # about 40 s (enumeration 23 s + 5000 keys); also args 12, 13
```
Proposed (not run) full scan: adapt `workshop/rounds/042/maverick_n12.py scan i 24` with n = 15 and targets TG for h = 5 (and 4); 4 concurrent shards x 6 waves, each under `timeout 10m`.

## Prior record

E-144 (n = 12 law, "K0 = 5 at 15 predicted, not run"), E-139 (n = 13 K0 = 4 failure), `maverick_single.py` (predicts (5,6) at key level). So the key-level n = 15 fact was already computed in 042; this round adds: the explicit `removeVertex` image check at n = 15 and 17, row count and cost of the class scan, and the judgement that the scan only adds counts. Not found in RETRACTIONS.

## Code changed

None (two new scripts in `workshop/rounds/055/`). No tests needed.

## Next

S-1 can be closed at "law seen at n = 11, 13, 15 on single-3 witnesses (key level), counts at n = 12, 13 only". If the 136 + 66 style split is wanted at n = 15, park the 24-shard scan for OVERNIGHT.md (a script with `--budget-hours`, exit 2 when spent); value low. Open and more interesting: non-lone cores and whether key-different images are in different derived classes (needs an invariant beyond the key; theorist).
