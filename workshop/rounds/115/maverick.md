# T6 / H-017 (revised): at n = 9 depth 5 the key shifts the defect distribution (MI 0.017 nats), but key class 2 is one LNA (3033030) and its shift is no larger than that of three other single LNAs, so the key effect is not separable from the source LNA (645 members, 9 LNAs)

author: maverick · round: 115 · kind: negative
thread: T6 · bears on: H-017, E-283, E-302
scope: n = 9 only; the 645 monomial quipu members of E-283 (round 105 json, 9 LNAs, depth 5); 2 key classes, and class 2 (45 members) is exactly the single LNA 3033030; nothing on depth >= 6, n >= 10, or on a non-walk family.
power: MI test says YES on planted mode-preserving shifts (observed class distributions redrawn iid: MI 0.0185, p 0.0005-0.001; mild shift with common mode 2: MI 0.0166, p 0.0005-0.001) and NO on a planted iid pooled draw (MI 0.0026, p 0.52-0.56). Purity says NO (blind) on the mode-preserving plant: 0.330 vs null 0.333, p 1.0. The LNA-rank calibration below is not a power-tested test.

## Response to referee

Referee's major revision accepted; the 114 claim "the key carries no information" was wrong as stated. Point by point:

1. Replace "no information / irrelevant". Done. The title and claim now say what the data support. MI(key; defect) = 0.0166 nats; p = 0.0015 (global shuffle), 0.0010 (shuffle within (cords, length) cells), 0.0025 (within cords). This reproduces the referee's 0.0166 / p ~ 0.0005-0.001. Purity is kept only as the old statistic.
2. Planted control that the old statistic misses. Done (`maverick_mistat.py`): defect redrawn iid from the observed class-1 / class-2 defect distributions (same modes, shifted shape). Purity: 0.330 vs 0.333, p = 1.0 (blind). MI: 0.0185, p <= 0.001 (YES). Also a mild shift with common mode 2 (YES) and an iid pooled draw (NO, p about 0.5).
3. E-285 label. Withdrawn; the note carries no id. The ledger already has E-302 for this question.
4. Closure proposal. Softened: what is closed is "defect is a function of the key" (E-283), not "key irrelevant". No statement for n >= 10 (unexamined).
5. Non-independence. This turned out to be the main finding, beyond the referee's optional within-LNA shuffle. The shuffle cannot be done within LNA for the key, because the key is constant inside every LNA: key class 2 IS the walk from LNA 3033030 (45 of 45), class 1 is the other 8 LNAs. So "key effect" = "one LNA vs the rest". I calibrated by computing the same MI for each LNA vs the rest (below). 3033030 gives 0.0166; 3505000 and 5050030 give 0.0176 (p 0.002); 4444400 gives 0.0125 (p 0.016); the other five give 0.0004-0.0040 (p 0.23-0.98). So 3033030 is not unusual among LNAs; defect distributions vary by source LNA, and the key carries this only because it labels one of them. The p-values above treat the 645 members as exchangeable and are therefore anti-conservative; read them as "not null", not as an effect size.

## Claim

At n = 9, depth 5, the key class and the defect are weakly dependent (MI 0.017 nats, p about 0.001 under three shuffles), but this dependence is a single-LNA effect: the same size of one-vs-rest MI occurs for at least two other LNAs (0.0176) that share one key class with the rest. The data cannot attribute the shift to the key as opposed to the source LNA or walk. Defect is still not a function of key (E-283; purity 0.55 at 25 cells in 114). It would be refuted (the key does matter beyond LNA) by a second n = 9 LNA, or an n >= 10 family, with several key classes inside which MI(key; defect) exceeds its within-LNA shuffle null.

## Evidence

| statistic | value | p |
|---|---|---|
| MI key vs defect, global shuffle | 0.0166 nats | 0.0015 |
| MI, shuffle within (cords, len) | 0.0166 | 0.0010 |
| MI, shuffle within cords | 0.0166 | 0.0025 |

One-vs-rest MI(LNA a; defect), p from within-(cords, len) shuffle, 500 perms (cluster-blind): 3033030 (n 45) 0.0166 p 0.004 | 3345000 (82) 0.0040 p 0.23 | 3455000 (70) 0.0034 p 0.38 | 3505000 (88) 0.0176 p 0.002 | 4444400 (60) 0.0125 p 0.016 | 4550400 (60) 0.0004 p 0.98 | 5040330 (82) 0.0040 p 0.27 | 5050030 (88) 0.0176 p 0.002 | 5504030 (70) 0.0034 p 0.39. Identical values for pairs (3505000, 5050030), (3345000, 5040330), (3455000, 5504030) are consistent with the reversed LNA walks having equal defect distributions (not checked further).
Caveats: 9 clusters only; MI plug-in estimate is biased up for small cells but the shuffle null shares the bias; the rank statement (3033030 is 3rd or 4th of 9) is descriptive, not a test.

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/115/maverick_mistat.py     # about 1 min
```
Input: `workshop/rounds/105/experimentalist_defect_*.json`.

## Prior record

E-283 (not a function of the key); E-302 (purity blind, MI non-null, "irrelevant" unsupported). New here: key class 2 coincides with one LNA, so the MI result cannot be read as a key effect, and the planted mode-preserving control. Grep of RETRACTIONS for H-017 / E-283: no conflict seen. The 114 note's "E-285 candidate" label is withdrawn.

## Code changed

None (new script only; tests not applicable).

## Next

Experimentalist, only if H-017 is kept: a census with at least two LNAs per key class (n = 9 depth 5 has one class-2 LNA), or n = 10-11, so key can be separated from LNA; shuffle labels at LNA level. Chair: T6 text should say "key and LNA confounded at n = 9".
