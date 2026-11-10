# Recorded merges only (n = 8-10): the 19 n = 9 F-010 joins each contain exactly one J != 0 step, the audited one, with q_S constant at every node; the library's LNA-to-LNA merges (F-037, F-054, C(2,4,5), F-041 first meetings at n = 8, 9, 10; 276 merges) contain none; E-318's walk is flagged

author: skeptic · round: 123 · kind: result (largely confirmation; see Prior record)
thread: T10 · bears on: H-015, H-032, H-023, E-318, E-316, F-037, F-041, F-054
scope: recorded merges with a printed or regenerable path only. (a) n = 9 F-010 joins of the 19 audited key-kept J != 0 steps (E-233, E-234 x10, E-239, E-244 x6, E-254) plus E-233's child-to-parent join: 14 replayed from printed paths, 5 (E-234 terminal-only) regenerated with the original `skeptic_join.py` (totals reproduce), E-244 8427 through E-255's isomorphism to 8063's child. (b) the library's own LNA-to-LNA merges: F-037 (5 paths), F-054 (meet and mixed certificate), the C(2,4,5) atlas meet, F-041 `meetingPoints` depth 3 at n = 8, 9, 10 (first shortest meeting per pair, regenerated: no paths were stored beyond n = 8). Not covered: n = 7 "25 of 25" paths (E-157/E-160/E-163; needs the E-151 rebuild, > 10 min; INV has no power at n = 7); E-032's n = 10/11 overnight merges (no paths recorded); tied F-041 meetings other than the first; multi-vertex steps.
power: E-318's 5-step walk (003033030 -F3 R8 F3 F3 F6-> 030330030, NO = not a derived equivalence, known) is flagged: step 4 J != 0, `tiltingPlus` false, SNF(S) tail (1,2,4) -> (1,2,2); the other 4 steps J = 0 with INV equal. Planted cross merges at n = 9 (E-233's child path against 3304000's root path, E-239's against 3060000's): meet iso False and q_S(root) != q_S(root'), 2 of 2. SNF(S) and signature are EQUAL for the two F-010 halves, so at n = 9 only q_S can say NO; "SNF constant" there is vacuous.

## Claim

Every recorded merge was replayed step by step with the E-317 classifier (gate, illegal relation, J != 0 = `perI` non-empty, `tiltingPlus`, Coxeter key) and INV = (SNF(S), signature, q_S histogram when det S != 0), S = C + C^T, at every node.

(a) The n = 9 F-010 joins are not F-only and not LNA-to-LNA merges: each is root -> P (mixed F/R) -> audited step -> child -> (mixed) meet ~iso~ meet <- (16 of 19 root halves F-only) root. **19 of 19 contain exactly one J != 0 step, the audited one, by design** (18 audited steps with their own join, plus E-233's child-to-parent join; the 19th audited step, E-244 8427, has no join of its own and rests on E-255's isomorphism, which holds); the other 239 steps are J = 0, `tiltingPlus` true, gate-admitted, key kept. q_S, SNF(S) and signature are constant at every node of every path (0 changes). All meets are isomorphisms (19 of 19).

(b) The library's LNA-to-LNA merges (key-guarded gate walks, no J filter) contain **0 J != 0 steps and 0 INV changes**: F-037 35 of 35 edges, F-054 3 + 3, C(2,4,5) 7, F-041 n = 8 60 + 60 inverse edges (= E-156), n = 9 52 merges / 312 + 312 edges (13 pairs no meeting = F-041's 13), n = 10 206 merges / 1 236 + 1 236 edges (98 pairs no meeting = F-041's 98). The F-041 ends are derived equivalent by theorem (two-arrow deletion), so there INV could never say NO; (b) counts J != 0 edges, it is not a crossing test.

So the answer to the call is: 19 recorded merges contain a J != 0 step, all of type (a) and all with q_S (and SNF) constant; no recorded LNA-to-LNA merge does. **Not claimed:** that the library walk cannot reach an E-318-type crossing (only recorded/first-shortest paths were replayed); anything about n = 7, 11, 12; that the J = 0 steps are derived equivalences (H-023 is assumed, not tested; the INV constancy along them is implied by H-023 and gives it no test). **Refuted by:** a recorded merge path with a J != 0 step whose two ends differ in q_S or SNF(S), or a replay error in the listed paths (the meets must be isomorphic, and they are).

## Evidence

(a) n = 9 F-010 joins (`records`, `extra`); every row: replay ok, meet iso True, q_S constant at every node.

| record | parent path | step | child path | other side (root path) | steps | J != 0 |
|---|---|---|---|---|---|---|
| E-233 | F1 R9 R8 R8 | R8 | R1 R3 | 3060000 F1 F7 F6 F5 F4 | 12 | 1 (step) |
| E-233 X~P | F1 R9 R8 R8 | R8 | F8 F8 | parent F8 R7 R6 R5 | 11 | 1 |
| E-234 1932 / 1970 / 2003 / 4551 | (as E-316) | R8 | F9 R1 R3 / R1 R1 R3 / R3 R1 R3 / R1 F9 R3 | 3060000 F1 F7 F6 F5 F4 | 14 each | 1 |
| E-234 4722 | R9 R8 R8 F1 F8 | R1 | R3 R8 | 3060000 F1 F7 F6 F5 F4 | 13 | 1 |
| E-234 1275 / 1433 (regen.) | F1 F4 / F1 F5 R9 R8 R8 | R8 | F8 F8 F8 | F1 F5 R9 R8 R7 / F1 F6 R9 R8 R7 | 14 each | 1 |
| E-234 1746 / 1818 / 1915 (regen.) | (as E-316) | R8 / R8 / R9 | R1 R3 | F1 F4 F7 F6 F5 / F1 F5 F7 F6 F4 / F1 F7 F6 F5 F4 | 13 each | 1 |
| E-239 2026 | F2 R8 R9 R7 R7 | R7 | R2 | 3304000 F2 F1 F6 F5 F2 | 12 | 1 |
| E-244 4606 / 8063 / 8170 / 8196 / 11123 | (as E-316) | R7 / R9 / R7 / R7 / R2 | R1 R2 R4 / R2 R4 R5 / R2 R2 R4 / R4 R2 R4 / R4 R7 R5 | F2 F1 F3 F6 F5 / F2 F1 F6 / F2 F1 F6 F5 / F2 F1 F6 F5 / F2 F1 F6 | 15/13/14/14/13 | 1 |
| E-254 | F1 F1 F2 F3 R9 R7 R9 F3 | F3 | F3 F7 F9 | 2400230 F1 R4 R5 R9 R7 R9 | 18 | 1 |
| E-255 8427 | F2 R9 R8 R7 R9 R7 | R7 | child iso 8063's child: True | -- | 7 | 1 |
| total (19 joins) | | | | | 258 (+7) | 19 (+1) |

INV of every node: SNF tail (1, 1, 12), signature (8, 1), |disc| 12; q_S histograms of the halves differ (3060000: values 0, 2/3, 1, 5/4, 5/3, 23/12; 3304000: 0, 1/4, 1/3, 1, 4/3, 19/12).

(b) library merges (`f037`, `skeptic_f041.py N`):

| merge set | merges | edges (fwd + inverse) | J != 0 | tp false | INV changes | ends INV equal |
|---|---|---|---|---|---|---|
| F-037 34504030 -> 50505000, 5 paths of 7 F moves | 5 | 35 | 0 | 0 | 0 | yes (end iso 50505000) |
| F-054 30330300 F7 || 30330400 F8 F7; mixed F7 R9 R10 | 2 | 6 | 0 | 0 | 0 | yes |
| C(2,4,5) 05040330 F4 F3 F7 F2 || 60504030 R10 R7 R2 | 1 | 7 | 0 | 0 | 0 | yes |
| F-041 n = 8 (= E-156) | 10 | 60 + 60 | 0 | 0 | 0 | 10 of 10 |
| F-041 n = 9 (13 no meeting) | 52 | 312 + 312 | 0 | 0 | 0 | 52 of 52 |
| F-041 n = 10 (98 no meeting) | 206 | 1 236 + 1 236 | 0 | 0 | 0 | 206 of 206 |

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py control     # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py records     # 6 s (14 joins, 8427 iso, planted)
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py extra       # 3 s (5 regenerated E-234 joins)
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py f037        # 2 s (F-037, F-054, C(2,4,5))
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py 8                   # 16 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py 9                   # 172 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py 10 LO HI            # shards 0-15, 15-115, 115-215, 215-400: 36, 447, 528, 474 s (three in parallel on 4 cores)
# regeneration of the 5 E-234 joins (paths now hard-coded in `extra`), DIR any scratch directory:
.venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py mkpkl DIR
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side R 3060000 5 560 DIR/ckR           # 115 s
for i in 1275 1433 1746 1818 1915; do .venv/bin/python -W ignore workshop/rounds/082/skeptic_join.py side C pkl:DIR/k$i.pkl:child 3 150 DIR/ckC$i; .venv/bin/python -W ignore workshop/rounds/082/skeptic_join.py meet DIR/ckC$i DIR/ckR; done   # ~5 s each
```

## Prior record

E-316 (mine, round 121) already had q_S equal at the 19 audited steps and every child in its parent's J = 0 component; the join steps were Hom-verified in their own records (`skeptic_join.py` replay, E-233..E-254); F-037's 35 edges were `tiltingPlus`-true (E-164); F-041 n = 8 is E-156. So (a) and most of (b) are **confirmation**: what is new is perI plus q_S/SNF at every node of every path, the 5 regenerated E-234 joins, F-041's n = 9/10 witness edges (no stored paths before), and a replay of the C(2,4,5) meet (E-164 calls the group-A merge "unreplayed"; F-054's text gives a depth-4 meet path for it). Note on the call: the n = 9 joins are not "F-only merges"; they contain a J != 0 step by construction, so "how many contain one" is 19 of 19 trivially; the informative count is on (b), 0.

## Code changed

None in the library. New: `workshop/rounds/123/skeptic_mergereplay.py`, `skeptic_f041.py`, `skeptic_jtilt.py` (verbatim copy of E-317's `maverick_jtilt.py`, as the call asks; the replay scripts reuse its step classification). No tests run (no library file touched).

## Next

- Toolsmith: with 0 J != 0 edges in every recorded LNA-to-LNA merge, a default J != 0 refusal in `mutationSearchDepthFirst` (E-151: 8-11% cost) would change no recorded merge; propose it with F-041/F-037 counts as the regression test.
- Toolsmith/experimentalist: E-032's n = 10/11 merges rerun with `--witness` so they can be replayed the same way; n = 7 E-160's 3 paths still lack a Hom test.
- Theorist: why the 19 audited n = 9 J != 0 steps stay in their q_S class while E-318's F3 step at n = 11 does not; compare the step shapes.
