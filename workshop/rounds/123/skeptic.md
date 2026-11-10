# Recorded merges only (n = 8-10): the 19 n = 9 F-010 joins each contain exactly one J != 0 step, the audited one, with q_S constant at every node; the library's LNA-to-LNA merges contain none, and that is forced: the depth-3 key-guarded balls of all F-041 pairs at n = 8, 9, 10 contain no key-kept J != 0 edge; E-318's walk is flagged

author: skeptic · round: 123 · kind: result (largely confirmation; see Prior record)
thread: T10 · bears on: H-015, H-032, H-023, E-318, E-316, F-037, F-041, F-054
scope: recorded merges with a printed or regenerable path only. (a) n = 9 F-010 joins of the 19 audited key-kept J != 0 steps (E-233, E-234 x10, E-239, E-244 x6, E-254) plus E-233's child-to-parent join: 14 replayed from printed paths, 5 (E-234 terminal-only) regenerated with the original `skeptic_join.py` (totals reproduce), E-244 8427 through E-255's isomorphism to 8063's child. (b) the library's own LNA-to-LNA merges: F-037 (5 printed paths of 19 hits), F-054 (meet and mixed certificate), the C(2,4,5) atlas meet (one path beyond E-169's 5), F-041 `meetingPoints` depth 3 at n = 8, 9, 10 (regenerated: no paths were stored beyond n = 8; every meeting pair has exactly one meeting key), and the whole depth-3 key-guarded ball of every F-041 pair at n = 8, 9, 10 (all 10, 65, 304 pairs, both ends, LNA and opposite). Not covered: n = 7 "25 of 25" paths (E-157/E-160/E-163; needs the E-151 rebuild, > 10 min; INV has no power at n = 7); E-032's n = 10/11 overnight merges (no paths recorded); F-037's other 14 hits (no recorded path); balls of other starts (E-318's n = 11 start has a J != 0 edge); multi-vertex steps.
power: E-318's 5-step walk (003033030 -F3 R8 F3 F3 F6-> 030330030, NO = not a derived equivalence, known) is flagged: step 4 J != 0, `tiltingPlus` false, SNF(S) tail (1,2,4) -> (1,2,2); the other 4 steps J = 0 with INV equal. Planted cross merges at n = 9 (E-233's child path against 3304000's root path, E-239's against 3060000's): both fail at the meet `iso`, not at a step (0 INV changes at their steps); one INV comparison, q_S(3060000) != q_S(3304000), shared by both. They test the meet check, not the step detector. For (b) the detector had **nothing to find**: the F-041 balls contain 0 key-kept J != 0 edges (n = 8, 9, 10, all pairs), so the 0 in (b) is forced and carries no information about the walk's choice of path. SNF(S) and signature are EQUAL for the two F-010 halves, so at n = 9 only q_S can say NO; "SNF constant" there is vacuous.

## Response to referee

All six required items done; the claim, title, *scope* and *power* lines are changed where the answers change them. New script `skeptic_ballJ.py`; `skeptic_f041.py` gains the key check; `skeptic_mergereplay.py` prints the planted controls honestly.

1. **E-169.** Accepted: E-169 (round 054) already replayed the C(2,4,5) / group-A merge (5 paths, 35 edges, J = 0, key kept), and its Consequences supersede E-164's "unreplayed". The C(2,4,5) row is now one more path (F-054's 4 + 3 atlas meet); Prior record corrected.
2. **Ball-level count.** Done for **all** pairs at n = 8, 9 and 10 (the full n = 10, not a sample). `skeptic_ballJ.py` walks the same depth-3 key-guarded gate walk as `quiversReachedFrom(., 3, alsoDual=True)` on both ends, with a visitor, and classifies every walk edge by `perI` and `tiltingPlus`. Per-pair distinct edges (parent `quiverKey`, vertex, side): n = 8 **2 396**, n = 9 **18 518**, n = 10 pairs 0-15 **3 065** (all three equal the referee's `ballJ` counts), all n = 10 **117 642**. **0 J != 0, 0 `tiltingPlus` false.** So the referee is right: the 0 in (b) is a property of the balls, not of the merges. (b), the call answer and *power* are reworded as Scope suggests. The "informative count is on (b), 0" sentence is gone.
3. **Planted controls.** Accepted. There is one INV comparison, INV(3060000) vs INV(3304000); the old line printed it twice. Both plants fail at the meet `iso`, not at a step (0 INV changes at their steps). They test the meet check, not the step detector, and they are not independent. The script now prints this.
4. **F-037.** Now "5 printed paths of 19 hits"; the other 14 have no recorded path and were not replayed.
5. **Tied meetings.** There are none. Every pair that meets has exactly one shared key: n = 8 10 of 10, n = 9 52 of 52, n = 10 206 of 206, from `skeptic_ballJ.py` over all pairs. Both replayed halves end at that key, `quiverKey(end) == key` on both sides, at n = 8, 9, 10 (all) in `ballJ` and at n = 8 in `skeptic_f041.py`. The exclusion is deleted.
6. **(a) and H-023.** Accepted: "gives it no test" becomes "a weak necessary-condition test (INV varies in the key group: the halves differ in q_S); none failed".

## Claim

Every recorded merge was replayed step by step with the E-317 classifier (gate, illegal relation, J != 0 = `perI` non-empty, `tiltingPlus`, Coxeter key) and INV = (SNF(S), signature, q_S histogram when det S != 0), S = C + C^T, at every node.

(a) The n = 9 F-010 joins are not F-only and not LNA-to-LNA merges: each is root -> P (mixed F/R) -> audited step -> child -> (mixed) meet ~iso~ meet <- (16 of 19 root halves F-only) root. **19 of 19 contain exactly one J != 0 step, the audited one, by design** (18 audited steps with their own join, plus E-233's child-to-parent join; the 19th audited step, E-244 8427, has no join of its own and rests on E-255's isomorphism, which holds); the other 239 steps are J = 0, `tiltingPlus` true, gate-admitted, key kept. q_S, SNF(S) and signature are constant at every node of every path (0 changes). All meets are isomorphisms (19 of 19).

(b) The library's LNA-to-LNA merges (key-guarded gate walks, no J filter) contain **0 J != 0 steps and 0 INV changes**: F-037 35 of 35 edges (5 printed paths of 19 hits), F-054 3 + 3, C(2,4,5) 7, F-041 n = 8 60 + 60 inverse edges (= E-156), n = 9 52 merges / 312 + 312 edges (13 pairs no meeting = F-041's 13), n = 10 206 merges / 1 236 + 1 236 edges (98 pairs no meeting = F-041's 98). **This 0 is forced:** the depth-3 key-guarded balls of the F-041 pairs at n = 8, 9, 10 (all pairs) contain no key-kept J != 0 edge (per-pair distinct walk edges 2 396, 18 518, 117 642; 0 J != 0, 0 `tiltingPlus` false), so no F-041 witness can use one; F-037 (5 of 19 hits printed), F-054 and C(2,4,5) (with E-164, E-169) likewise use none. The F-041 ends are derived equivalent by theorem (two-arrow deletion), so there INV could never say NO; (b) counts J != 0 edges, it is not a crossing test.

So the answer to the call is: 19 recorded merges contain a J != 0 step, all of type (a) and all with q_S (and SNF) constant; no recorded LNA-to-LNA merge does, and for F-041 none could, since the balls they are drawn from have no J != 0 edge. **Not claimed:** that the library walk cannot reach an E-318-type crossing from other starts or depths (E-318 itself is one); anything about n = 7, 11, 12; that the J = 0 steps are derived equivalences (H-023 is assumed; the INV constancy along them is a weak necessary-condition test of it, since INV varies in the key group (the halves differ in q_S); none failed). **Refuted by:** a recorded merge path with a J != 0 step whose two ends differ in q_S or SNF(S), or a replay error in the listed paths (the meets must be isomorphic, and they are).

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
| F-037 34504030 -> 50505000, 5 printed paths of 19 hits, 7 F moves each | 5 | 35 | 0 | 0 | 0 | yes (end iso 50505000) |
| F-054 30330300 F7 || 30330400 F8 F7; mixed F7 R9 R10 | 2 | 6 | 0 | 0 | 0 | yes |
| C(2,4,5) 05040330 F4 F3 F7 F2 || 60504030 R10 R7 R2 (one path beyond E-169's 5) | 1 | 7 | 0 | 0 | 0 | yes |
| F-041 n = 8 (= E-156) | 10 | 60 + 60 | 0 | 0 | 0 | 10 of 10 |
| F-041 n = 9 (13 no meeting) | 52 | 312 + 312 | 0 | 0 | 0 | 52 of 52 |
| F-041 n = 10 (98 no meeting) | 206 | 1 236 + 1 236 | 0 | 0 | 0 | 206 of 206 |

Ball-level count (`skeptic_ballJ.py N`): every walk edge of the depth-3 key-guarded gate walk from both ends of each F-041 pair, LNA and opposite. Distinct = (parent `quiverKey`, vertex, side), per pair, summed over pairs.

| n | pairs | meeting (exactly 1 key) | halves end at key | walk edges | distinct edges | J != 0 | `tiltingPlus` false |
|---|---|---|---|---|---|---|---|
| 8 | 10 | 10 (10) | 10 | 3 120 | 2 396 | 0 | 0 |
| 9 | 65 | 52 (52) | 52 | 24 852 | 18 518 | 0 | 0 |
| 10 | 304 (shards 0-15, 15-65, 65-115, 115-165, 165-215, 215-260, 260-304) | 206 (206) | 206 | 162 902 | 117 642 | 0 | 0 |

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py control     # 1 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py records     # 6 s (14 joins, 8427 iso, planted)
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py extra       # 3 s (5 regenerated E-234 joins)
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py f037        # 2 s (F-037, F-054, C(2,4,5))
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py 8                   # 16 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py 9                   # 172 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_f041.py 10 LO HI            # shards 0-15, 15-115, 115-215, 215-400: 36, 447, 528, 474 s (three in parallel on 4 cores)
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_ballJ.py 8                   # 27 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_ballJ.py 9                   # 307 s (alongside n = 10 shards)
timeout 10m .venv/bin/python -W ignore workshop/rounds/123/skeptic_ballJ.py 10 LO HI            # shards above: 67, 357, 452, 520, 431, 428, 474 s (3-4 in parallel on 4 cores)
# regeneration of the 5 E-234 joins (paths now hard-coded in `extra`), DIR any scratch directory:
.venv/bin/python -W ignore workshop/rounds/123/skeptic_mergereplay.py mkpkl DIR
timeout 10m .venv/bin/python -u -W ignore workshop/rounds/082/skeptic_join.py side R 3060000 5 560 DIR/ckR           # 115 s
for i in 1275 1433 1746 1818 1915; do .venv/bin/python -W ignore workshop/rounds/082/skeptic_join.py side C pkl:DIR/k$i.pkl:child 3 150 DIR/ckC$i; .venv/bin/python -W ignore workshop/rounds/082/skeptic_join.py meet DIR/ckC$i DIR/ckR; done   # ~5 s each
```

## Prior record

E-316 (mine, round 121) already had q_S equal at the 19 audited steps and every child in its parent's J = 0 component; the join steps were Hom-verified in their own records (`skeptic_join.py` replay, E-233..E-254); F-037's 35 edges were `tiltingPlus`-true (E-164); F-041 n = 8 is E-156. So (a) and most of (b) are **confirmation**: what is new is perI plus q_S/SNF at every node of every path, the 5 regenerated E-234 joins, F-041's n = 9/10 witness edges (no stored paths before), the ball-level J != 0 count (response item 2), and one more C(2,4,5) path (F-054's 4 + 3 atlas meet). The C(2,4,5) / group-A merge itself is **not** new: E-169 (round 054) replayed it with 5 explicit 7-step paths, 35 of 35 edges gate-admitted, J = 0 and key kept, superseding E-164's "unreplayed"; E-169 also noted that its 1 202 gate-admitted n = 10 edges had no J != 0, "not shown able to fail here", which the ball count extends to every F-041 pair at n = 8, 9, 10. Note on the call: the n = 9 joins contain a J != 0 step by construction (19 of 19 trivially), and the 0 in (b) is forced by the balls; neither count is informative about whether the walk prefers J = 0 routes.

## Code changed

None in the library. New: `workshop/rounds/123/skeptic_mergereplay.py`, `skeptic_f041.py` (now also counts halves ending at the meeting key and pairs with > 1 meeting), `skeptic_ballJ.py` (response), `skeptic_jtilt.py` (verbatim copy of E-317's `maverick_jtilt.py`, as the call asks; the replay scripts reuse its step classification). No tests run (no library file touched).

## Next

- Toolsmith: a default J != 0 refusal in `mutationSearchDepthFirst` would prune **0 edges** in the F-041 balls at n = 8, 9, 10 and in every recorded LNA-to-LNA merge; it costs only the perI check (E-151: 8-11%), and E-318 is the only known case it changes. Regression test: F-041 counts unchanged, and the E-318 walk refused.
- Toolsmith/experimentalist: E-032's n = 10/11 merges rerun with `--witness` so they can be replayed the same way; n = 7 E-160's 3 paths still lack a Hom test.
- Theorist: why the 19 audited n = 9 J != 0 steps stay in their q_S class while E-318's F3 step at n = 11 does not; compare the step shapes.
