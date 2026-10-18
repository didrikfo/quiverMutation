# Eight of the nine n = 7 class-2 key-kept J != 0 steps have parent and child joined by Hom-verified tilting chains with isomorphism at the meets; the NO cases tested are mostly legality, isomorphism and forced key rejections (1 Hom rejection among 195 corruptions), no key-equal inequivalent control exists at n = 7, and the n = 9 F-010 NO control (4+4) found no meet

author: toolsmith · round: 073 · kind: result
thread: T10 / T5 · bears on: H-015, E-157, E-161, E-204
scope: n = 7, key class 2 only, 8 of 9 steps (record 6 undecided), 0 key-equal inequivalent controls at n = 7, the 9 E-151 failing steps (rebuilt at 20 000 expansions), not the 36 class-2 edges of E-204; class 1 not run; generation of K^b(proj) by T assumed (E-168); chains are the printed E-157 paths, search caps quoted below; n = 9 NO control at depth 4 + 4 only

## Response to referee

1. Title: reworded as the referee proposes (NO cases are legality, iso, forced key rejections; 1 Hom rejection in 195; no n = 7 control; F-010 4+4 no meet). The Claim already said the Hom test cannot be shown to see inequivalence at n = 7; it now also says the "NO" power rests on the table, not on a key-equal control. Done.
2. `parallel-at-meet` at both meets, from `toolsmith_certificate_out.txt` (parent side / child side), all False: rec 0: False/False; rec 1: False/False; rec 2: False/False; rec 3: False/False; rec 4: False/False; rec 5: False/False; rec 7: False/False; rec 8: False/False. Seed-to-seed links (`seedlinkiso`: #4~#0, #7~#0, #9~#11) print `parallel False`. Printed key: grep `parallel-at-meet` in the output file gives 16 lines, all False. Done, no re-run needed (the lines were already in the file; the referee's one-line grep only showed some).
3. Prior record: E-157 already joins parent and child separately to an LNA of their class; the new content is the seed-to-seed link plus isomorphism (not labelled key) at meets. Added below. Done.
4. n = 9 F-010 5+5 NO control: not done this sitting (it exceeds the 10-minute cap; YES control alone took 507 s). Proposal for `OVERNIGHT.md` (I may not edit it): `timeout 6h .venv/bin/python -W ignore workshop/rounds/073/toolsmith_f010_overnight.py /tmp/ts73/c2.pkl --budget-hours 5 --depth 5` (exit 0 finished, 2 budget spent with message, 1 error; script written and its exit 2 path tested with a 18 s budget). Verdict is a bound unless a meet is found.
Not re-run: referee reproduced the chain mode; keymoving/f010 numbers remain from my output files.

## Claim

For eight of the nine class-2 failing steps (records 0-5, 7, 8; record 6 has no child path, E-157's miss), the parent A and the child B of the J != 0 step are derived equivalent: A and B are each joined to a seed LNA by a chain of F/R tilting edges, every edge accepted by an independent Hom(T,T[m]) test, and the two LNAs are equal or isomorphic (records 2, 3, 4), isomorphic as bound quivers (records 0, 1, 5: LNA 04000 and its relation dual are the same algebra), or joined by a short Hom-verified chain (record 7: 5 edges, seeds 9 ~ 11; record 8: 2 edges, seeds 7 ~ 0). So at these 8 the J != 0 step lands in the parent's derived class, in spite of Hom(T,T[-1]) != 0 (dim 1 at all 8) for the failing T itself; this is not a claim about the other 28 class-2 edges of E-204 (distinct children not counted), nor class 1, nor the Hom test being able to see inequivalence at n = 7: it cannot, because no key-equal inequivalent pair exists there to test it on. It would be refuted by an edge of a chain failing the Hom test or an iso test at a meet, none of which happened.

## Evidence

Verifier (`toolsmith_certificate.py`): an edge passes iff Hom(T,T[-1]) = Hom(T,T[1]) = 0 (`skeptic_tilt.py`, no `tiltingPlus` / `perI` / gate / key) and dim Hom(T_i,T_j) = Cartan of the child; a meeting passes iff an explicit networkx isomorphism of the structure graph (vertices, arrows, relation terms) exists. This replaces `canonicalKey`, which is LABELLED: seed pairs (10,11), (8,9), (2,3), (4,7), (5,6), (1,12) of class 2 are isomorphic yet have different keys (the E-164 point, now seen on the seeds).

Table (record = index among the 9 failing parents; edges = parent side + child side):
| rec | step | parent->seed | child->seed | seeds | edges |
|---|---|---|---|---|---|
| 0, 1 | F5, F2 | #10 | #11 | iso (04000 is self-mirror) | 18 |
| 2 | F5 | #10 | #10 | same | 15 |
| 3, 4 | F6, F3 | #11 | #11 | same | 18, 15 |
| 5 | F3 | #11 | #10 | iso | 18 |
| 7 | F1 | #9 | #11 | Hom-verified chain, 5 edges (depth 5+5, 60 s) | 17 + 5 |
| 8 | F5 | #7 | #0 | #7 iso #4; #4 to #0 by R6 R6 (2 edges) | 17 + 2 |
| 6 | F? | - | MISS (E-157) | undecided | - |
Record 4 in full: parent R3 R2 || LNA #11 F6 F5 F5 F7 F7; child F6 F1 F1 || #11 R2 R3 R4 F3 R4; all 15 edges and both meets pass.

Power (what can say NO):
- Single-edge corruptions of record 4's chain (other vertex or kind at each of the 13 positions): 195 of 195 rejected. Causes: 125 no legal step, 70 non-isomorphic meet, 1 Hom failure, so this mostly tests legality and the iso, not Hom. Swapping the LNA end for each of the other 13 seeds: 13 of 13 rejected (the moves are labelled, so this is weak for isomorphic seeds).
- Hom test on ground truth that uses no tilting code: Coxeter key moved => child certainly not derived equivalent. Complete class-2 ball to parent depth 8, sharded 4 ways (`toolsmith_keymoving.py`, each shard 590 s budget, not stopped): key-moving children rejected 14/14, 15/15, 10/10, 13/13 per shard (shards repeat the first 7 levels, so do not add); key-kept `tiltingPlus`-true children accepted 37 956, 37 951, 37 991, 37 940 per shard, 0 disagreements; key-kept FAIL rejected 4, 1, 1 in the last-level shards 1-3 (all Hom(T,T[-1]) = 1). Rejecting a key-moving step only says T is not tilting, which is forced; the information is the 0 disagreements over ~150k verdicts.
- Inequivalent key-equal control, the only kind that tests the chain finder: n = 9 F-010 pair 3060000 vs 3304000 (Ladkani Cor 3.15, E-057), same search with meets up to isomorphism (WL hash then iso): NO control depth 4 + 4, balls 1584 and 1088, no meet, 116 s. Paired YES control: end of a random 8-step walk from 3060000, meet found at 5 edges, Hom-verified (55 s). Depth 5 + 5: YES control meets at 7 (507 s); the NO control hit the 10 minute cap and was not finished. CAP: absence of a meet at 4 + 4 is a bound, not a proof; E-204 had J = 0 7 + 7 with labelled meets and also no join.
- Cartan-level invariants have no NO power on that pair: det, SNF(C), SNF(C+C^T), SNF(C-C^T), signature of C+C^T, Coxeter polynomial are all equal for 3060000 and 3304000 (`toolsmith_f010inv.py`; answers the STATE request on Smith form of C + C^T). So any certificate that only reads the Cartan matrix cannot say NO here.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/049/toolsmith_collect.py 7 2 20000 /tmp/ts73/c2.pkl 100   # 311 s (needed by all modes)
.venv/bin/python -W ignore workshop/rounds/073/toolsmith_certificate.py /tmp/ts73/c2.pkl 2 4 chain      # <1 s; idx 0..8 (not 6)
... certificate.py /tmp/ts73/c2.pkl 2 4 corrupt   # 17 s;  edges 2 s;  isoseeds 5 s
... certificate.py /tmp/ts73/c2.pkl 2 4 seedlinkiso 9 11 5   # 60 s (also 4 0 5, 7 0 5: 33 s)
... certificate.py /tmp/ts73/c2.pkl 2 4 f010 4 1             # 116 s;  `f010 5 2`: YES control 507 s, NO control not finished at the cap
timeout 10m .venv/bin/python -W ignore workshop/rounds/073/toolsmith_keymoving.py 7 2 8 590 K/4   # K = 0..3, ~400 s each
.venv/bin/python -W ignore workshop/rounds/073/toolsmith_f010inv.py   # seconds
```
Outputs kept: `toolsmith_certificate_out.txt`, `toolsmith_keymoving_out.txt`, `toolsmith_f010_d4.txt`.

## Prior record

E-157/E-160/E-163 already join each parent and each child separately to an LNA of their class under the J = 0 premise (so the new content is the seed-to-seed link plus isomorphism meets); they joined the children and parents to LNAs under the J = 0 premise and E-161 Hom-tested those paths; E-174/E-177/E-186 showed End(T) comparisons cannot see J; E-204 left "distinct children of the 70, and the join of those not in E-157" open. New: the parent-to-child certificate assembled per step with the iso at the meets (not labelled keys), the seed-to-seed links (the E-157 "LNA of their class" for parent and child can be different seeds; this was not stated), and the finding that the labelled `canonicalKey` separates isomorphic seeds. Not new: the Hom edge verdicts themselves (E-161, E-184) and F-010 no-join (E-057, E-204).

## Code changed

New scripts only, in `workshop/rounds/073/`: `toolsmith_certificate.py`, `toolsmith_keymoving.py` (copy of 070's cyccount with the Hom test), `toolsmith_f010inv.py`, `toolsmith_f010_overnight.py` (budgeted wrapper). No library file touched, no tests run (none apply). Caveat of the iso test: parallel arrows are counted but arrow identity inside relation terms is collapsed (necessary condition there); no parallel-arrow bundle occurred at any meet of these 8.

## Next

- The 36 class-2 edges of E-204 are not decided: the same verifier needs their parents (rebuild balls to depth 8 with the collector; its walk caps at 20 000 expansions) and a path search; class 1 likewise (8 parallel-arrow parents). Experimentalist/toolsmith.
- Record 6 stays open: no child path at 6 + 5.
- A key-equal inequivalent control at n = 7 does not exist; the live one is n = 9 F-010 at depth 5 + 5 (about 20 min with iso meets): proposal for `OVERNIGHT.md`: `timeout 6h .venv/bin/python -W ignore workshop/rounds/073/toolsmith_f010_overnight.py /tmp/ts73/c2.pkl --budget-hours 5 --depth 5` (exits 2 when the budget is spent).
- Skeptic: the labelled `canonicalKey` misses isomorphic meets; `meetingPoints` relabelling-aware (E-164) would also shorten the E-157 paths.
