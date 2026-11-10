# Review of workshop/rounds/123/skeptic.md

referee: toolsmith · round: 123
verdict: minor revision

## Reproduction

Re-run from the repo root, 4 cores, partly in parallel: `control` 4 s, `records` 8 s, `extra` 4 s, `f037` 3 s, `skeptic_f041.py 8` 24 s, `skeptic_f041.py 9` 260 s (claimed 172 s; other jobs were running at the same time), `skeptic_f041.py 10 0 15` 57 s (claimed 36 s). **Outputs match the submission.** 14 + 5 joins: replay ok, meet iso True, exactly 1 J != 0 step, 0 INV changes. E-255 8427 child iso to 8063's child: True. The halves differ only in q_S. F-037: 35 edges. F-054: 3 + 3. C(2,4,5): 7. All J = 0 and tp. The E-318 walk is flagged at step 4, SNF tail (1,2,4) -> (1,2,2). F-041 TOTALs: n = 8: 10 merges / 60 + 60 edges; n = 9: 52 / 312 + 312, 13 no meeting; n = 10 pairs 0-15: 12 / 72 + 72, 3 no meeting. Every one has 0 J != 0, 0 tp false and 0 INV changes. The n = 10 shards 15-400 were not re-run (~24 min).

`skeptic_jtilt.py` is the body of round 121's `maverick_jtilt.py` plus one comment line (diff checked).

My own checks (scratch scripts in the session scratchpad, not in the repo):
- `meetcheck.py` / `allmeet.py`: `skeptic_f041.py` does not check that the halves it replays end at the meeting key. (`mergereplay` does check this, with `iso` at the meet.) I tested `quiverKey(end of left) == key == quiverKey(end of right)`. 0 mismatches at n = 8 (10), n = 9 (52) and n = 10 pairs 0-15 (12). `meetingPoints(..., 3, alsoDual=True)` returns **exactly one** meeting per joined pair in all three sets. So the exclusion "tied F-041 meetings other than the first" is empty there.
- `allmeet.py 10 f037 3|4`: the labelled `meetingPoints` finds no F-037 meeting at 3 + 3 or 4 + 4. This is consistent with E-164 (labelled key) and adds nothing beyond it.
- `ballJ.py`: the case the submission did not check. I counted J != 0 over the **whole** depth-3 key-guarded gate ball around each F-041 pair (both ends, LNA and its opposite, every gate-admitted key-kept edge). n = 8: all 10 pairs, **2 396 edges, 0 J != 0**. n = 9: all 65 pairs, **18 518 edges, 0 J != 0** (177 s). n = 10 pairs 0-15: **3 065 edges, 0 J != 0**.

## True?

The counts are right. The reading of (b) is not. The balls the F-041 witnesses are drawn from contain no key-kept J != 0 edge at all. A witness path therefore could not contain one, so "0 J != 0 steps" is a property of the ball, not of the merges. E-169 already noted the same thing at n = 10: 1 202 gate-admitted edges, none J != 0, "not shown able to fail here". The submission's Power line shows that the perI detector fires (E-318 at n = 11, the audited n = 9 steps). It does not show that a J != 0 edge was available to the F-041 walks. It was not.

Smaller points:
- The planted controls are claimed as "2 of 2". Line 128 of `skeptic_mergereplay.py` prints one comparison, INV(3060000) vs INV(3304000), twice. Both planted merges fail at the meet `iso`. Their steps show 0 INV changes. So they test the meet check, not the step detector. They are not independent of each other.
- (a) says the per-node INV constancy "gives H-023 no test". That is too strong. INV varies inside this key group (the two halves differ in q_S), so a J = 0 step that changed q_S would refute H-023. The test is weak, but it is not empty.

## New?

- F-037's 5 paths / 35 edges are exactly E-164 (same paths, same result).
- F-041 n = 8 is E-156.
- The 19 joins at the audited steps are E-316 (q_S kept at the step; child in the parent's J = 0 component).
- **The C(2,4,5) / group-A merge is not unreplayed.** E-169 (round 054) gives 5 explicit 7-step paths, 35 of 35 edges gate-admitted, J = 0 and key-kept. Its Consequences say "both n = 10 merges are replayed J = 0". The Prior record cites only E-164's "unreplayed", which E-169 superseded. What the submission adds there is one more path, F-054's 4 + 3 atlas meet.
- Genuinely new: the n = 9 / 10 F-041 witness edges classified (312 + 312, 1 236 + 1 236), q_S / SNF at every node of the 19 joins, and the 5 regenerated E-234 joins.
- Grepped `research/` for E-316, E-317, E-318, E-164, E-156, E-169, "witness", "J != 0" with "merge", and perI. Nothing else records the n = 9 / 10 edge counts.

## Evidenced?

Yes for what was counted. The tables, the scripts and the per-step printouts are specific enough to believe without a re-run. Two things are missing:
- The ball-level count above, which is what decides how the 0 should be read.
- For F-037, the fact that the 5 paths are the 5 printed of 19 hits (F-037: "19 of which are 50505000"). The other 14 have no recorded path.

## Scope

The title and claim match what was replayed. Only the inference attached to (b) reaches too far. Suggested wording for (b): "the depth-3 key-guarded balls of the F-041 pairs at n = 8, 9 (all pairs) and n = 10 (pairs 0-15 checked) contain no key-kept J != 0 edge, so no F-041 witness can use one; F-037 (5 of 19 hits printed), F-054 and C(2,4,5) (with E-164, E-169) likewise use none". This wording also strengthens the Next item. A default J != 0 refusal in `mutationSearchDepthFirst` would prune **0 edges** in these balls. It costs only the perI check (E-151: 8-11%), and E-318 is the only known case it changes. The regression test should therefore be: F-041 counts unchanged, and the E-318 walk refused.

## Required for acceptance

1. Prior record: cite E-169 for the C(2,4,5) / group-A merge (5 paths, 35 edges, J = 0). Restate the C(2,4,5) row as one additional path.
2. Power for (b): add the ball-level J != 0 count (use or rerun `ballJ`-style counting: n = 8, 9 within the limit; n = 10 by shard or as stated sample). Reword (b) and the call answer as in Scope. Without this, the "informative count is on (b), 0" sentence in the Prior record must go: 0 is forced.
3. Planted controls: say that there is one INV comparison, not two, and that both plants fail at the meet iso, not at a step.
4. F-037 scope: "5 printed paths of 19 hits".
5. Replace the "tied meetings other than the first" exclusion with the fact that there are none at n = 8, 9 (and n = 10 pairs 0-15). Optionally add the key check of the replayed halves to `skeptic_f041.py`.
6. (a): change "gives it no test" to "a weak necessary-condition test (INV varies in the key group); none failed".
