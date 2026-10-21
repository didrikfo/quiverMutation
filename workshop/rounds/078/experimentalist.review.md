# Review of workshop/rounds/078/experimentalist.md

referee: theorist · round: 078
verdict: minor revision

## Reproduction

Full `run 5 4 30000` (6 m 49 s) not re-run. Ran the same script at `run 3 3 30000`: 26 s. A levels [1,16,126,710], B levels [1,16,124,674], shared 0, MISS. These match the submission's table and `experimentalist_join10_out.txt` at levels 1-3 exactly, and the key-equality line matches too. Levels 4-5, the shard-plan timings (62.5 min/shard, 33.7 ms/row) and the 6 + 6 / 7 + 7 extrapolations were not re-run.

## True?

No error found at the depths I ran. Reading the code, the meet test (intersection of the two `seen` dicts, which include the roots and all levels) does bound the total labelled-key path length by 5 + 4 = 9, as stated. The cap (30 000) is not reached (18 776 keys), so no level is truncated.

Gaps:
- The claim hangs on the moves being the right J = 0 F/R set. `moves` is exec'd from rounds/049 with the base key reset. The "stats" line (gate 15 082, ok 15 082, J 0, key 0) only says no step was rejected. It does not say the move set is complete. E-212 checked Hom-verification only at n = 7 and records the verifier's NO-power as thin. At n = 10 no step was verified.
- There is no positive control at n = 10. A MISS means nothing until the same code returns a HIT on some pair it should join (for example a pair at known short tilting distance, or 90000000 against a ball point it generated). The author concedes this in Next, but then the title's "no meet" is a null with unknown sensitivity. The labelled `canonicalKey` loses to E-164 relabelling (E-212 itself missed a class-1 child at 6 + 5 this way), so this MISS is weaker than the same MISS at iso level.
- "Every gated step was admissible, no key-changing step" is a property of these balls. It is not evidence that the tilting graph at this pair stays inside J = 0.

## New?

Grepped `research/` for "tilting-only", "90000000", "00000030", E-190/E-210/E-212/E-218. Nothing records this join. E-218 (line 37) already gives the pair as a "usable long-path control for the n = 10 join search", so the object is anticipated and this is its first run. The E-173 key-only scan sizing (5 500 s) is confirmed by the submission's own number, as it says.

## Evidenced?

The join part is specific (radii, key counts per level, cap, one run, command) and would be believed without re-running at levels 1-3. Levels 4-5 are a single run with no log of the killed `run 5 5`. The shard plan is sizing only. The 16x deletion-key cost rests on one 400-row sample (summary only; no output file is cited). The "x4.3 per level" growth is from levels 3-5 (710, 3 376, 14 547: ratios 4.76, 4.31), while the 6 + 6 figure of ~50 CPU-min is an extrapolation from them and is labelled est. Fine, but it does not support "a join to total depth 12 costs ...".

## Scope

Title says "no meet at depth 5 + 4 (MISS, a bound)". That matches. The word "bound" is accurate only for the labelled-key sense, which the Claim (1) text already says. Narrowed wording: "no shared canonicalKey between the J = 0 F/R balls of radius 5 and 4 around the two roots, with the rounds/049 move set, unverified against Hom at n = 10." Claim (3) is a design for a 24-shard plan, nothing run; it belongs under "proposal", and the title's "(proposal)" does flag it.

## Required for acceptance

1. Give the n = 10 positive control: a pair the same `moves`/`canonicalKey` machinery must join (or show 90000000 reaches a node, then run from that node), with the depth at which the HIT appears. If it cannot be done in one command, mark it `[next round]` and keep the "MISS" labelled as uncalibrated.
2. Answer your own theorist line: count how many of the 646 LNAs with the target key lie in the A-ball (18 776 keys) and in the B-ball. This is a set intersection on keys, seconds. It tells whether the hereditary side is reached.
3. Save the output of `shardplan15.py plan 400` to a file and cite it, so the 33.7 ms/row and 62.5 min/shard numbers can be checked without re-running (about 2.5 min to regenerate).
4. State in the Claim that the move set is the rounds/049 one, not Hom-verified at n = 10, and that the labelled meet can hide a meet at iso level (this is already half there; make it part of the claim line).
