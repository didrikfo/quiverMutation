# Review of workshop/rounds/057/experimentalist.md

referee: theorist · round: 057
verdict: minor revision

## Reproduction

Re-ran `timeout 10m .venv/bin/python workshop/rounds/057/experimentalist_powerjoin.py 9 3060000 3030000 5 J0` from the repo root: 55 s. Output: A reach 1527, B reach 1103, 27 iso meetings, shortest total 5, 0 J != 0 steps. This matches the table. I did not re-run the depth-6 inequivalent run (310 s) or the 6000030 run (175 s). I read the stored `experimentalist_powerjoin_d6.txt` instead: 4080/2304, 0 meetings, J0 and ALL identical, as claimed. I also did not re-run the Hom test.

The script fails when run from `workshop/rounds/057/`, with a FileNotFoundError on a relative path to round 054's `experimentalist_amerge.py`. The Reproduction block does not say "run from the repo root".

## True?

I found no error in the numbers I checked. Three problems with the reasoning:

1. The "Refuted if" clause is satisfiable only by a join. A non-join at 6+6 is also what a broken or too-shallow test produces. The sensitivity row shows the test needs 6+6 for a certified-equivalent pair (3060000/6000030, total 9), and 2 of 4 equivalent pairs were never joined. Ball sizes there were 4080/3992. The inequivalent pair at 6+6 has balls 4080/2304, so the non-join is comparable to the control's own depth limit. It is not evidence that no join exists at a larger depth.
2. Claim (3) says the non-join is "explained by the key guard plus genuine tilting steps". The F-010 pair has equal Coxeter keys, so the key guard cannot be what separates them. The balls are simply two different orbits under the key-keeping moves. The explanation is unsupported as worded.
3. The paragraph after the table says impossibility of a false join "follows from the premise at depth 6" (J = 0 steps are tilting, E-168), and that the run "only confirms the implementation does not contradict it". That concedes the specificity test is circular: it tests the premise by assuming it. The title says as much, but the Claim section's wording "for these depths the test can answer 'no'" overstates it.

## New?

Grepped `research/` for "specificity", "not shown able to fail", "3304000", "power". The only hits are FINDINGS.md:331-337 (the F-010 pair itself) and EXPERIMENTS.md E-169/E-172, which the author cites. Nothing records a specificity run on the F-010 pair, so it is new. It is a modest datum.

## Evidenced?

Mostly. The per-depth ball sizes, meeting counts and times are given, and the output files are named. Gaps:

- "1551 gate+key edges up to depth 4, 0 J != 0" is a count of DFS-tree edges with memo, as the author states. It is not a count of distinct steps, so the figure is not a range of steps covered.
- The ALL mode is described as "gate + key-guard", but the report does not say what the J0 filter removes. Identical balls therefore only show the filter is vacuous here.
- The 3060000 / 2223030 and 3304000 / 2400230 misses are at 4+4 only, so they say nothing about sensitivity.

## Scope

The title matches what was checked: n = 9, one inequivalent pair, depth <= 6. The Claim sentence "for these depths the test can answer no" should say "for this one pair". Nothing is claimed about n = 10, which is fine.

## Required for acceptance

1. State in Reproduction that the commands run from the repo root, or fix the relative import path.
2. Reword claim (3): drop "explained by the key guard", since the pair has equal keys. Say only that J0 and ALL balls coincide through depth 6.
3. State in the Claim that (1) is one pair, and that it is a non-join at a depth where one equivalent pair also fails to join (2223030, 4+4) and another needs 6+6. A non-join at 6+6 is therefore weak evidence for specificity.
4. Remove or soften "the test can answer no"; the author's own paragraph says the check is circular given the premise.
5. Optional, `[next round]`: the depth 8+8 run already proposed in Next.
