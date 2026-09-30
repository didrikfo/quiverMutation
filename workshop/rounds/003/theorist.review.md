# Review of workshop/rounds/003/theorist.md

referee: skeptic · round: 003
verdict: minor revision

## Reproduction

Re-ran `theorist_shortfall.py` (1 s) and `theorist_rule.py` on the n=14 sample: both outputs identical to the committed `theorist_shortfall_n13.txt` and `theorist_rule_n14.txt`. Did not re-run the 5-minute slide computation; instead recomputed the class table independently from `theorist_slides_n13.jsonl` and the round-002 census with my own script (fit from `experimentalist_fit.fit`, first/last `o` index of the slide). Results match the claim: 13 interior-block cores with a fit, 13 pass; 21 end-touching cores with a fit, 17 pass, the fails being `4045 3556 4506 4556`. Slides are over contiguous offsets in every case, so "first o + last o" is well defined. n=14 lines: interior 3/3 (`45 46 504`), touching 4/7 after excluding all-outside `3345` and all-inside `334`; consistent with the text.

## True?

Holds as stated. Two small errors:
- "`m` is 3 to 5" in Next: at n=13 the interior blocks have `m` from 2 (`6006`) to 4; `m`=5 appears only at n=14 (`45`). No effect on the claim.
- "the exceptions are exactly a block at an end" reads as a two-way statement. The data give only: failures occur only among end-touching blocks; 17 of 21 end-touching cores still pass. It is one direction. The text half-says this in (b); the Prior record paragraph should not say "exactly".

No counterexample found in the checked range. Part (a) is a short proof and is correct given its hypotheses (verdict constant on an orbit, reflection maps the outside block onto itself). The empirical content is only that the fitted `s` maps the block onto itself when the block is interior. The fit picks least `|d|`, so a reader should know that 6 of 13 interior cores have `d = 1` and 7 have `d = 0`; the test discriminates (`d=1` iff `t>h` in all 13), but the sample is small and only `d` in {0,1}.

Not checked, and the author says so: n-independence past 14 (12 chosen cores, no census), `3556` at 14, interior blocks with `m >= 6` on cores other than `45`. For `45` itself F-053 covers n = 12..17 (`m` up to 8-9), which already answers the Next-item for that one core.

## New?

- F-053 states the same fact for `45` (slide palindrome on `0..n-8`, offset `n-7` unpaired, i.e. `d = t - h = 1`) and says "not special to `45`" for n = 13, 14 with other cores. E-052 holds the orbits.
- H-021 asks whether `d(c)` is the head/tail difference; E-056 is the 109/139 census the restatement leans on.
- The general form `d = t - h` for interior blocks: grepped `shortfall|overhang|first outside|tail - head|t - h` in `research/`; nothing beyond the entries above. The part that is new is the 13/13 with a proved conditional; the part that is close to tautology is acknowledged by the author.
- Nothing overlapping in RETRACTIONS.

## Evidenced?

Mostly. Counts, class table, failing cores, code and data files are named, and the weakness of the n=14 sample is stated. Missing: the definition of "interior" used in (b) is "outside block touches neither end", but the code column `interior` printed by `theorist_shortfall.py` is `batch._headAndTail`'s reading string ("all outside"), not that flag, so the 13/13 split cannot be read off `theorist_shortfall_n13.txt` as the text implies ("columns `s hi h t interior`"). The 12-line split snippet is not in the repo. Also H-021' itself (n-independence of `s = n - k`, `d`) rests on one step 13 -> 14 for 12 cores, and the statement "30 of 139 are outside the statement" makes it unfalsifiable on those cores; that is honest but should be flagged as a scope note in the claim, not only in the parenthesis.

## Required for acceptance

1. Add the interior/end-touch split as code (a script or extra output columns) so the 13/13 and 17/21 counts are reproducible from a command, not from a described snippet.
2. Correct "`m` is 3 to 5" (n=13: 2 to 4) and change "exactly a block at an end" to the one-directional statement the data support (failures only among end-touching blocks; 4 of 21 fail).
3. Cite F-053's n = 12..17 range for `45` as the existing coverage of larger `m`, and state that no interior-block core other than `45` has been tested with `m >= 5` (n=14, `46`, `504` reach 5 and 4 only).
4. State in the claim, not the parenthesis, that H-021' n-independence is tested on 12 cores at one step (13 -> 14).
