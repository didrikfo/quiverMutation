# Review of workshop/rounds/081/experimentalist.md

referee: skeptic · round: 081
verdict: minor revision

## Reproduction

Re-ran both commands (`experimentalist_jkey.py`, `timeout 600`). 3304000: 2m49s, TOTAL gate 10564, J 69, all 69 `J:moved`, 0 kept, ball sizes 1/11/64/282/1088/4046. 3060000: 3m20s, TOTAL gate 16573, J 99, 98 moved, 1 kept (level 4), ball sizes 1/13/83/386/1584/6354, same KEPT line. Every number in the table matches. (Runs were quicker than the author's 300 s / 560 s; not an issue.) Did not run `join97.py plan`, which is sizing text only.

## True?

The counts hold. Two points the write-up does not pin down.
- The kept step prints "child tiltingPlus-J ok", and the TOTAL line shows tp 0, key 0 failures. So under the E-212 `moves` the step is gate-admissible, and "J != 0" is relative to `perI`, not an independent Hom test. The author says "not Hom-tested", which is honest. The consequence is that the headline "contains a key-kept J != 0 step" is a statement about the `perI` flag; whether it is a real J != 0 step (as opposed to a flag artefact) is open.
- "1 of 99" cannot support a rate. The "6e-5 rate consistent with E-192's 0 in 52 190" comparison is one event against one event of the opposite sign; drop it or say it is uninformative.

No counterexample to the stated counts found.

## New?

Partly. Greps of `research/EXPERIMENTS.md` for E-219, E-204, E-226 and "J != 0 step" give:
- E-219 already has the 0.6% figure and flags it as needing reconciliation; E-226 repeats it.
- E-204 already records 70 key-kept `tiltingPlus`-false edges at n = 7/8, and E-151/E-192/E-213 cover key-kept J != 0 steps. So key-kept J != 0 steps existing is not new in general; only the split of the n = 9 0.6% into 98+1 and 69+0, and the one n = 9 instance, are new.
- No R- or F- entry that I found covers the key-moved/key-kept split at n = 9.

## Evidenced?

Yes for the counts: stated by level, reproducible in 3-4 minutes, and matching E-219's totals. Not evidenced: the "7 + 7 sized at about 2 h" rates (self-reported, not rerun, and the author's own run was under load with 2.5-2.7 expansions/s versus the quoted 5-9). The title's "2 h wall" is therefore a loose estimate and should say "estimated". The 7 + 7 was not run, which the title does say.

## Scope

n = 9, one pair, E-212 `moves`, balls to depth 5. Title is faithful. The reworded E-204 proposal ("no key-kept J != 0 step up to parent depth 3 (n = 9 F-010)") is right for 3060000 (first candidate at level 3 is moved; first kept at level 4) but E-204's own balls are n = 7/8 plus a 6 + 6 n = 9 statement; the chair should check which E-204 sentence is actually about n = 9 before rewording. For 3304000 the correct wording is "none through parent depth 4", not 3.

## Required for acceptance

1. State that "J != 0" here means `perI` true under the gate, and that the kept step passes `tiltingPlus`-J; one sentence in Claim.
2. Remove or qualify the "6e-5, consistent with E-192" comparison (one event).
3. Mark the 2 h / 0.5-1 h figures as estimates from a loaded run, and give the unloaded rate range once.
4. Give the E-204 sentence being corrected by line or quote, and give 3304000 its own wording (depth 4, not 3).
5. [next round] Hom-test the kept step (already proposed); this decides whether the finding matters.
