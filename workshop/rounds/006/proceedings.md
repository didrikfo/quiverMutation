# Round 006 -- proceedings (ordinary)

Step 0.5: round 005 left no open questions; nothing to decide.
Worked: toolsmith (T3/T8), theorist (T2/T4), scholar (T5). Referees: skeptic (scholar, toolsmith), experimentalist (theorist).

## toolsmith -- orbit+mirror vs key over the catalogue; the 20300 pairs
Claim: orbit-plus-mirror refines the key in all 139 placed cores at n = 10, 12..16 (key finer or incomparable: 0; equal in 129-132); the 20300 pairs of `4056`, `46`, `3355`, `3445` at n = 16 are one orbit and its mirror.
Referee (skeptic): **accept**; all numbers reproduced, `mirrors` loose reading is exactly orbit-of-mirror. Decision: **accept**. Promoted: **E-066**; H-021 status line; GLOSSARY (orbit-plus-mirror class).

## theorist -- mechanism for `k(33x) = 2x`
Claim: drift `33x@o -> 33(x-1)@(o+1)` conserves `c = x + o`; the seed `333` identifies `c` with `n - c`; hence `s = n - 2x`, `d = x - 3`. Extension to `44x` fails; `34x` only at x = 5.
Referee (experimentalist): **minor revision**; 33x reproduces at n = 14 and n = 17 (x = 7); mechanism new vs E-063; "conservation law" overstates; test covers only 33y rows; 34x `k = x + 3` law has no power beyond x = 5.
Decision: **accept with the referee's wording** (I applied it myself rather than call a revision: retitled as a conserved label along a drift chain, range stated as n = 14..16, x <= 8 plus the referee's n = 17, 34x demoted to an observation). Not done: the full-orbit size comparison (the referee's point 3) and an output file for the n = 18 rule-table claim; both go on the next theorist's list. Promoted: **E-067**, H-021 status line, GLOSSARY (drift chain).

## scholar -- literature on E-032 step 7
Claim: the step-7 rejection is correct; explicit kernel element `c = [8,6,4] + [8,10,4]`; AI 2.32(b) = Ladkani 2.3(c) = `tiltingPlus` as one map; CHZ Cor 3.6 path-wise form would pass step 7.
Referee (skeptic): **minor revision**; witness reproduces; mostly known (E-057, E-032, F-038); "n = 10 first size", the one-map identity and the CHZ caveat are unsupported.
Decision: **note**, promoted only what survived: **E-068** (the witness and side check; the rest marked untested/UNVERIFIED) and an UNVERIFIED flag in `literature/2509.12983`. I did not add the "monomial only" caveat as fact.

## Questions for the steering committee
1. The proposed agenda from round 005 is still unapproved; ordinary rounds keep working from it. Recommend: approve as is (rank 1, orbit-vs-key, is now answered by E-066; T3 is done up to the parity cores).
2. Overnight: `toolsmith_orbitclass.py` takes 10+ min per slice at n = 16 and does not yet accept `--budget-hours`; a `--max-word 5` comparison at n = 14 needs `--plan` first. Recommend: do not propose an overnight run yet; ask the toolsmith to size it.

## Decisions taken for the steering committee
(none this round)
