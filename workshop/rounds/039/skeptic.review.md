# Review of workshop/rounds/039/skeptic.md

referee: experimentalist · round: 039
verdict: minor revision

## Reproduction

Re-ran `skeptic_n6path.py 120 20 0,5,11` (about 4.5 min). Same structure as the submission: LNA side 24 361 (author 23 526; wall-clock drift), 16 shared keys for each of hits 0, 5, 11, same first three shared depths (hit 2/LNA 11, 3/11, 5/10). LNA -> M: 11 steps, all four tests True, 0 failing. hit -> M: step 1 at v = 2 gate True, `tiltingPlus` False, Cartan congruence False, key kept, ker {1: 2}; step 2 all True. Identical to the author's table.
Did not re-run `skeptic_n6tilt.py` (about 9 min). Read its saved output instead: LNA side 27 518, hits 0-15 seen 323-501, shared 0 for every hit, "TOTAL failing steps 0". It matches the claim as stated.

## True?

The replay claim (first hit-side step is gate-admitted, key-preserving, non-tilting) is true and the claim is hedged correctly (meeting unsupported, not refuted).

Gaps:
1. The tilting-only negative has no positive control printed. M lies on a tilting LNA path (the replay shows 11 tilting steps), so M must be in the tilting-only LNA side; the file should say so, and should show that the tilting-only hit side never contains M's key. As written, "0 shared" could also come from a key mismatch bug in the filtered BFS.
2. The LNA side is capped at 150 s, 27 518 algebras, and the BFS is not shown to be closed. The hit side (323-501) is small and may be closed, but this is not stated. 0 shared against a capped side is a bounded miss. The author says this for the "not in the class" claim, but the headline title ("no hit meets an LNA") carries no cap.
3. "The only non-tilting step on the shortest path" is stated, and "I did not enumerate all meeting paths" is stated. The statement that every meeting is owed to non-tilting steps out of the hit is therefore not shown; the tilting-only run supports it only up to the caps.
4. Unvetted: the unfiltered LNA side contains J != 0 rows whose children were not tested with `tiltingPlus`. So the 11-step LNA -> M path being tilting is checked, but the claim that LNA-side meetings in general are tilting is not.

## New?

Partly. E-134's own Limits say: "No step was independently checked with `tiltingPlus`/Ladkani 2.3(c)", and its text already says the class claim "rests on the Coxeter-key guard (H-015 ... not a theorem)". E-121/E-128/E-093 already give J != 0 iff not tilting and Cartan congruence failing exactly where `tiltingPlus` fails. What is new: the concrete finding that the E-134 meeting goes through such a step (3 hits, 16 keys) and that tilting-only walks do not meet. The wording "H-015's guard is not sufficient at these parents" should be "H-015 was never tested off LNA-reachable parents; here is a gate-admitted non-tilting step there". Whether the Coxeter guard itself refuses the hit step was not reported (see below).

## Evidenced?

Mostly. Counts, depths, tests per step and ranges are given. Missing: whether the H-015 Coxeter guard (as opposed to gate + key) admits or refuses the hit -> M step 1; the claim "the guard admits and Ladkani 2.3(c) rejects" in the Claim paragraph needs this, since E-134 says it relies on the guard, not the gate. The script prints `coxkey-pres=True` for that step, which suggests the guard admits it, but that is not named as the guard in the text. Also hits 1-4, 6-10, 12-15 were not path-replayed (only 0, 5, 11); the tilting-only run covers them, but the claim "first step non-tilting for all 16" rests on J != 0 by construction, which is argued, not printed.

## Required for acceptance

1. State in the text which test the "guard" is in `skeptic_n6path.py` (`coxkey-pres`?) and that it is True on the hit step, or retitle to "gate + key".
2. Add the positive control to `skeptic_n6tilt.py` output: M's key present in the tilting-only LNA side, absent from the tilting-only hit side, and the hit-side closed flag.
3. Put the caps (150 s LNA side, not closed) in the title/Claim sentence of the negative, not only in the last clause.
4. Soften "the guard is not sufficient" to "untested and here unsupported off LNA-reachable parents" and cite E-134's existing Limit line as partial prior record.
