# Toolsmith, round 132 (conference, phase 3 research)

Read: notebook; my R130 submission and its review (the only toolsmith files in 129-131; none in 129 or 131); E-329..E-340; STATE; plan. Board not read before forming this.

## Most promising question
Does `refuseNonTilting` ever fire, and what does it cost, on a population where it can fire? Run the option-on counter (one walk, wrapped `isTiltingAtVertex`, exact shortcut of R130) over the walks that contain d = max_i dim e_i A e_v >= 2 steps: E-332's depth-6 pure-F crossing, E-316's 19 mixed F/R F-010 joins, and the n = 9 merges of `merges.py` at depth 5-7. E-337 is a null that E-338 makes a theorem for thin balls, so the default-flip decision (off, R130) currently rests on a population that cannot answer it. One counted pass, no new machinery, decides whether the option loses joins that matter or only crossings.

## One question outside the main thread
T6/H-017b (dormant since 119): per-LNA defect statistic with delta* fixed first, n = 10/11. It is the only open line that is a pure compute job I can size; the q_S cache and orbit census make it cheaper than when it was parked. Alternatively E-149's 154 lost edges rechecked with `sameIdeal` (an artefact-of-code check, cheap).

## Weakest claim
"The option costs nothing" / "default stays off" (E-337 as read in R130 notes). Zero refusals in 3 634 balls and 172 047 steps, none with d >= 2: no power, as E-338 shows. Second, the J = 0 premise for walks outside thin balls; E-340 says "no LNA walk reaches a cyclic node at n <= 8, depth <= 5", which is a capped sample, and the library gate errs on cyclic parents (39 over-admissions in 18 145 rows).

## Closest to a finding
E-338: at a gate-admitted vertex with all dim e_i A e_v <= 1, J = 0 (and the simple-path gate is exact even on cyclic quivers). It is proved and checked per row (29/29, 92/92) but is a rediscovery of E-128 L1 in part. The one check it still needs: a positive control population, i.e. confirm on the E-318 planted step and the E-332 depth-6 step that d >= 2 occurs exactly at the J != 0 steps and at no J = 0 step in those walks (the converse direction, "d >= 2 plus gate => J != 0" is not claimed and may be false; count it).

## What I need
- Theorist: the exact statement of when d >= 2 at an admitted step forces J != 0 (converse of E-338), so the counter output can be read as a decision, not a sample.
- Skeptic: decide, before the run, the threshold at which `refuseNonTilting` flips default (what loss of F-010 joins is acceptable); and the non-monomial cyclic Hom engine for T5.
- Chair: if the n = 11 depth-4 pass is wanted (4.2 h CPU, ~75 min wall), a go for the background job; otherwise I do the n <= 10 and E-332/E-316 walks only.
- Toolsmith's own owed items (not blocking): make `isMutable`/`isTiltingAtVertex` refuse on a cycle upstream of i (E-340), with the 8 n = 4 + 175 long-glue regression.
