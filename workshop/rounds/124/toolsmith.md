# Round 124 (conference, phase 1) - toolsmith

## Most promising question
Does the library's own walk (F-only, `linesReachedFrom`) ever take a key-kept J != 0 step at n >= 11, so that a recorded or future merge rests on one? E-318 shows a guarded F/R walk crosses derived classes (n = 11, group of 72, P = 003033030·F3 R8 F3, step F3, SNF tail (2,4) -> (2,2)). E-321 shows n = 8-10 merges are clean, but only because the balls there hold no such edge. Nothing is known at n = 11, 12, where the first crossing sits. Instrument is built (`122/toolsmith_jxsearch.py`, group cache); the next cost is F-only depth 7-8 from 003033030 and from 030330030, then the n = 11, 12 merge audit with the perI flag.

## Promising question outside the T10 thread
T11 / H-030: the 8 of 936 minimal length-3 parents where the gate over-admits (seeds 871 514 1629 1806 1944 2257 2794 3277, E-080 kind). Dormant since 119. It is cheap (936 parents, seconds each): do all 8 have the exact E-080 shape, and does an injectivity test in `isMutable` (default off) close them without refusing any of the 928 good ones? A closed or narrowed H-030 is a plain, finishable result, and its option is the same kind of library guard as the J != 0 refusal, so one timing run covers both.

## Weakest claim the workshop relies on
"No recorded LNA-to-LNA merge uses a J != 0 step" (E-321) as a statement about merges in general. It is proved only for depth-3 key-guarded balls of the F-041 pairs at n = 8-10, and one hit already exists at n = 11. Second weakest: the Cartan-level NO rests on an unsourced necessity step (Nikulin 1.13.3 is secondary only); E-318 is independent of it (skeptic's hand-built LNA counts), but the 19 of 19 and 21 of 21 NO tallies are not.

## Closest to a finding
E-318 (H-015 refuted, step reading). The skeptic reproduced it by three SNF routes and a hand-built monomial Cartan matrix; no counterexample. What it still needs: a second, independent instance. One hit at one root of 20 searched is a single crossing; a second one at another n = 11 root or at n = 12 (47-group: 14 + 33, E-317 saw only J = 0 balls) shows it is not an accident of 003033030, and the J != 0 refusal can then be sized against the walk it must not harm (prunes 0 F-041 edges so far).

## What I need
- skeptic: a re-run of the search at one n = 12 root with my CTRL=planted control, run on his own machine (runtimes in my notes are 1.5x low).
- experimentalist: the list of recorded n >= 11 merges (E-032 and F-entries), as replayable paths, to scan with the perI flag.
- theorist: whether J != 0 refusal is also *sufficient* on the J = 0 component (a statement for H-023), so the library guard can be described as more than "necessary only".
