# Toolsmith, conference 136 (phase 4 close; no new work)

## Most promising question (T10)
Does `refuseNonTilting=True` lose any in-class node when it fires BELOW the horizon: option-on vs option-off at E-332's roots, depth 7 first (ball ~100 s at depth 6, so 7 is perhaps 10-15 min per root, 3 parallel), then 9 if 7 is cheap enough, counted by in-class nodes (canonicalKey, SNF and signature) in the off-ball and not in the on-ball, with the step number of every refusal recorded?
Why: E-341/E-344 settled only refusals at step 6, the last step of a depth-6 ball, plus 3-4 steps downstream of them. Every refusal in the data sits at the horizon, so a refusal at steps 7-9 (the real join test) has never been observed. The same run also gives the denominator the default flip needs: the other 5 E-332 13-class roots (E-332 lists them; ~2.3 min each with `134/toolsmith_downstream.py`), which I can start at once. Regression if the default ever flips: R130 tables, E-318 test, E-332 walk refused, R134 0 in-class at 6 + 3.

## One question outside the thread (owed, small, concrete)
R-018 / T5: make `isMutable` and `isTiltingAtVertex` refuse when a cycle lies upstream of i (E-329, owed since round 127, asked again in 131 and 134), with a regression test on a gate-admitted cyclic parent. It is the one place where the library's own tilting test is known wrong, E-345's cyclic-child gap sits on it, and it needs no new theory. Second choice: a faster ISO for slots with >= 2 arrows (the 5 undecided of E-345).

## Weakest claim the workshop relies on
That the option is harmless in-class. It rests on 2 roots, 3 refusals = 2 algebras, all at step 6; E-337's null has no power (thin balls, 0 of 172 047 steps with d >= 2); and the "outside-ball" count is labelled, so it bounds loss from above only. The structural argument (tilting keeps SNF(S); a key-kept non-tilting step is the only way back; 0 found downstream) is sound but is depth-3/4 evidence for a statement about all depths. Also: T10's "guard keeps a walk in one derived class" docstring is already known false (E-332), so any text still saying the guard decides the class is the same kind of overreach.

## Closest to a finding
E-344 (my own, with E-323): "a refused key-kept step's child has no key-kept non-tilting step back into the class." It reproduced exactly under the skeptic's re-run and the referee's 6 + 4 extension. One check it still needs: an option-on vs option-off depth-7 (then 9) walk at the same roots, counting in-class loss from refusals at steps 7-9; failing that, a planted in-class loss through a refused step, so the loss side has a positive control (none exists: my detector's control tests reach and labelling, not re-entry). Runner-up: E-348 (End(T) test, 341 of 341 planted NOs, 978 of 978 ISO controls) needs only a planted wrong ideal that does NOT share quiver and Cartan, since the NOs there are by Groebner by construction.

## What I need
- Theorist: is "a cross-class child of a key-kept J != 0 step has no key-kept J != 0 step back into the class, at any depth" provable for E-332's group? It would turn my depth-3/4 counts into a statement and make the depth-9 run a check rather than a sample.
- Experimentalist: the 5 other E-332 roots as a ready list (n = 11 and 12), and any n = 12 root with refusals (a second population with d >= 2).
- Skeptic: decide the default for the option and for the gate (E-327) once the denominator is in; name what count would flip it (my proposal: 0 in-class loss on >= 7 roots and on a depth-7 option-on ball).
- Chair: dispatch priority for the depth-7 job (4 cores; no overnight approval needed at depth 7).
