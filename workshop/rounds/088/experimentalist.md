# Experimentalist, round 088 (conference)

Read: notebook, STATE, DIGEST entries 084-087, STEERING suggested questions. Nothing run.

**Most promising question (next few rounds):** does the n = 9 F-010 pair 3060000 vs 3304000 get a J = 0 meet at 7 + 7, the total where YES controls join (E-233, E-234)? It is the only NO the premise has (E-219, E-226), and both runs so far stopped at 5 + 5 and 6 + 6, below the YES range. A meet at 7 + 7 would make the NO test comparable to the YES controls; no meet there would be the first real evidence against the J = 0 join test. Sized at about 2 h (E-229, estimate), not approved.

**One question outside the thread that took most of rounds 084-087 (dormant line):** S-1, vertex deletion (E-218, E-236). Last worked in 083. Item 2 asks whether the change under deletion follows a rule. The fit |orbit B_n| = n 2^(n-5) (E-218) rests on n = 6..12 and has not been tested on a new length. Prediction for n = 13, checked by orbit counting before any class work, is a cheap test of a rule rather than a fit.

**Weakest claim the workshop relies on:** "25 of 25 failing n = 7 children are in the derived class". It rests on the J = 0 premise, which no test has shown to say NO across a J != 0 step. The one NO (3060000 vs 3304000, E-219, E-226) is one-sided: 2400230 does not join 3304000 at 6 + 8, and the 3304000 move-set gap (E-238) is open. The YES controls join at totals 7-8 (E-233, E-234), and the NO has not been run there.

**What I need:**
- toolsmith: `--plan` for the 7 + 7 F-010 meet and a checkpointed depth-6 build of the 3304000 ball (E-244 is at depth 6 with 11 355 parents).
- skeptic: which move set joins 2400230 to 3304000 (E-238 says the move set, not a certificate, is missing); an independent inequivalence certificate for the F-010 pair with power.
- theorist: the step-4 hypotheses check of E-235 (whether D = 0 on the spectators holds after gr_G) and the loops of degree >= 1 at i. Without it the cancellation claim stays conditional.
- maverick: a Cartan-incongruent key-equal control at n = 13 if the T3/T8 slot is taken (n = 12 has none, E-242), sized with `--plan`; its ~40 min is not approved.
- human: approve or refuse the 7 + 7 run (about 2 h) in `OVERNIGHT.md`; or say whether the S-1 n = 13 orbit check is in scope.

**Fallback if none of the above is approved:** my own item 1 from the notebook, the 12 remaining m = 2 E-230 cells at k = 4 plus the (10,4)(5,2) control (~4 h core time, no timeout). It tests H-010 further (T7), not the premise, so it ranks below the NO question.
