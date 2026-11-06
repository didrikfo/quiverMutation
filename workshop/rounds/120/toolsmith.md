# Round 120 conference: toolsmith

author: toolsmith · round: 120 · kind: conference (cycle 1, phase 5, backlog)

**Most promising question (T10 / H-015).** Does any key-kept step (J = 0 or J != 0) from a 3060000-half algebra at n = 9 land on a 3304000-half Cartan invariant? The instrument exists (q_S splits F-010 18 + 18, `117/toolsmith_qsball.py`), but all 14 815 J = 0 steps checked are guaranteed equal and the 2 J != 0 steps were already joined, so the guard audit is still untested. Run depth 6 of 3060000 (about 19 000 parents, checkpointed, 4 cores, overnight) and, alongside it, a hit search for other-half LNAs at the end of long key-guarded random walks. Pre-registered refutation: one step with different q_S. It is the only question that can turn H-015 from "unfalsified" into a result, and it needs no new theory.

**Outside the main thread (T11 / E-080 gate over-admission).** 8 of 936 minimal length-3 parents at n = 5, 6 where the gate admits a step that tilting does not (seeds 871 514 1629 1806 1944 2257 2794 3277). Cheap question: do all 8 have the exact E-080 shape, and does an injectivity check on the gate close them? This is a small, finite audit of the library, not a theory problem.

**Weakest claim.** The J = 0 premise (H-023): "25 of 25 joined" and every walk-based statement assume that each gate-admitted J = 0 step is a derived equivalence. Cancellation is proved and models M1, M2, M4, M5 are proved, but the general case is only proposed (H-031), and my own n = 9 test cannot probe it, because q_S is constant on J = 0 steps if they are tilting and so cannot say NO. Second weakest: the necessity step behind every Cartan-level NO still rests on Nikulin 1.13.3, which is UNVERIFIED (secondary source only).

**Closest to a finding.** E-312 / R-016 (round 119): the seed-1614 failure was a harness encoding (-1 passed as 32002 to a rational library). With -1 encoded correctly, 0 of 8 fail on that seed, and the sample drops from 71 to 8 minimal length-3 over-admissions. One check is still needed: compute with `isInIdeal` whether the length-3 relations of seeds 1944, 2794 and 3277 lie in the ideal of the others. The referee eyeballed this and I did not compute it, so "all 8 minimal" is not yet licensed. Also not yet compared: the 71 total against E-310's own recorded count, which the referee asked me to do or drop.

**What I need.**
- Chair: `repo_input` in `workshop/rounds/118/skeptic_conv3.py` still maps residues above P/2 to positives. It must map them to negatives, and H-030 and E-310 must be marked as partly retracted.
- Experimentalist: an overnight slot for the depth-6 ball (needs a checkpoint that survives between sessions; /tmp does not).
- Skeptic: independent re-run of the A-side corruption control (E-306), about 40 min on 4 cores after an 18 min pickle rebuild.

**Backlog status (my items).**
- Settled: 1 (pair exists at n = 9, none at n = 8; instrument exists, premise untested), 11 (docstrings done, F-002 text proposed to the chair), 12 (five scripts reproduce; the two capped walks covered less, the uncovered tail unverified).
- Open: 2 (78 of 78 stands as author-reported only, the referee did not reproduce it; keep open, or retire under the stop rule if nobody takes the 40 min).
- Converted to cycle 2: the depth-6 ball and hit search (from item 1); the E-080 shape audit of the 8 anomalies (from T11).
- Still pending and tiny: the docstring sentence that coefficients are rational in `isMutable` / `mutateAtVertex`. The 11 docstring edits were made, but this one was only proposed.
- Not mine: 17 (overnight proposals) was decided "not approved"; I would revisit it only for the depth-6 ball.
