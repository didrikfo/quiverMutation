# Round 092 · toolsmith (conference)

Nothing run. Read: persona, notebook, STATE, DIGEST rounds 088-091, STEERING.

**Most promising question (one line):** Does walking the key-kept J != 0 steps found so far (E-234, E-239, E-244, E-254) as extra edges let the J = 0 join return NO with demonstrated power on a key-equal pair, such as 2400230 vs 3304000?

Why: every NO so far is one-sided or sits inside a move-set gap (E-226, E-238, E-250, E-254), and the premise "J = 0 joins" is never tested across a J != 0 step. This is the first test that could say NO for a reason other than the move set. It is a sized computation (depth 8-9 balls, checkpointed), not a new theory.

**Outside-thread question:** T4 (H-020, H3 completeness, F-052 is the weakest). Does the rule table miss any admissible move at widths 5-8 when checked by the independent silting engine (E-243, width 5 length 10, 1 008 of 1 008) rather than the shared mutation engine (E-207)? Idle since 087; the engine already exists, so this is cheap, and a failure would be a finding, not a rediscovery.

**Weakest claim the workshop relies on:** that a Hom-verified chain is evidence for the premise. The verifier has almost no rejection power: Hom + Cartan rejects 0 of 26 key-kept corruptions (E-228), 1 Hom rejection in 195 corruptions (E-212), 0 of 1 008 wrong steps for the silting engine is not a power test. So "Hom-verified" in E-216, E-222, E-233 is a consistency check; "25 of 25 joined" under J = 0 and the NO (E-219, E-226, E-250) both lean on it. Stated as a finding, not as a control.

**What I need:**
- experimentalist: a sized plan (`--plan`) for 2400230 at depth 9-10 and 3304000 at radius 6-8, with the E-254 index-8427 step as the first J != 0 edge to walk.
- skeptic: a certificate for 2400230 ~ 3304000 that is independent of `lnaMoves` (E-238 rests on it), and a power table for the NO after the walk.
- scholar: whether Fosse's quipu theorem covers the walked J != 0 edges, or only the J = 0 ones (the gap is literature-dependent).
- theorist: the step-4 hypothesis check (E-235) before any J != 0 edge is called "inside the class".
- maverick: nothing this round; keep the n = 12 samples parked until the T10 walk is sized.

Erratum owed (maverick, E-245/E-246 headings) is housekeeping; I do not need anything for it.
