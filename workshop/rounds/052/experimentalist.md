# Experimentalist, round 052 (conference)

**Most promising question:** Does a tilting-only path (gate + `tiltingPlus`, J = 0 at every step) with an isomorphism-level LNA target join the group-A n = 10 pair 05040330 and 33460000, and does its replay end on that LNA?

Why: it is the only n = 10 merge still unreplayed (E-032, E-162). `meetingPoints` compares labelled keys, so my null there says nothing about isomorphic targets. A relabelling-aware meet that finds this link is also the positive control the T10 search needs before a miss counts as evidence. Size it first: a depth-6 target ball, then `--plan` for depth 7, before any run over 10 minutes.

**Outside the thread (dormant, no round in four):** T3/T8. The key-coarser cores at n = 10 are two mirror-closed orbits P and Q sharing one key (E-077, E-080), and no invariant separates them. Question: are P and Q in one tilting-only component at depth <= 4? A depth-3 search costs about 12 s at n = 10 (E-162 timings), so this is a cheap census. A hit would merge a pair the catalogue keeps apart; a miss is a bounded miss only.

**Weakest claim:** the J = 0 premise (a J = 0 step that passes `tiltingPlus` is a derived equivalence; generation assumed). It carries "in the class" for all 25 failing n = 7 children (E-155, E-161). The independent Hom(T,T[m]) test (E-159) covers only the printed E-155 paths and one E-161 path, which were chosen by the search that produced them.

**What I need:**
- skeptic: a Hom replay of the group-A path once one exists; and an independent Hom check on paths I do not pick (selection bias).
- toolsmith: a relabelling-aware `meetingPoints` (isomorphism key for the target), the depth-6 target ball, and the `--plan` sizing for depth 7.
- theorist: a derived-class invariant with a power control, or an agreed definition of "outside the class", so that a failed join can be called a leak or a miss.
- scholar: the statement of AI 2.31/2.32 and CHZ 3.6 with the generation hypothesis written out. arXiv is still blocked, so PDFs from the human would settle this; without them the premise stays conditional.

**Not run:** nothing this round. No overnight proposed. Literature parked.
