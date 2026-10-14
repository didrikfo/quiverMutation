# Maverick, round 064 (conference)

**Most promising question for the next few rounds:** at n = 11, do the four
periodic groups (535/64/24/15, all signature (2,0,9)) each contain two or more
derived classes, and do those classes merge under vertex deletion to the two
classes at n = 10 (D_10 and W)? Why: signature separates nothing there
(E-190), so only a class count inside each group, and a deletion check that
ties it to n = 10, can say whether Candidate C is really open. Level: idea.
Cost is one script at n = 11 (about 40 s for the periodic set).

**Outside the thread that got most of the last 4 rounds (T10, the J = 0
premise):** the Euler-form signature of C + C^T (H-017, T6, parked and
dormant). Question: does the Smith form of C + C^T separate the F-010 pair at
n = 9 from the E-154 pair, and does it agree with the Phi-group split at n = 10
(D_10 vs W)? An invariant that is not Cartan-determined could test the premise
from a different side, which the End(T) and Hom checks cannot do (E-181,
E-184). Level: idea. One line of output is enough to decide.

**Weakest claim the workshop relies on:** "25 of 25 failing n = 7 children are
in the derived class" (T10). It rests on the J = 0 premise: generation is
proved per loopless step (E-168), but G2 holds only as a span statement
(E-183), and no join test has returned NO across a J != 0 step (E-189, E-190).
Every power control so far is circular or weak. Second-weakest: E-147's
class-1 count, reproduced in 12 of 16 steps only (E-151).

**What I need:**
- theorist: is the uniform H^T the intended left/right convention for G2, and
  does the cancellation step (Oppermann Prop., l. 1027) need the generator
  count of step 7, not just the span?
- toolsmith: the derived class count inside each of the four n = 11 periodic
  groups (535/64/24/15), with the n = 10 deletion of one representative per
  class; and Smith form of C + C^T for the F-010 pair and E-154 (one line).
- skeptic: the sign convention for [S] = -Phi against the repo's `tiltingPlus`
  (Phi^a = (-1)^(a+b)); a check that the n = 10 D_10 / W split is not a
  convention artefact.
- scholar: whether 0911.5137 or 1310.1557 already gives the class count for
  Phi-periodic LNAs (so the n = 11 question may have a known answer).

**Notes:** I ran nothing this round and edited no other file. My earlier
claim that Candidate C is moot at n = 10 is kept; no change to any H-, F-, or
R- status is proposed.
