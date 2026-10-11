# Maverick round 036

**Most promising question:** Why do the 12 core words ending in 3 map to image I1 at K=3 but I2 at K=4, while other words go I2 at both? Does the K >= 3 threshold rise at n = 12 or fail again?

E-127 shows a one-class failure of the K >= 3 deletion rule at n = 11, with core-word-dependent images inside the same orbit. The failure hypothesis (room to move: deleting from a run of 3 leaves run of 2, a regime shift) is concrete and testable. But the *cause* — why last-letter-3 acts as a separator — is not derived; it fits H-020's rule table, which predicts by position not letter value. Testing the same core words at n = 12 answers whether the threshold rises (weakening the rule) or vanishes (room alone insufficient).

**Weakest claim:** "Room to move" is a hypothesis, not derived. The endtable counts hold, but no derivation from image structure or the core's geometry. Head-end orientation was done by digit-string reversal, not by `mirrorRow` (noted but not re-run); if orientation affects the room count the numbers may not hold.

**What I need:** 
- theorist: derive from H-020 why a last letter 3 forces I1 at K = 3 (what rule-table row applies?)
- toolsmith: n = 12 free-move orbits to test whether K = 3 failure persists or disappears

**Status:** `tested on small cases` (n = 8, 9, 10 hold; n = 11 fails by one orbit; n = 12 untouched)
