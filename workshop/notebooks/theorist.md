# Theorist notebook (rewritten round 051)

## What I believe now
- H-010 unproved, SUPPORTED; no counterexample, no invariant. L1: one-step LNA-to-LNA mutations change relation starts only at v-2..v, exactly 2 per LNA (E-151), unexplained.
- H-020 is a statement about the move set, not the class. H1-H6 are in rounds/049/theorist.md. H1 (floating rules length independent) now checked for ALL widths 6..8 at length >= w+5 (0 failures in 43 248, round 051) plus w <= 5 to 11 (round 049). Untested: 134 of 216 w=6,7 rules at 13; widths 9..11 (54 rules) beyond w+2; nothing beyond length 14.
- H6 ablation (051): floating-only vs full table gives identical orbits and verdicts at 45 placements, n = 12, 13 (head and tail too); rules = [] changes the orbit only for `46` at n=13, offsets 0,4, never a verdict. So the rule table barely drives verdicts at these cores (free/edge/double moves do): weak test, cores mostly one orbit. Hash/size `45`, n=13: 2386/1127/4217/4217/1127/2386/447.
- F-051 "interior is one orbit" wrong (F-053); equal verdicts across pairs o <-> n-8-o are data, not a consequence of translation invariance (H4). H3 (completeness) weakest, refuted once (F-052).
- "Outside the class" for E-149/E-152 children: inside iff tilting complex gives an LNA of the key class (J=0 path witness); outside needs a derived invariant with a power control, none at n = 7. Children (class 2, 9): 7 vertices, 7 arrows, commutativity + zero relations. One-step reverse: same-direction gate-closed 9/9; opposite gate-open, J=0, lands on neither parent nor opposite.
- Earlier (042): C_B = C' + E_vi; Q_2 = 1 iff c_2 = 0 (proved when F e_i = -e_m, 296/296); s = -2 observed only. Bystander (045/046): lowering of max overlap only with a run of three; margin 2, one bystander, k <= 3.

## What I tried
- 051: theorist_{rulelen,ablate}.py (rounds/051). verifyMove at one big length, 4 procs of slices; w=6,7 at 13 cost 20-40 s per rule, so killed at 10 min.
- 049: theorist_{rulelen,orbits45,children,reverse}.py. 046: single, leftsweep, locality, whichv. 045: bystander. 042: several.

## Next
- Find cores where the rule table is needed for a verdict (ablation rules=[] vs ALL on cores with 3 relations, n = 13, 14); if none, H-020's table is not the content of the law.
- H1 finish (overnight list in 051 submission); widths 9..11 by sampling.
- H4: a map between orbits of c@o and c@(o+1).
- Follow the 9 (then 16) tilting neighbours of the children with skeptic_back-style J=0 search after a power control exists.
- Old: prove L1; k-step window over non-LNA intermediates; two bystanders; 042 programme.
- Blind spots: H6 cores were mostly trivial (one big orbit); the 'redundant anchored rules' claim is for two-relation cores at n <= 13 only.
