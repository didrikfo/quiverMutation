# Skeptic's notebook (after round 118)

## Round 118 (backlog 13-15)
- Item 14 (T4/H-020): first-step (arrows, relations) delta of class-B rules is a per-rule constant (14/14, 54/54, 94/94; width 7 vacuous, one match per rule; length transfer 82/82 at widths 5, 6), two values (1,1)/(1,2). 120 simple rule features: best 0.872 on width 7 = label-shuffled max 0.872, so no formula. Converted to proposed H-029; H3/H4 retired as stated for B rules.
- Item 15 (E-258): new engine `rounds/118/skeptic_conv3.py` (quotient Hom, vector composition). Parallel-arrow monomial: 0 failures (~5 800 tilting steps). Cyclic monomial: 0 failures (~1 350 steps with a repo result; ~2 450 more refused by CyclicDescendantsError, untested). Non-monomial: repo gate/F/R fail on NON-MINIMAL relation input with length-3 relations (74 anomaly lines/206 parents; same algebra minimal: 0/8); minimal + length-2 commutativity: 0/360 parents; minimal + length 3: 7/33, all gate True & tilting False (E-080 kind). Proposed H-030.
- Item 13 referee: theorist_steps.py/allJ.py reproduce E-305 exactly (13/13, 75 steps, 19 kept equal, 2 of 56 refused differ, corruptions 69/30/13). Another process shared /tmp/skeptic_c1_13_ck1.pkl at once.

## Believe now
- r113: HH/Hom-only tests lack power for generation; K0 det needed. r110 E-269 NO on 19/19 q_S pairs. r106 + r118: B-rule delta is rule-level, no local formula. r102 E-269 bound exact. r099 q_S no power in periodic groups. r094 F/R convention holds acyclic monomial; r118 extends to parallel, cyclic (arrows). r091 F = Oppermann LEFT, child = End(T)^op. r087 H1 on LNAs. r085 2400230 ~ 3304000. r082 E-229 real J != 0. r074 70 failing edges = 12 classes. r069 E-147 13,16,9 cap artefacts. Unit of evidence = orbit.

## Tried
- r118 ruledelta, rulefeat, conv3. r113 homctl. r106 midrule. r102 bound. r099 qs_periodic. r097 omega. r094 convention2. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r069 c1_13. r046 collect/inv/power/iso.
- Lessons: generator artefacts (redundant relations, a generator can make non-minimal input); best-of-120 feature fit needs a shuffle null (matched baseline); tag every counter by the property that explains the failures (minimality, relation length) before reading it; isolated vertices vanish; ap.cartanMatrix wrong on cyclic quivers; `| tail` buffers; sleep >= 100 s is blocked, use Monitor/until; shared /tmp checkpoints can be written by another session; "first step leaves LNA" definitional for class B.

## Not done
- Locate the failing step inside procedure for the redundant-presentation case; n = 7 and arrow length 4 non-monomial; child multiplication (not just dimensions); HH^2 = 1 for ab = 0 square; tree status of the 2 NO children; window-pattern is dead (null); q_S orbit controls; key-equal inequivalent pair n <= 9 for a chain-search NO control; primary Lehmer/FGR/MRW; cross-pair isometry; n >= 7 random parents.

## Next
1. Referee toolsmith's reply on the non-minimal input (regression test seed 1614). 2. H-029 test at widths 8-9 (rule table needed). 3. Key-equal inequivalent pair for a T10 NO control. 4. Primary sources for degree-44 claim.

## Habits
- Pre-register hit/miss; vacuity by definition; control as deep as the case; independent route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; grep EXPERIMENTS before "not recorded"; reject generator artefacts first; mark memory UNVERIFIED; a junk-feature control is easy, say so.
