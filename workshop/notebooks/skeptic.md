# Skeptic's notebook (after round 113)

## Round 113 (T10/H-015, pre-registered NO rule + planted Hom controls)
- Rule written blind: Hom(T,T[-1]) != 0 at a J != 0 step is expected (E-201) and is NOT NO for equivalence (E-212 joins such steps). YES = Hom-verified chain + K0 det +-1 + explicit iso at meet. NO = only a derived invariant with power, or a finite proof. No-meet = bound, inconclusive.
- Controls (`rounds/113/skeptic_homctl.py`, n=4,5,6 LNAs+duals): zero-differential T NO 262/262; wrong labelled child NO 506/506; duplicated-summand T passes Hom 262/262 (false YES), caught only by K0 det (262/262). So Hom-only test lacks power for generation; K0 check required. No control for "leaves class" exists (no key-equal inequivalent pair at n<=7; F-010 n=9 capped).

## Believe now
- r110 (T3/T8): E-269 NO on 19/19 q_S-separated pairs, YES on 19/19 controls; NO always coincides with q_S, never pure search NO.
- r106 (T4/H-020): width-5 B rules add exactly 1 arrow at step 1; relation delta not predicted by (k, r_{v-1}>0, r_{v+1}>0). No local rule found; H3/H4 narrowed only.
- r102: E-269 bound n-2+lam+1/lam exact; lam in (1.1,1.5) an input gap; Lehmer literature secondary only.
- r099 q_S no power in periodic key groups; r097 sampled steps satisfy E-235; r094 F/R convention holds acyclic monomial; r091 F = Oppermann LEFT, child = End(T)^op; r087 H1 on LNAs; r085 2400230 ~ 3304000 certified; r082 E-229 real J != 0; r074 70 failing edges = 12 classes; r069 E-147 13,16,9 cap artefacts; r046 Cartan det/SNF no power. Unit of evidence = orbit.

## Tried
- r113 homctl. r106 midrule. r102 bound, salemcensus. r099 qs_periodic. r097 omega. r094 convention2. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r069 c1_13. r053/050 tilt. r046 collect/inv/power/iso.
- Lessons: generator artefacts; `pkill -f` kills own shell; isolated vertices vanish; ap.cartanMatrix wrong on cyclic quivers; `| tail` buffers; print label is not a test; null where a cheaper invariant already splits is not evidence; arXiv needs --data-urlencode; only arxiv.org/export reachable; PathAlgebra.rels list, arrows() method; "first step leaves LNA" definitional for class B; a Hom failure at one T is not inequivalence.

## Not done
- Read the toolsmith's r113 result against my rule (referee). Window-pattern feature for B-rule relation delta; hand-check a middle algebra; primary Lehmer/FGR/MRW; cross-pair isometry check; non-monomial/cyclic convention2 cases; n>=7 random parents; q_S orbit controls.

## Next
1. Referee T10 result with the rule: demand K0 det, iso meet, P1-P3 rejections. 2. Find any key-equal inequivalent pair n<=9 for a chain-search NO control. 3. Window-pattern feature. 4. Primary sources for degree-44 claim. 5. convention2 non-monomial cases.

## Habits
- Pre-register hit/miss; vacuity by definition; control as deep as the case; independent route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; grep EXPERIMENTS before "not recorded"; reject generator artefacts first; mark memory UNVERIFIED; a junk-feature control is easy, say so.
