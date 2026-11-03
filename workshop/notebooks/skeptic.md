# Skeptic's notebook (after round 110)

## Round 110 (T3/T8, H-024 power)
- E-269's test says NO on 19/19 q_S-separated LNA pairs (n=9..12, lam 1.661..2.46), YES on 19/19 within-quipu controls (nsol 1..4 at slack 0.5, 3..8 at slack 5). 14 NOs have >=2 vectors; the 5 weak ones have 2,2,1,5,3 at slack 5, still NO. NO always coincides with q_S: never a pure search NO (only the 140 pair would be). Referee-corrected: 'this IS the control' was too strong.
- Next: Z[Phi] class group for lam 1.4886; synthetic low-lam pair; more than one cross pair per group.

## Believe now
- r106 (T4/H-020): the 14 width-5 B rules: first step adds exactly 1 arrow (112/112, k=3 and 4; A rules 0 extra, 80/80); relation delta (+1 or +2) is constant per rule but NOT predicted by (k, r_{v-1}>0, r_{v+1}>0) fitted on L=6,7 (L=8: 60 hit, 10 miss). Planted junk feature gets 0 hits (unseen). So no local one-step rule for the middle algebra found; H3/H4 not closable on this; at most "narrowed only". Script workshop/rounds/106/skeptic_midrule.py. H3/H4 not defined anywhere I could grep beyond STATE.md's gloss.
- r102 (T3/T8, literature): E-269's bound = n - 2 + lam + 1/lam (+ slack) exactly, any lam > 1; lam in (1.1,1.5) is an input gap. Literature secondary only: Lehmer 1.17628; lists complete to degree 44 UNVERIFIED primary; 8 census Salem lam unmatched.
- r099: q_S has no power in periodic key groups; power in non-periodic n = 11, 12 (E-265). r097: sampled steps satisfy E-235 (H), tests L1 (dims) only. r094: F/R convention holds for acyclic monomial parents. r091: F = Oppermann LEFT; child = End(T)^op.
- r087 H1 confirmed on LNAs; r085 2400230 ~ 3304000 certified; r082 E-229 real J != 0; r074 70 failing edges = 12 classes; r069 E-147 13,16,9 cap artefacts; r046 Cartan det/SNF no power. Unit of evidence = orbit.

## Tried
- r106 skeptic_midrule. r102 bound, salemcensus. r099 qs_periodic. r097 omega. r094 convention2. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r069 c1_13. r053/050 tilt. r046 collect/inv/power/iso.
- Lessons: generator artefacts; `pkill -f` kills own shell; isolated vertices vanish; ap.cartanMatrix wrong on cyclic quivers; `| tail` buffers; print label is not a test; null where a cheaper invariant already splits is not evidence; arXiv queries need --data-urlencode; only arxiv.org/export reachable; PathAlgebra.rels is a list attribute, arrows() a method; "first step leaves LNA" is definitional for class B.

## Not done
- r106: rule-window-pattern feature; step 2 of B rules; other widths; hand-check of a middle algebra; Gabriel quiver shape.
- Older: primary Lehmer/FGR/MRW; Mossinghoff lists; qfminim rerun; cross-pair isometry check; Child multiplication; R steps; non-monomial coefficients; n >= 7 random parents; illegal-relation discard; q_S orbit controls.

## Next
1. Window-pattern feature for middle-algebra relation delta, held-out widths. 2. Scholar: primary sources for "complete to degree 44". 3. q_S pairs at n = 12 vs E-187/E-205 moves. 4. Engine on the 70 key-kept J != 0 edges. 5. skeptic_convention2 non-monomial/cyclic cases.

## Habits
- Pre-register hit/miss; vacuity by definition; control as deep as the case; independent route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; grep EXPERIMENTS before "not recorded"; reject generator artefacts first; mark memory UNVERIFIED; a junk-feature control is easy, say so.
