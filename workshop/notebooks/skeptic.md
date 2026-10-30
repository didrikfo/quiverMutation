# Skeptic's notebook (after round 099)

## Believe now
- r099 (T3/T8): q_S of C + C^T is isomorphic for 90000000 vs 00000030 (n = 10, equivalent per E-190): MISS (pre-registered hit = non-isomorphic). The Phi^18 key group (323 rows) splits into 321 (D_10, sig (10,0)) + 2 (34504030, 50505000, sig (8,2)): exactly the E-190 classes, but signature gives the same split, so q_S has no power in periodic groups. Its power is in non-periodic groups n = 11, 12 (E-265). Submission workshop/rounds/099/skeptic.md. maverick_cartanZ imports cypari, absent here; stub it.
- r097 (T10/T5): sampled search steps (n 4..6, acyclic parents) satisfy E-235 (H) and dim End(T) = child; E-248 pattern discarded by `isIllegalRelation`. Tests L1 only (dims), not L2.
- r094: F/R convention holds for acyclic monomial parents (E-258). r091: Fosse convention (P_i = paths from i): F = mutateAtVertex = Oppermann LEFT; child = End(T)^op.
- r087 silting engine on LNAs confirms H1; r085 2400230 ~ 3304000 certified, so 6+8 no-meet is a move-set gap; r082 E-229 key-kept step real J != 0; r074 70 failing edges = 60 keys = 12 classes; r069 E-147 13,16,9 cap artefacts; r046 Cartan invariants (det, SNF) have no power. Unit of evidence = orbit.

## Tried
- r099 skeptic_qs_periodic.py (pair + whole key group partition). r097 skeptic_omega*.py. r094 convention2. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r069 c1_13. r053/050 skeptic_tilt.py (independent Hom(T,T[m])). r046 collect/inv/power/iso.
- Lessons: generator artefacts; `pkill -f` kills own shell; isolated vertices vanish; ap.cartanMatrix wrong on cyclic quivers; `| tail` buffers; print label is not a test; "length 3-5" ambiguous; a null on a group where a cheaper invariant already splits is not evidence for the new invariant.

## Not done
- Child multiplication not compared; R steps; non-monomial coefficients; n >= 7 random parents; illegal-relation discard losing J = 0 edges; q_S orbit control only on small subsets (132..645 members); 140 group Z[Phi] iso; Fosse sentence unquoted.

## Next
1. q_S vs real E-187/E-205 moves: take q_S-separated pairs at n = 12 (15 groups) and search for a move path (a hit refutes). 2. Count PAT steps lost by pre-reduce discard. 3. Engine on the 70 key-kept J != 0 edges.

## Habits
- Pre-register hit/miss; vacuity by definition; control as deep as the case, both sides; caps and depth; independent route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; grep EXPERIMENTS before "not recorded"; reject generator artefacts first.
