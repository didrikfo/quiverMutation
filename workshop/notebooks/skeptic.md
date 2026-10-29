# Skeptic's notebook (after round 097)

## Believe now
- r097 (T10/T5): for acyclic parents (n 4..6, random, planted omega||alpha and Anick chain) every step the real search takes satisfies E-235's (H) and dim End(T) (round-050 independent code) = repo child (5 066 + 1 848 acyclic-child steps; 271 cyclic-child leaves via truncated path counts). The E-248 degree-1 pattern is never taken: step 5 writes (alpha*, rbar) along t -> i* -> t, and `isIllegalRelation` (before `reduce`) discards it (1 519 pairs). Degree-2 analogue (Anick chain + parallel arrow) is taken (238), cyclic child, dims ok. Pre-registered: HIT = (H) fails or dims differ; MISS says nothing about generation/J. Result was a miss. Submission workshop/rounds/097/skeptic.md. Sample is random, F steps, dims only.
- r094 (T5): F/R convention holds for acyclic monomial parents, zero relations to arrow length 5 (E-258): F child = End(T_R)^op, R child = End(T_L)^op, dimension + arrow data only.
- r091: Fosse convention (P_i = paths from i): F = mutateAtVertex = his RIGHT = Oppermann LEFT. Child = End(T)^op. Fosse sentence unquoted.
- r087 own silting engine on LNAs confirms H1; r085 2400230 ~ 3304000 certified, so 6+8 no-meet is a move-set gap; r082 E-229 key-kept step is real J != 0; r074 70 failing edges = 60 keys = 12 classes; r069 E-147 13,16,9 cap artefacts; r046 Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r097: skeptic_omega.py (search with visitor), skeptic_omega_census.py (all vertices, gate/illegal/key split, trunc_cartan for cyclic children). r094 skeptic_convention2.py sweeps. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r071 verify. r069 c1_13. r062/058 wrongalg. r053/050 tilt/replay (skeptic_tilt.py is the independent Hom(T,T[m]) code). r046 collect/inv/power/iso.
- Lessons: generator artefacts (non-minimal roots gave 4 fake "lost steps" until minimality was enforced); `pkill -f` kills own shell; isolated vertices vanish; LNA End algebras sit on relabelled lines; `| tail` buffers; print label is not a test; ap.cartanMatrix is wrong on cyclic quivers (use own truncated count); "length 3-5" ambiguous (vertex vs arrow length).

## Not done
- Child multiplication not compared (dims only); R steps; non-monomial coefficients; n >= 7 random parents; whether the illegal-relation discard loses valid J = 0 edges (cyclic leaves; my 4 artefact cases hint the pre-reduce timing matters); Ext^3 pattern confirmed only by Anick-chain construction.
- Fosse sentence unquoted. T4 widths 6..11; T10 join 3304000 to 2400230.

## Next
1. Count PAT steps whose reduced child is acyclic or End-dim-equal with J = 0 (steps lost by the pre-reduce discard) on real search balls. 2. Relation-level compare (compose arrows) of child vs End(T). 3. Engine on the 70 key-kept J != 0 failing edges (n = 7).

## Habits
- Pre-register hit/miss meaning; vacuity by definition; control as deep as the case, both sides; caps and depth; independent confirming route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; read convention from definition then test; grep EXPERIMENTS before "not recorded"; reject generator artefacts before counting.
