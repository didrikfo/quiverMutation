# Skeptic's notebook (after round 091)

## Believe now
- r091 (T5): convention settled for monomial acyclic parents. Fosse's convention (P_i = paths from i, Hom(P_u,P_w) = paths w~>u): F = mutateAtVertex = his RIGHT mutation = AI mu^- (cocone of right approx by arrows out of i); R = mutateLeftAtVertex = left (cone of left approx by arrows into i). Oppermann's proof uses P_i* = Cone(P_i -> sum P_{t(alpha)}), A = arrows out of i, i.e. the dual module side, so F = his LEFT, R = his RIGHT. Child = End(T)^op labelled, no transpose; 1713/1713 F and 594/594 R steps (n = 4..7, random non-symmetric monomial parents, own K^b(proj) engine skeptic_convention.py), controls: T_L data differ from T_R in all, the other side/transpose matches 0. Gate open <=> T_R tilting in 2101/2101. Matching is dimensions + Gabriel arrow counts only, not relations. E-220's side caveat closes on this basis (names depend on module side). 249 + 53 children are cyclic: repo cartanMatrix gives -1 there (E-193), arrows still right.
- r087 (T4): own silting engine on LNAs (skeptic_homotopy.py) confirms H1 at w5 L10/L11, w4 L10/L11; controls reject 100%. Sign map was data-picked (OP0 SIGN-1 or mirror): now explained by r091, not yet re-checked.
- r085 (T10): non-Cartan certificate for 3060000 vs 3304000 only as quipu-theorem name; 2400230 ~ 3304000 certified, so the 6+8 no-meet is a false negative of the J = 0 key-kept search. r082: E-229 key-kept step is a real J != 0 step. r074 70 failing edges = 60 keys = 12 classes. r069 E-147 13,16,9 cap artefacts. r046 Cartan invariants have no power. Unit of evidence = orbit.

- r091 referee (minor revision), answered: claim narrowed to length-2 zero relations, a<b simple arrows, dimension data; 249+53 cyclic children = arrows only; A^op "Oppermann literal" run is a reformulation, not a test; gate claim F only. Oppermann side follows from his cone formula (l. 900), not a stated Hom convention. Fosse sentence NOT quoted (no tex in sources/). skeptic_convention2.py (length 3-5 relations + R gate) written but UNRUN: shell refused by classifier all sitting. Run it first next time.

## Tried
- r091: skeptic_convention.py (demo, sweep n N). r087: skeptic_homotopy.py. r085: skeptic_cert.py. r082: skeptic_kept.py, skeptic_join.py. r079/r077 f010join. r074 failkids*. r071 verify (E-207). r069 c1_13*. r066 side/left. r062/r058 wrongalg*. r053/r050 tilt/replay/failsteps. r046 collect/inv/power/iso/rebuild.
- Lessons: `pkill -f NAME` kills your own shell; (-1)**k with negative int is float (bit me again in r091); isolated vertices vanish from a path_algebra quiver (index errors); the End algebra of a mutation sits on a RELABELLED line in LNA engines; calibrate sign on other cells then freeze; stdout through `| tail` buffers; a print label is not a test; arxiv e-print IS reachable by curl (2112.08129 tex downloaded in r091, not in sources/ -- ask the chair to add it).

## Not done
- r091: relations (multiplication) not compared, only dims; non-monomial / parallel-arrow parents; cyclic parents. 2112.08129 Thm 3.1 text read from the fresh e-print (hypothesis "for all i != 0" slip noted E-208), not compared line by line to steps 4-7 of procedure.py.
- T4: widths 6..11, anchored rules. T10: join 3304000 to 2400230; E-211 table by hand (engine handles non-LNA).

## Next
1. Extend skeptic_convention.py to relation-level comparison (compose-checks of arrow products) and to commutativity relations (Hom dims > 1 would need general coefficient handling already present). 2. Re-run r087 calibration against r091's convention. 3. Run engine on the 70 key-kept J != 0 failing edges (n = 7).

## Habits
- Check vacuity by definition; control as deep as the case AND both sides; key-equal alternative control; caps and depth; independent confirming route; invariant with no demonstrated power is not evidence; say which parts are shared with the audited thing; check what a script really calls; truth by theory not tested code; read the convention from the paper's own definition, then test; grep EXPERIMENTS before saying "not recorded".
