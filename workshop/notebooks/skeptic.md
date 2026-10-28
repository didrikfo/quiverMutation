# Skeptic's notebook (after round 094)

## Believe now
- r094 (T5): skeptic_convention2.py (R1 script from r091) ran unmodified. n = 5..7 acyclic monomial parents with zero relations up to arrow length 5 (458 parents with length >= 3): F child = End(T_R)^op (1522 steps), R child = End(T_L)^op (1529), other side/transpose 0 matches; repo gate = tilting of that side in all ~16 400 pairs per side, R gate included. Dimension + arrow data only. Extends E-252; not a new law. (Submission workshop/rounds/094/skeptic.md.)
- r091 (T5): convention for monomial acyclic parents. Fosse's convention (P_i = paths from i): F = mutateAtVertex = his RIGHT = AI mu^-; R = mutateLeftAtVertex = left. Oppermann: F = LEFT, R = RIGHT (his cone formula, l. 900; Hom convention not stated). Child = End(T)^op, no transpose. Cyclic children: arrows only (E-193). Fosse sentence NOT quoted (no tex in sources/; 2112.08129 reachable by curl from arxiv e-print).
- r087 (T4): own silting engine on LNAs (skeptic_homotopy.py) confirms H1 at w5 L10/L11, w4 L10/L11; sign map was data-picked, now explained by r091, not re-checked.
- r085 (T10): 3060000 vs 3304000 non-Cartan certificate only as quipu-theorem name; 2400230 ~ 3304000 certified so 6+8 no-meet is a false negative of the J = 0 search. r082: E-229 key-kept step is real J != 0. r074 70 failing edges = 60 keys = 12 classes. r069 E-147 13,16,9 cap artefacts. r046 Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r094: ran skeptic_convention2.py sweeps (45 s, 105 s, 122 s). r091: skeptic_convention.py. r087 homotopy. r085 cert. r082 kept/join. r079/r077 f010join. r074 failkids. r071 verify. r069 c1_13. r066 side/left. r062/058 wrongalg. r053/050 tilt/replay. r046 collect/inv/power/iso.
- Lessons: `pkill -f NAME` kills own shell; (-1)**k with negative int is float; isolated vertices vanish from path_algebra quivers; End algebra of a mutation sits on a RELABELLED line in LNA engines; calibrate sign on other cells then freeze; `| tail` buffers stdout; print label is not a test; "length 3-5" is ambiguous (vertex vs arrow length) -- state which.

## Not done
- Relations (multiplication) of the child not compared, only dims; non-monomial / parallel-arrow / cyclic parents; arrow-length-5 relations seen on 1 parent only; gate counts not split by relation length.
- Fosse sentence unquoted. T4 widths 6..11; T10 join 3304000 to 2400230.

## Next
1. Extend skeptic_convention.py to relation-level comparison (compose-check arrow products), commutativity relations. 2. Re-run r087 calibration against r091's convention. 3. Run engine on the 70 key-kept J != 0 failing edges (n = 7).

## Habits
- Check vacuity by definition; control as deep as the case and both sides; key-equal alternative control; caps and depth; independent confirming route; invariant with no demonstrated power is not evidence; say what is shared with the audited thing; check what a script really calls; truth by theory not tested code; read the convention from the paper's definition, then test; grep EXPERIMENTS before saying "not recorded".
