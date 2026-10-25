# Skeptic's notebook (after round 087)

## Believe now
- r087 (T4): own silting-mutation engine (complexes of projectives over the LNA, Hom_K mod 32003, minimal approximation, cone, cancellation; skeptic_homotopy.py) confirms H1 at w5 L10 (1008/1008), plus w5 L11 3168, w4 L10 528, w4 L11 1716: T tilting, End(T) = predicted LNA (dims + arrow composites). Controls reject 100%: wrong result, wrong sequence, Coxeter-equal alternative LNA (this last has power where E-207 Coxeter has none). Sign map calibrated on w3/w4 L<=8 then frozen; mirror convention also passes. Takes seconds, widths 6..8 and 9..11 should fit too. H-020 status unchanged. Rule table itself was read, not re-derived.
- r085 (T10): non-Cartan certificate for 3060000 vs 3304000 exists only as the quipu-theorem name: 3060000 -> P^(1,0,1)_(1,4); 3304000 not theorem-nameable but moves to 3004000/3024000 -> P^(1,1,2)_(1,2). Power table n = 9: 1002 key-equal same-name pairs, 32 different-name, 7149 unnamed, 0 conflicts. 2400230 ~ 3304000 certified, so the 6+8 no-meet is a false negative of the J = 0 key-kept search. No certificate for non-LNA algebras.
- r082 (T10): E-229's key-kept step is a real J != 0 step; child joins by J = 0 at total 7. r079 E-219 repaired. r074 70 failing edges = 60 keys = 12 classes. r069 E-147 13,16,9 cap artefacts. r066 mutateAtVertex = Oppermann left mutation. r046 Cartan invariants have no power. Unit of evidence = orbit.

## Tried
- r087: skeptic_homotopy.py (modes run, ctl_after, ctl_seq, ctl_cox, calib). r085: skeptic_cert.py. r082: skeptic_kept.py, skeptic_join.py. r079/r077 f010join. r074 failkids*. r071 verify (E-207). r069 c1_13*. r066 side. r062/r058 wrongalg*. r053/r050 tilt/replay/failsteps. r046 collect/inv/power/iso/rebuild.
- Lessons: `pkill -f NAME` kills your own shell; (-1)**k with negative int is float in numpy; the End algebra of a mutation sits on a RELABELLED line, match along the arrow line not by label; calibrate the left/right sign on other cells, then freeze; orbitOf is forward closure; stdout through `| tail` is buffered and truncates; a print label is not a test.

## Not done
- T4: widths 6..11 with the new engine, anchored rules, step-by-step comparison with mutateAtVertex, rule table re-derivation. T10: join 3304000 to 2400230; premise test needing an equivalent pair crossing a J != 0 step; hand-check E-211 table (the new engine can do this: it handles non-LNA intermediate algebras). T5: op/left-right convention of Thm 3.1 (2112.08129) vs repo; the engine's convention (OP0 SIGN-1 or mirror) is a data-picked match only.

## Next
1. Run skeptic_homotopy.py on w6..w8 cells and floating w9..11 (shards). 2. Use the engine on the 70 key-kept J != 0 failing edges (n = 7): does it also call them non-tilting, independently of Hom in repo? 3. E-211 cyclic-D(i) counterexample with the engine.

## Habits
- Check vacuity by definition; run a control as deep as the case AND on both sides/classes; a key-equal alternative control; caps and depth; labelled key vs iso class; independent confirming route; invariant with no demonstrated power is not evidence; say which parts are shared with the audited thing; check what a script really calls; truth by theory not tested code; say how scarce non-iso variants are; grep EXPERIMENTS for a root name before saying "not recorded".
