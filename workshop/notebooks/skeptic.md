# Skeptic's notebook (after round 126)

## Round 126 (T11 / H-030: the 8 E-312 over-admissions)
- Shape of the 8 (rounds/126/skeptic_inj.py seeds): exact E-080 (one side arrow, dim e_uAe_i = 2, kernel p - q, input relation (p-q)beta):
  514, 1629, 1944, 2257, 2794 (1629/2257/2794 on R side = dual). E-080 + one extra surviving path (dim 3): 871, 1806.
  3277: two in-arrows, p, q with equal NONZERO products on both = E-112's "nn" circuit, needs a consequence relation.
- All 8: Hom(T,T[-1]) total 1, nothing else. Injectivity (Ladkani 2.3(c) = tiltingPlus = mutation.isTiltingAtVertex) refuses all 8,
  two implementations (library copy over Q, own mod-p kernel).
- Agreement with Hom tilting: 114 372 of 114 372 (vertex, side) cases, N n = 5 (30k seeds), 6 (40k), 7 (12k), P n = 6 (3k). 0 gate-F/Hom-T.
- Not new: the test is recorded theory (E-059, E-080, E-112, E-161). New: shape split 5/2/1; nn circuit exists among random parents;
  more "other" shapes at n = 6, 7 (two-arrow; consequence-relation identifications, seeds 7926, 696).

## Believe now
- H-030 is explained: the gate checks "no PATH in the kernel", tilting needs "no ELEMENT in the kernel". Fix exists (refuseNonTilting, opt-in).
- r123: no recorded merge uses a J != 0 step except audited ones; depth-3 key-guarded balls of F-041 pairs n = 8-10 have 0 J != 0 edges (forced).
- r121: n = 9 key-kept J != 0 children land in parent's J = 0 component. r113: HH/Hom-only tests lack power. r099 q_S no power in periodic groups.
  r094 F/R convention; r091 F = Oppermann LEFT. Unit of evidence = orbit.

## Not done
- refuseNonTilting default decision (T10 request to me): measure cost on F-041 counts + E-318 refused; cost E-151 8-11 %.
- n = 7 "25 of 25": E-160's 3 paths not Hom-tested (needs E-151 c1 rebuild > 10 min).
- E-032 n = 10/11 merges need --witness reruns. H-029 widths 8-9. n = 12 47-group second instrument.
- Are nn / consequence circuits reachable from an LNA walk (E-112 says not seen)? Would need walk rows, not random parents.

## Next
1. refuseNonTilting default: regression = F-041 counts unchanged AND E-318 walk refused AND seed 514/3277 refused.
2. Referee any n = 12 T10 run: pre-registered YES, q_S/SNF per kept child.
3. Hunt: injectivity vs Hom at n = 8, coefficients != +-1 (E-112 "balanced" caveat: a circuit with gains != 1 gives no kernel).

## Habits
- Grep the TEST first: "injectivity test" here was tiltingPlus from round 001. Pre-register hit/miss; count the pool before reading a 0;
  exec'd scripts overwrite sys.argv (save MYARGV); control as deep as the case; no-power invariants are not evidence; mark memory UNVERIFIED.
- Classifiers on literal input relations miss consequence relations (2794 first misread): classify mod the ideal.
