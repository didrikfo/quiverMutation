# Skeptic's notebook (after round 129)

## Round 129 (T10, E-323 blocker), submission rounds/129/skeptic.md
- Hom-level (E-327 engine, rounds/129/skeptic_hom323.py, 3 s): the step P·F3 at P = 0040033030·F3 R9 F3 has Hom(T,T[-1]) = 2,
  so it is non-tilting by Hom, not only by injectivity. The dual mirror R3 gives the same. P·R5 is also J != 0 (Hom -1: 1) and changes the key.
  35 gate-admitted steps at the root and at P: J, isTiltingAtVertex and Hom agree, 0 disagreements. Hom dims of T = the library child's Cartan, even when not tilting.
  P's relations have E-080's shape (one arrow 3->8, e_6Ae_3 = <6-7-3, 6-9-3>, relation (6738)+(6938)).
- Invariants (skeptic_inv323.py, 5 s): ranks over F_p of f(Phi), 23 f, p <= 13. 15 (f, p) separate the pair. Every one is a multiple
  of (x-1) mod 3 or 5, i.e. SNF(S) = coker(Phi-1) seen again. x+1 (skew) and the rest agree. So there is no second argument, as in E-324 at n = 11.
- The NO chain is: derived equivalence => C_B = r C_A r^T (Ladkani 3.5 + AI 2.27 + Rickard, E-308; E-272) => SNF(S) and Phi-conjugacy invariants equal.
  Nikulin 1.13.3 is only a YES-side (class = genus) input. The conference-128 blocker's premise misread F-055's wording.
- Controls: E-080 Hom NO, abd=acd YES, 33 J=0 Hom YES; invariants YES on 5 Hom-tilting walk steps and on 10 r C r^T.
  No independent NO control for the invariant exists (every certified NO goes through the same chain). The det-2/3 plants are caught only by unimodularity.

## Believe now
- E-323's crossing is as solid as the Cartan necessity chain. The step is genuinely non-tilting (Hom). Asking for a non-Cartan invariant
  blocks forever (HH of an LNA is trivial, E-324). Recommend the F-candidate with "inequivalent = Cartan not GL_n(Z)-congruent" in its scope.
- Unverified: whether End(T) is isomorphic to the library child (only Cartan dims match). If it is, the library computes End of a silting-not-tilting
  complex, which is the theorist's question (agenda item 1).
- r126: H-030 explained (path vs element in the kernel); the injectivity test agrees with Hom on 298 337 random cases and now on 35 walk cases.
- r123: no recorded merge uses a J != 0 step except audited ones. r121: n = 9 key-kept J != 0 children stay in the parent's J = 0 component.
- r113: HH/Hom-only join tests lack power. r099: q_S has no power in periodic groups. r094: F/R convention. r091: F = Oppermann LEFT. The unit of evidence is the orbit.

## Not done
- Hom check of the other 49 E-323 steps (asked of the toolsmith; skeptic_hom323.check, seconds each). It also tests the R-candidate on walk parents.
- refuseNonTilting default decision (phase 5): regression = F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused.
- n = 7 "25 of 25": E-160's 3 paths not Hom-tested. E-032 n = 10/11 merges need --witness reruns. H-029 widths 8-9.
- Are the nn / consequence circuits reachable from an LNA walk? Needs walk rows.

## Next
1. If the toolsmith does not: run check() over toolsmith_jxclass_out.jsonl rows (50 steps) myself. Pre-register: 50 of 50 Hom NO, Hom(T,T[-1]) total 1-2.
2. End(T) vs child as algebras at P·F3: compare quiver and relations of End(T) (arrows = irreducible maps; sc.invariants gives Ar).
3. Referee F-candidate text: scope must name the chain (E-308), not Nikulin; say "at two sites".

## Habits
- Grep the TEST first ("injectivity" = tiltingPlus = isTiltingAtVertex, none of them Hom). Pre-register hit or miss. Count the pool before reading a 0.
- exec'd scripts overwrite sys.argv (save MYARGV). Make the control as deep as the case. No-power invariants are not evidence. Mark memory UNVERIFIED.
- Before attacking a blocker, check its premise against E-entries (r129: the Nikulin premise was already answered by E-272/E-308).
- Classifiers on literal input relations miss consequence relations: classify mod the ideal. Rank-mod-p invariants: check what f reduces to mod p.
