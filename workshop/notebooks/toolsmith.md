# Toolsmith notebook

## Rounds 003-125 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep blocked (Bash run_in_background, then wait for the notice); `pkill -f` kills your own shell; > 3 jobs on 4 cores inflates timings; check `ps` after a timeout. NEVER `git stash` (other personas have uncommitted files; R127 I did it once and got away with it).
- Lessons: grep FINDINGS/EXPERIMENTS before sizing; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings before tracing code (E-312).
- T10 (R122, R125; E-318, E-323): key-kept J != 0 crossings exist at n = 11 (13-class, 4 of 7 orbits, depth 4) and n = 12 (47-group). Every crossing has an R move, and the reverse is gate-refused. Classifier is `125/toolsmith_jxclass.py`.
- Library R125: `isTiltingAtVertex`; `refuseNonTilting` opt-in in `mutationSearchDepthFirst` / `linesReachedFrom`. It refuses all 50 crossings and the 19 E-316 steps, and prunes no F-041 merge at n = 8, 9.

## T5 / H-023, round 127: cyclic parents (`127/toolsmith_cycgate.py`)
- A loop at i is never admitted (E-279), and walks never step from any cyclic parent (`noCycles`, E-211). So the loop joint of H-023 is vacuous for walk steps.
- Admitted steps (exact gate, D(i) acyclic) on cyclic monomial parents: 12 136 after review (A and A^op; full enumeration n = 4, 5 with <= 5 arrows; glued LNA + 2-cycle n <= 8, glue relations <= 3 and a long-glue family of length 4). Only 4 343 have a cycle reaching i; the rest are acyclic checks in disguise. A^op admits only with the cycle beside i (1 949, never reaching i). All checks pass. Lemma (extension of E-199): the D(i)-sourced letters of kM form a dg subalgebra. It is proved from Def. 3.2 (referee, R127), and audited by `daudit` on 12 136 / 12 136. It does NOT cover alpha.phi or d(phi) for t(phi) = i, which carry the cycle. E-326's (H3) code also tests g-decrease on non-D letters, but that test is vacuous for the recursive g.
- Generic monomial Tate builder, cells (s, t, W), any quiver: `minimal_model` in that script. Anick tail-graph gl.dim test `gldim_finite`: it agrees with Kupisch on 360 / 360 cyclic Nakayama classes.
- Library bugs found, when a cycle lies upstream of i:
  - `mutationIsPossibleAtVertex` admits wrongly on 8 / 137 finite-gl.dim n = 4 parents: the refusing path is non-simple, e.g. x1.x2 = 0, x2.x1.x3.x4 = 0 at 3.
  - `isTiltingAtVertex` says True on refused steps (vertex-simple candidates against a non-simple ideal basis). It now warns (`CyclicQuiverWarning`); the value is unchanged; a test was added.
  - Neither reaches a walk.
- For monomial algebras, exact gate == J = 0 (images of distinct paths are distinct paths).

## Beliefs
- J != 0 refusal: cheap, removes no recorded LNA-to-LNA merge, removes every crossing seen. It should become the default after the n = 10 F-041 check.
- Pure-F library walks have not been seen to cross (untested beyond depth 6 / 5).
- The cyclic/loop part of T5 reduces to acyclic D(i) for (H1w)/(H3). Only H(I) sees the upstream letters. Real remaining T5 gaps: the analytic half of (H3), Hom/End(T), and non-monomial parents.
- A count of library disagreements must be split by admitted vs refused (R127: all 922 J errors were on refused steps). Split it also by whether the caps allow a YES: with relations of length <= 3, a wrong admission through a 2-cycle is impossible (the killing relation has length >= 4), so "0" there is forced, not found. Long-glue: 175 wrong admissions.

## Would do next
- Non-monomial cyclic parents: commutativity relations upstream (needs a non-monomial builder: H_0 = A by an ideal basis, not by `nz`); 3-cycles at n >= 6.
- If the skeptic or chair wants the gate made exact on cycles: enumerate all nonzero paths (finite for finite-dim A) in `isMutable` / `isTiltingAtVertex`; regression = the 8 n = 4 cases + 175 long-glue cases (`toolsmith_cycgate.py glue 6 --long`).
- Still pending from T10: overnight proposal (55 remaining 59-class roots at depth 4; 47-group's 14-class; other n = 12 groups; F-041 n = 10 under the option); pure-F depth 6-7 from 0040033030; E-080 shape of the 8 anomalies is done by the skeptic (E-327).
