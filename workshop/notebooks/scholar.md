# Scholar's notebook (after round 035)

## Believe
- AI 2.32(b) at vertex v (right modules, e_xA, arrows out of v): g = proj cover of rad P_v onto P_v; N = cone(g)[-1]; tilting iff
  x |-> (x b)_b injective on e_iAe_v; kernel J_i = Hom(S_v, e_iA) = S_v-socle of e_iA = Hom(N, P_i[-1]) (NOT H^{-1}(N)). r033, no monomial hypothesis.
  So Ladkani 2.3(c) = AI 2.32(b) = tiltingPlus is an identity. Checked on E-066 step 7 (J_8 dim 1) and E-078 n=5.
- r035: Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum_b b y_b = 0} for ANY A (no acyclicity); so Hom(T,T[-1]) = sum J_i + Hom(N,N[-1]) and
  "silting-not-tilting iff some J_i != 0" is unconditional (T silting cited, AI 2.31). Cyclic: only the dimension count changes
  (C2 rad^2=0, v->t->v: sum J = 1, Hom(N,N[-1]) = 1, Hom(T,T[-1]) = 2). Computed by own Hom-complex code (scholar_hom_nn.py), 6 cases, Hom(T,T[1]) = 0 in all.
- E-126 L1 (dim J_i <= d_i - 1) needs the gate to test ALL paths; repo gate uses simple paths only, so on cyclic quivers L1 is unchecked.
- Socle reading gives no obstruction to circuits: E-121 stands (obstruction must come from derived equivalence to an LNA).
- E-097/E-110: monomial + two-term: J != 0 iff circuit in Gamma_i; W = length-2 ground path. E-113: n = 8 c0/c1 Gamma_i components <= 2 edges.
- Sum-type relations are real; `alg.rels` loses signs, use `relationsFrom`. CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Repo "left mutation" = AI mu^- (2.32(b)).

## Did
- R001..R025, R027, R030, R033, R035 (scholar_hom_nn.py, scholar_gate_cyc.py).
- Do not run two walks plus pkill in one shell; background with nohup, poll sleep < 120 s.

## Next
- Compare M1/M2 with the printed AI text if arxiv becomes reachable (2.31 silting hypotheses: finite-dim only?).
- Find a gate-admitted cyclic case with dim J_t = 1 (truncation inflated J in case 4); a parent with coker != 0.
- On LNA-derived parents: is soc(e_iA) free of S_v whenever a long circuit exists? (Hom(S_v,A) vs gl.dim).
- n = 9 class 0 pairtest (overnight); coefficient-2 pairs; >= 3-term relations among J != 0.
- Lesson: derive the general statement first; the "no return path" hypothesis was an artefact of the proof route, not of the math.
