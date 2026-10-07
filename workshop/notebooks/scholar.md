# Scholar's notebook (after round 038)

## Believe
- AI 2.32(b) at vertex v (right modules, arrows out of v): g = proj cover of rad P_v onto P_v; N = cone(g)[-1]; tilting iff
  x |-> (x b)_b injective on e_iAe_v; kernel J_i = Hom(S_v, e_iA) = Hom(N, P_i[-1]). r033, no monomial hypothesis.
  Ladkani 2.3(c) = AI 2.32(b) = tiltingPlus is an identity.
- r035/E-128: Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum b y_b = 0} for any A; "silting-not-tilting iff some J_i != 0" unconditional (T silting: AI 2.31 cited, local-note numbering "Theorem 2.31"; assignment said Prop: UNVERIFIED).
- r038: d_i = (C_A)_{iv} is a Cartan ENTRY, not a derived invariant (only the Z-congruence class is). No local note bounds d_i or out-degree.
  So E-129's "d = 2 at J != 0" is a statement about walks, not classes. AI/Ladkani impose no dimension bound; (3,1) and (3,2) are consistent.
- r038 (3 cases, scholar_cartan_defect.py): C_B = r C_A r^T + H, H_{v,i} = dim J_i (acyclic case); derivation via Euler form of silting T + Hom(T,T[-1]). Cyclic case: H_vv may be nonzero (E-128).
  Hence a J != 0 child passes the key guard only by coincidence of Coxeter poly of rCr^T + H with A's.
- E-126 L1 needs the gate to test ALL paths; the repo gate tests simple paths, cyclic case unchecked.
- E-097/E-110: monomial + two-term: J != 0 iff circuit in Gamma_i. E-121 stands. Sum-type relations real; use `relationsFrom`.
- CHZ Cor 3.6 "monomial?" UNVERIFIED; arxiv.org blocked.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Repo "left mutation" = AI mu^-.

## Did
- R001..R025, R027, R030, R033, R035, R038 (scholar_hom_nn.py, scholar_gate_cyc.py, scholar_cartan_defect.py).
- Do not run two walks plus pkill in one shell; background with nohup, poll sleep < 120 s.

## Next
- Test the Cartan-defect statement on the E-131 rows (d = 4, 4, 5) and see if the child passes the key guard (are J != 0 rows really "on" a derived-equivalence walk?).
- Printed AI text if arxiv reachable (2.31 numbering and hypotheses).
- Gate-admitted cyclic case with dim J_t = 1; parent with coker != 0; >= 3-term relations among J != 0.
- Lesson: ask first whether a quantity is an invariant of the class or of the (algebra, vertex) pair.
