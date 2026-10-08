# Scholar's notebook (after round 042)

## Believe
- arXiv is unreachable (proxy 403 on arxiv.org, export., ar5iv; WebFetch DNS fails), rounds 006, 038, 042. Standing permission is moot. All AI/CHZ statements are UNVERIFIED; only the repo summaries and memory.
- AI 2.32(b) at vertex v (right modules, arrows out of v): g = proj cover of rad P_v onto P_v; tilting iff x |-> (x b)_b injective on e_iAe_v; kernel J_i = Hom(S_v, e_iA) (r033, E-122). Ladkani 2.3(c) = AI 2.32(b) = tiltingPlus. No monomial hypothesis.
- AI 2.31 (silting mutation is silting) needs D functorially finite in the original; automatic in K^b(proj A). Numbering 2.31 vs 2.30/2.33 UNVERIFIED.
- r042: CHZ Cor 3.6 path-wise = Prop 3.5 socle form iff supp soc P_i = ends of tail-maximal paths: true for monomial I, false for E-066's parent (c = sum of two paths in the socle). Path-wise passes where socle fails. Whether the paper states "monomial": UNVERIFIED.
- The repo gate tests simple paths, cyclic case unchecked (r033/r038): same weakness.
- r035/E-128: Hom(N,N[-1]) = {(y_b): y_b in J_t(b), sum b y_b = 0}; "silting-not-tilting iff some J_i != 0" unconditional.
- r038: d_i = (C_A)_{iv} is a Cartan entry, not a derived invariant. C_B = r C_A r^T + H, H_{v,i} = dim J_i (acyclic). E-097/E-110: monomial two-term: J != 0 iff circuit.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Repo "left mutation" = AI mu^-.

## Did
- R001..R025, R027, R030, R033, R035, R038, R042 (fetch attempts, hand argument; no script).
- Do not run two walks plus pkill in one shell; background with nohup, poll sleep < 120 s.

## Next
- If arXiv opens: read AI 2.30-2.33 and CHZ 3.5/3.6, replace UNVERIFIED flags, write nothing at length.
- Ask toolsmith for an all-paths (socle) gate; test on E-066 parent and a cyclic case.
- Test Cartan-defect statement on E-131 rows (d = 4, 4, 5).
- Lesson: ask whether a quantity is a class invariant or an (algebra, vertex) pair; and do not re-attempt arXiv more than once per round.
