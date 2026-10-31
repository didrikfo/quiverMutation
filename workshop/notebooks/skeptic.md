# Skeptic's notebook (after round 102)

## Believe now
- r102 (T3/T8, literature): E-269's bound = n - 2 + lam + 1/lam (+ slack) exactly, any lam > 1 (hand derivation: H = alpha^2+beta^2 on the Salem plane, H^-1 S = swap + identity; checked on 450 LNA matrices, err 1e-13). So lam in (1.1,1.5) is an input gap, not a bound gap; bound is n+0.01..n+0.17 there. Literature (secondary only): Lehmer 1.17628 smallest known; lists complete to degree 44 (2409.11159 intro, citing FGR 1999 / MRW 2008, UNVERIFIED primary); 2008.10965 Prop. 1: 1.230391. Census (n<=12, 10 Salem lam < 1.5) matches those two; 8 others unchecked. Mossinghoff web lists unreachable (DNS). Submission workshop/rounds/102/skeptic.md.
- r099 (T3/T8): q_S has no power in periodic key groups (321+2 split = E-190 classes); its power is in non-periodic n = 11, 12 (E-265). maverick_cartanZ imports cypari, absent here; r102 scripts are numpy/sympy only.
- r097 (T10/T5): sampled steps (n 4..6) satisfy E-235 (H); tests L1 (dims) only. r094: F/R convention holds for acyclic monomial parents (E-258). r091: Fosse convention F = mutateAtVertex = Oppermann LEFT; child = End(T)^op.
- r087 H1 confirmed on LNAs; r085 2400230 ~ 3304000 certified (6+8 no-meet = move-set gap); r082 E-229 real J != 0; r074 70 failing edges = 60 keys = 12 classes; r069 E-147 13,16,9 are cap artefacts; r046 Cartan det/SNF no power. Unit of evidence = orbit.

## Tried
- r102 skeptic_bound.py, skeptic_salemcensus.py. r099 qs_periodic. r097 omega. r094 convention2. r091 convention. r087 homotopy. r085 cert. r082 kept/join. r074 failkids. r069 c1_13. r053/050 skeptic_tilt.py. r046 collect/inv/power/iso.
- Lessons: generator artefacts; `pkill -f` kills own shell; isolated vertices vanish; ap.cartanMatrix wrong on cyclic quivers; `| tail` buffers; print label is not a test; null where a cheaper invariant already splits is not evidence; arXiv API queries need --data-urlencode (quoted phrases silently return nothing); only arxiv.org/export reachable, not wayback/cecm.

## Not done
- Primary Lehmer/Boyd/FGR/MRW not read; Mossinghoff fixed-degree lists; 8 census lam unmatched; qfminim not rerun (no cypari); cross-pair isometry check of the identity (only R = Phi^k tested); Child multiplication; R steps; non-monomial coefficients; n >= 7 random parents; illegal-relation discard; q_S orbit controls; Fosse sentence unquoted.

## Next
1. Scholar-level: primary sources for "complete to degree 44". 2. q_S-separated pairs at n = 12 vs real E-187/E-205 moves. 3. Count PAT steps lost by pre-reduce discard. 4. Engine on the 70 key-kept J != 0 edges. 5. skeptic_convention2 cases (non-monomial, cyclic, parallel arrows).

## Habits
- Pre-register hit/miss; vacuity by definition; control as deep as the case; independent route; no-power invariants are not evidence; say what is shared with the audited thing; check what a script really calls; grep EXPERIMENTS and research/literature before "not recorded"; reject generator artefacts first; mark memory as UNVERIFIED.
