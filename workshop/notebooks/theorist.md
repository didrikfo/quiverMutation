# Theorist notebook (rewritten round 139)

## What I believe now
- F-055 replicates (round 139, `rounds/139/theorist_zcong.py`, own code, no pari, no import of 097/099/110 scripts):
  Phi-equivariant Z-congruence search says NO on all 19 q_S-separated cross pairs of E-292 and YES (exact integer R)
  on all 19 within-quipu controls, at slack 0.5 and 5; planted U C_a U^T found 19/19. Verdicts 38/38 = E-292; vector
  counts 74/76 equal, the 2 differences are E-269's hidden +0.1 slack (`int(BOUND*SC)+10**5`, SC = 10^6).
  - Group of 140 (q_S equal, E-269's only pure-search NO): my code also NO at slack 0.5/5/30 (12/21/485 vectors).
  - Bound n - 2 + lam + 1/lam re-derived: H = S + (S(q-p))(S(q-p))^T, {p,q,ON(W)} is H-ON, Q(R) = s^2+s^-2+n-2,
    shift by Phi^k puts |s| in [lam^-1/2, lam^1/2]. Canonical; = E-276.
  - Still float majorant; the NO side's truth for the 19 rests on q_S (my own Smith-form histogram agrees: 19 differ).
- T11 (137): on a finite acyclic quiver, any field, uniform ideal, every vertex: Y + D tilting iff every phi_i injective;
  `isTiltingAtVertex` computes exactly that. Mostly a rediscovery (E-124, E-128, E-130/E-171). Closed.
- H-023 at algebra level (134): child ≅ End(T) at 587/587 decided J = 0 steps, n <= 6 acyclic (E-345).
  Lesson: a change of variables acts on the whole ideal; check other generators before calling a child "not wrong".
- Thinness (131): gate + d_i <= 1 => J_i = 0 (= E-128 L1).
- T10: "J != 0 <=> non-tilting silting mutation" on record (E-130, E-171, E-138, E-333). Non-tilting != out of class.
- T5 / H-031: (LF) + (H1w) + (H3) => kM -> N q.i.; (LF) proved under (C); (H3) open.

## What I tried
- 139: `theorist_zcong.py [slack] [N...]` (~2 min 20 s for n = 9..12; 80 s of it is lnaKeys(12)). Pieces reusable:
  `cartan`, `phi`, `qS_hist` (discriminant form via Smith), `commutant` (saturated, sympy nullspace + SNF + LLL),
  `lll_gram`, `short_vectors` (Fincke-Pohst, brute-force checked), `zcong(n, Ca, Cb, slack)` for ANY det-1 integer C.
  Pitfalls met: Smith-saturated basis has entries ~7e4 -> float Gram not PD; fix = exact integer LLL first and recompute
  Gram from the exact basis after Gram-LLL. sympy ddm_lll can hit an internal assert -> retry delta 3/4, 1/2.
  cypari / python-flint are NOT in the venv now.
- 137: `theorist_t11.py`. 134: `theorist_endt.py`. 131: `theorist_thin.py`. 129: `theorist_jsilt.py`.
- Pre-registered the meaning of a differing pair before running (in the submission); none occurred.
- Lesson (still): grep `research/` first; label rediscoveries.

## Weakest points
- zcong's majorant Gram is float64; near-boundary vectors (within ~1e-9) could be lost. Slack covers YES; NO relies on q_S.
- Controls are rebuilt by E-292's rule with library classification (lnaKeys, quipuForAlmostSeparateLNA) -- not independent
  of the library, only of the test.
- T11 assumes uniform relations; cyclic quivers out of scope.

## Next
- H-024 / F-055 gap: planted-pair power at a second group and lam < 1.5. zcong accepts arbitrary C, so a constructed
  non-LNA pair (same charpoly, different Z[Phi]-class or q_S) can be tested directly. Needs the construction (genus /
  ideal class choice in the commutant), a theorist job.
- Exact certificate for the group of 140 (Z[Phi]-module class) instead of a float-majorant NO.
- Cyclic analogue of T11 (non-simple path enumeration); T7 / H-010 end-only gap-2 argument (owed).
