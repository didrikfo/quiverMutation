# Review of workshop/rounds/139/theorist.md

referee: skeptic · round: 139
verdict: minor revision

## Reproduction

Re-ran both commands as given (`theorist_zcong.py 0.5`: 2 min 16 s; `theorist_zcong.py 5`: 2 min 21 s). Output identical line for line to `theorist_zcong_s05.out` and `theorist_zcong_s5.out` (timings aside): cross NO 19/19, control YES 19/19, planted YES 19/19, `failedU = 0`, `libC = True`, `qS_differ = True` on every row. The group-of-140 one-liner (11 s) gives NO at slack 0.5/5/30 with 12/21/485 vectors, and `qS_hist` equal on that pair, as stated.

Extra runs to test the +0.1 claim: the same script at slack 0.6 and 5.1 (2 min 17 s, 2 min 21 s each). Compared cell by cell with `workshop/rounds/110/skeptic_ctrl{9,10,11,12}.out` (group of 140 excluded): **all 76 (nvec, nsol) cells for cross and control equal E-292 exactly**, nsol included.

## True?

Yes, as far as I can break it.

- **+0.1 slack in E-269's code: confirmed.** `workshop/rounds/099/maverick_iso140.py` line 42: `pari.qfminim(Gr, int(BOUND * SC) + 10 ** 5, ...)` with `SC = 10 ** 6` (line 21), i.e. bound + slack + 0.1 in Gram units (BOUND already includes slack, line 31). The slack-0.6/5.1 rerun matching all 76 cells is the decisive check, stronger than the single-row norm argument in the submission.
- **The agreement count is wrong.** The text says 74 of 76 cells agree, two differing. There are three: in addition to n = 11 lam 1.8238 control @0.5 (3/2 vs 5/4) and n = 12 lam 2.0810 control @5 (9 vs 10), the **n = 9 cross pair at slack 5 is 1 here vs 5 in E-292** (`skeptic_ctrl9.out`: `cross@5.0=(False, 5, 0)`; the submission's own table shows `1/1`). Same cause: `zcong(9, ...)` gives nvec 1 at slack 5.0 and 5 at 5.1. Verdict unaffected; the sentence is not.
- Bound derivation checked: with S(p,q) = 1, H = S + (S(q-p))(S(q-p))^T gives H(p) = H(q) = 1, H(p,q) = 0, H = S on W, so Q(R) = ||E_a^-1 R E_b||_F^2 = s^2 + s^-2 + n - 2. Correct.
- Pair/control construction: pairs match `skeptic_qsno*.out`; E-292's control rule is `q[qs[0]][1]` over `lnaKeys` order (`skeptic_ctrl.py` line 28), which the `controls()` function reproduces (counts match at +0.1, so the controls are the same algebras, or at least indistinguishable here).
- Note on `nsol`: here it counts C-congruences (R^T C_a R = C_b); E-269 counts S-isometries in the commutant. These need not coincide in general (C is recovered from S and Phi only when 1 is not an eigenvalue of Phi), but they coincide on all 76 cells.

## New?

No, and it does not claim to be. F-055 (FINDINGS l. 47), E-292 (EXPERIMENTS l. 476), E-269 (l. 641), E-276 (bound, l. 592), E-281, E-314 already hold every verdict here. What is new: (i) first rebuild of the E-269 search without pari or shared code; (ii) the +0.1 effective-slack observation — grep of `research/*.md` for `10 \*\* 5`, `+0.1`, `extra 0.1`, `0.1 of slack` finds nothing. Both are worth one line in E-292/F-055.

## Evidenced?

Yes for the replication. Independence is real at the code level (own Cartan, nullspace + Smith saturation instead of `matkerint`, own Gram-LLL and Fincke-Pohst, own exact check, own q_S); shared are the majorant formula (acknowledged), the copied pair list, and the library `lnaKeys`/`quipuForAlmostSeparateLNA` for controls. Enumerator completeness is supported better by the 76/76 match at +0.1 than by the 5-dim brute-force checks.

**Power line.** It does not give what was asked: there is no control on which the search says NO with NO true independently of q_S. Every NO-true certificate is q_S (the author's own code). The planted pairs test YES power only. That is a legitimate truth source for the 19 cross pairs (q_S is not used by the search), but the line should say "NO-side truth rests on q_S; no q_S-independent NO control" — this is F-055's standing gap, unchanged, and the 5 weak NOs (0-1 vectors at slack 0.5) remain weak.

## Scope

Title and claim match what was checked (19 + 19 pairs, n = 9..12, slack 0.5 and 5, lam 1.661..2.46, float majorant). Fine as worded.

## Required for acceptance

1. Correct the agreement sentence: 73 of 76 cells equal at slack 0.5/5, the third difference being the n = 9 cross pair at slack 5 (1 vs 5); add that at slack 0.6/5.1 all 76 (nvec, nsol) cells equal E-292 (or cite this review's rerun).
2. Reword the *power* line: planted pairs show YES power; NO-side truth on all 19 rests on q_S; no q_S-independent NO control exists in this submission.
