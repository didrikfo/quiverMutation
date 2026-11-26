# Review of workshop/rounds/139/scholar.md

referee: experimentalist · round: 139
verdict: minor revision

## Reproduction

- `timeout 10m .venv/bin/python -I workshop/rounds/139/scholar_scherotzke.py 14`: 1.2 s, 11 of 11 OK, 10 of 10 controls refused, same as the submission.
- `timeout 10m .venv/bin/python -I workshop/rounds/139/scholar_isometry.py 7`: 1 min 50 s. The output is identical to `scholar_isometry.out` apart from the timing lines (4 FOUND, 6 NONE).
- One case further: a scratchpad script checks the stated n ≥ 5 witness shape directly, without searching. The shape is columns −(e_1+..+e_{n-1}), e_n, e_1−e_n, e_2..e_{n-2}, with V^T E_L V = E_D exactly and det V = −1. It holds at **every n = 4..20** (17 of 17), including n = 4, where the search returned a different witness. Run time is a few seconds.
- Sources re-fetched into new scratchpad directories: 0906.4987v1, 0902.0518 **v1, v2 and v3**, and 1310.1907v2 (LaTeX, plus the compiled PDF through pdftotext).

## True?

Every number and statement I checked matches.
- **0906.4987v1:** Thm 5.7 is at l. 1145 (section 5 counter: 5.1 l. 817 ... 5.5 Example l. 1049, 1→2→3→4 with I = ⟨αβγ⟩, 5.6 l. 1103, 5.7 l. 1145). The statement is quoted verbatim. The proof's last line (l. 1210) reads "kA_n", the typo the submission reports. The field is "a field k" (l. 113).
- **0902.0518v3:** Thm 4.9 is at l. 728 ("bounded components" ⟹ only component). Thm 4.15 is at l. 1016 (⟹ A ~der kD). A commented-out `%\begin{lemma}` at l. 829 does not count. 4.15's proof cites BGS Thm A, B.
- **1310.1907v2:** I matched 4.6 (l. 858), 4.7 (l. 883), 4.8, 4.9 (l. 906), 7.4 (l. 1245, Happel 5.12, alg. closed k, basic, simply connected), A.4 (l. 1397), A.7 (l. 1444, `co:complete`), A.8 (l. 1453), Lemma 4.3 (l. 796) and Thm 3.15 (l. 720). In the PDF, "Theorem 7.4", "Corollary A.7" and "Corollary A.9" appear as numbered. The A.7-versus-8.7 trap is real: a counter-order count gives S8.7.

Gaps, none fatal:
1. **Which version of [S] the citation means.** 0906.4987v1's bibliography gives [S] as "archiv: ." with no id and no version. In **0902.0518v1** (Feb 2009) the numbering differs: 4.9 is a corollary (finite component ⟹ semisimple) and 4.15 is `inj proj` (a Z[A_∞] statement). Under v1, "[S, 4.9]" and "[S, 4.15]" would point at the wrong results. **v2** (updated 2009-06-26 17:43, four minutes before 0906.4987v1's source timestamp) has 4.9 = `bounded` and 4.15 = `inj pro`, as v3 does. The submission checked v3 only. The match holds, but it should be stated against v2, the contemporaneous version.
2. **The dependency list is incomplete.** 4.15's proof (v2 l. 965, v3 l. 1027) also uses Vossieck [V, 1.1] (discreteness) and Happel's book [H1, I.5.5, I.5.6] (same tree class ⟹ derived equivalent kQ ~ kD). The label "resting on BGS A/B (UNVERIFIED primary)" understates this. The route rests on BGS A/B plus Happel I.5.5/I.5.6 (and Vossieck's definition), all unread.
3. **The control label is wrong.** "kA_n = N_n(2)" (scholar_scherotzke.py docstring, and "kA_4..kA_9" in the power line): N_n(2) is the radical-square-zero algebra, which is derived equivalent to kA_n but is not kA_n itself. The control is still valid, because the K0 test only sees the derived class. Only the wording is wrong.
4. **Claim (3) omits route 3's own hypotheses.** Thm 7.4 needs k algebraically closed and A simply connected. The literature file notes that "tree algebras are simply connected" is unchecked, but claim (3) does not carry that caveat.

## New?

- I grepped `research/` for "Scherotzke", "0906.4987" and "0902.0518": nothing found. The Scherotzke route and its BGS dependency are new.
- The result itself is not new. E-352 (R137) already proves Λ(n,n-1) ~ kD_n for all n ≥ 4 (Brüstle 1.1/1.2, PD Euler form), and the submission says so.
- E-352 also "retired" BGS. This route revives BGS as a dependency, which the submission flags.
- Hubery-Krause 1310.1907 is cited in E-308 and F-056, but there was no literature file. The new file is untracked (`git status`: `??`).
- The isometry route (3) is weaker than E-352 (finite n, an extra simply-connected hypothesis). It is an independent cross-check, not a new result.

## Evidenced?

Yes. Line numbers, versions and verbatim quotes are given, and every unread primary is marked (BGS, Happel 1987/1988, Igusa-Schiffler). The two gaps are the [S] version (item 1) and the missing Happel-book and Vossieck dependencies (item 2). The power statement for the isometry controls is honest: it notes they are decided by root count alone.

## Scope

The headline "Scherotzke Thm 5.7, all n ≥ 4" matches the source. Claim (3) is correctly limited to n = 4..7 as searched. The referee check above puts the explicit witness shape at n = 4..20, so (3) may be widened to "the explicit witness of (3) verified at n = 4..20". A general-n Gram proof stays with the theorist, as the submission proposes.

## Required for acceptance

1. In claim (1) and the evidence table, note that 0906.4987v1's [S] carries no arXiv id. State that 4.9 and 4.15 match in 0902.0518 **v2** (2009-06-26) and v3, and that **v1 numbers differently** (v1 4.9 = Cor., 4.15 = `inj proj`).
2. Extend the honest label in (1) and in "Next" to list all unread primaries: "resting on BGS A/B, Happel's book I.5.5/I.5.6 (and Vossieck 1.1), UNVERIFIED primaries".
3. Reword the control "kA_n = N_n(2)" as "A_n-type control N_n(2) (rad² = 0, ~der kA_n)", in the claim's power line and in the script docstring.
4. Add to claim (3) its hypotheses (k algebraically closed; N_n(n-1) simply connected, not checked in a source). Optionally widen the witness range to n = 4..20, citing this review.
