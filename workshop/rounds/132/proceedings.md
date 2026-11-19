# Round 132 -- proceedings (conference, phase 3 of cycle 2)
Six position statements (`rounds/132/<persona>.md`). Step 0.5: round 131 left no questions; nothing to settle. Merge of origin/main: already up to date.

## Ledger
- **F-057 written** (F-candidate E-323/E-332): key-guarded walks cross derived classes at n = 11 (pure F) and n = 12 (F/R) through a key-kept J != 0 non-tilting step. The blocker named at 128 (a Hom-level / non-SNF check) was done in E-333/E-334, so it is written at current scope; open points in its body: `sameIdeal` of E-332's end node against 403033030 and a Hom check of its nine J = 0 steps (skeptic), cost of the option on non-thin balls.
- **R-018 written** (R-candidate "gate decides tilting"): library gate over-admits on cyclic monomial parents (E-340), library scope, n = 5-7; the blocker (gate vs Hom) was done. Not a statement about LNA walks.
- H-015 status gets a pointer to F-057; H-032 status gets "written up as F-057". Both stay as they were otherwise.
- H-025 stays SUPPORTED: E-339 rests on one referee (skeptic and scholar both call it the weakest claim); its blocker is an object-level nu^m check on a simple module (scholar).
- T4 / H-029: closed (SUPPORTED, E-335); locality still unproved, no assignment.
- Not changed: H-010 (SUPPORTED), H-023/H-031 (OPEN), H-026 (SUPPORTED, odd n), H-017 (dormant).
- Library/docstring consequence still owed: `isMutable`/`isTiltingAtVertex` refusing when a cycle reaches i (toolsmith, R-018).

## Where the personas agree and disagree
Experimentalist, skeptic, toolsmith, maverick all name the same gap: E-337's "no join lost" is forced by thin balls (E-338), so the option's cost on balls with dim e_i A e_v >= 2 steps is unmeasured; theorist instead asks whether the child equals End(T) as an algebra (the premise under every "J = 0 => joined"); scholar asks what H-025's n = 10-12 groups are. Weakest claims named: SNF(S)/Phi necessity chain (maverick, scholar), E-339 (skeptic, scholar), J = 0 => derived equivalent child (theorist), E-337 reading (experimentalist, toolsmith). Disagreement: toolsmith offers an n = 11 depth-4 background run (4.2 h CPU); skeptic wants the flip threshold registered before any run.

## Proposed agenda for phase 4 (rounds 133-136)
1. **T10: option cost on non-thin balls** (toolsmith, maverick; skeptic first). First question: pre-register a threshold for "option on loses joins"; then one counted pass over walks that hold d >= 2 steps (E-332 depth-6 walk, E-316 F-010 joins, n = 9 merges), reporting the denominator of non-vacuous steps, with the planted E-318 edge as control; also E-332's `sameIdeal` + nine-step Hom check. No n = 11 depth-4 background job until the counted pass says it is needed.
2. **T5 / H-023: child = End(T) for J = 0 steps** (theorist; toolsmith normal-form emitter; skeptic planted wrong child). First question: End(T) multiplication table against the child on n <= 6 acyclic parents.
3. **T3/T8 / H-025** (scholar). First question: simple-module nu^m growth on kA_15/rad^3 (m = 144, 288, 432) and the type of the n = 10-12 groups; the owed correction to literature/2112.15587. Skeptic or theorist second-reads the E-339 step.
4. **T7 / H-010 overlap 2** (experimentalist; theorist argument on enterability). First question: engine at k = 4 on 10-20 clean (8,4)(10,3) cells.
5. **T6 / H-017b** (skeptic, pre-register only; run only if the baseline has power) -- dormant since 119. Maverick's grading question (why no LNA ball reaches a cyclic node) as a side question for the theorist if free.

## Questions for the steering committee
None. Standing: PDFs of Nikulin 1979, Happel-Seidel, BGS, Anick 1982 would remove UNVERIFIED tags (not needed to proceed).

## Decisions taken for the steering committee
- Agenda is the chair's (item list above). The toolsmith's n = 11 depth-4 overnight-sized run is not approved; `refuseNonTilting` default stays off; no OVERNIGHT.md change.
