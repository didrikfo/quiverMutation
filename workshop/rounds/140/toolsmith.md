# Round 140 closing conference: toolsmith

author: toolsmith · round: 140 · kind: conference position (nothing run, nothing edited)

1. **Most promising question (next cycle).** H-023 at algebra level: is the child of every J = 0 step literally End(T) (presentation, not Cartan)? Every "J = 0 step is a derived equivalence" claim (F-057, E-345, E-357, the 25-of-25 joins) hangs on it, and my Hom engine (`129/theorist_jsilt.py`, `TH.test`/`TH.build`, fast at n <= 12) already computes the compositions. Build End(T) from them, compare with the child's ideal on one walk step at n = 6, then the E-332 walk.
2. **Outside the T10 thread.** S-1 / H-026: even-n self-mirror failure (E-307) and the tail end. It is the one line whose law is stated at key level only and whose zero is near-definitional (E-349); a second, independent invariant than the key would say whether K0 = min(h, n-4-h) means anything.
3. **Weakest claim.** "The J = 0 premise": that J = Hom(T,T[-1]) = 0 steps yield derived-equivalent children. Evidence is Cartan level only (E-345 n <= 6, E-357 n = 11) and E-357's three verdicts compute one quantity (E-351 forces agreement). Second: the T10 harmlessness of `refuseNonTilting` rests on n = 11 roots where no in-class refusal occurs (E-350 A has no power for the realistic loss mode).
4. **Closest to a finding.** E-357 (Hom check of the nine J = 0 steps) plus E-354 (rebuild, `sameIdeal` closed). It is a record closure, not a discovery. One check it needs: the algebra-level End(T) ~ child comparison on those 14 steps (item 1), with a planted wrong-ideal control of its own (OWN had no NO control; the referee supplied 9 of 9).
5. **Need.** Theorist: a proof or counterexample that acyclic + legal relations + dim e_iAe_v <= 1 forces End(T) = child (H-023, Anick (H3)(ii) still UNVERIFIED). Skeptic: a planted same-Cartan wrong ideal at n = 11 for the End(T) check. Scholar: Anick 1982 PDF.

## Closing ledger (my items)
- T10 depth-8 sharded run: OPEN, not approved (see recommendation).
- T10 4 unrun duals, depth 5+ at E-316: OPEN, converted to the overnight menu entry.
- `refuseNonTilting` default: SETTLED, stays off (recommendation below).
- `isTiltingAtVertex` illegal-relation refusal (E-351 exception): OPEN code item, small (a guard plus a regression on an E-080/E-112 shape); not done this cycle.
- `isMutable`/`isTiltingAtVertex` refuse on a cycle upstream of i (R-018): OPEN, owed since 127; carry to cycle 3 first.
- `mutationSearchDepthFirst` docstring (l. 322, "guard keeps a walk in one derived class" false, E-332) and canonicalKey/guard docstrings: OPEN, one edit, should be done at once; not a research item.
- Option real-cost measurement (off walk vs on walk): OPEN, minor.
- E-149's 154 lost edges with `sameIdeal`: RETIRED unless a use appears (the n = 11 `sameIdeal` capability now exists, E-354).
- E-306 78/78 corruptions independent re-verification (~40 min, 4 cores): OPEN, the one real rebuild left from the plan.

## Recommendations
- **`refuseNonTilting` default: keep OFF.** At E-316's 19 n = 9 in-class J != 0 steps it loses 93-180 in-class nodes per depth-4 walk (the step's own child in 19 of 19); it saves nothing measurable in-class at n = 11. Off is the only default with no demonstrated in-class loss. Keep it as an explicit option and put the E-350 B number in its docstring.
- **Depth-8 overnight run: do NOT approve now.** It tests only the cross-class pruning at n = 11, where E-350 A already shows no in-class refusal through depth 7; a null at depth 8 cannot change the default (the default is settled by the n = 9 loss) and cannot reach the realistic loss mode. Approve only if a theorist argument predicts an in-class refusal at depth 8, or after the cheap items (docstring, R-018 refusal) are done. If approved, a sharded run with a planted in-class loss at depth 8 as its power control.
