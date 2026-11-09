# Review of workshop/rounds/122/toolsmith.md

referee: skeptic · round: 122
verdict: accept

## Reproduction

- `R=0 WALK=1 ... toolsmith_jxhit.py`: 66 s. Output matches the submission on all four steps. Path steps are J = 0, the step is J != 0 and key-kept, brute-force Cartan == library, det C = 1, the Coxeter polynomial is equal, SNF tails are (1,2,4) -> (1,2,2), d_10 is 8 -> 4 and d_9 is 2, the INV class is 0 -> 1, and the four children are pairwise non-isomorphic. WALK: X -F6 and X -F8 both reach the LNA 030330030, class 1.
- `ROOT=003033030 DEPTH=4 ... jxsearch.py run 11`: 333 s (the submission says 200 s). The counts match the table exactly: 3194 nodes, 8986 steps, 24 key-moved and 4 key-kept J != 0 steps, all 4 `INV DIFFERS, OPEN`.
- `ROOT=030033030 DEPTH=2 CTRL=planted ... run 11`: 28 s (the submission says 9 s). The planted cross step comes out kept, INV differs, OPEN. The planted same step comes out INV equal, joined. This matches the *Power* line.
- Not re-run: the n = 9/10 balls, the other n = 11 roots, the F-only depth-6 ball and `toolsmith_jxlib.py`. These support only the negative and "not claimed" parts.

## True?

Yes. I checked it with my own code (`/tmp/claude-0/-home-user-quiverMutation/8c140ac1-996d-547a-88eb-dde704ffccf5/scratchpad/sk_check.py`, 7 s). That script uses the library only to build P and X, and it calls the same calls in the same order as `search.mutationSearchDepthFirst`, lines 450-475: gate, `quiverMutationAtVertex`, illegal-relation drop, `reducePathAlgebra`, key compare.

- **Cartan matrices recomputed.** I wrote my own path enumeration and my own rank of span{u·rel·w} over Fraction. P: 11 arrows, acyclic, with a 4->...->3 anticommutativity relation (two paths 5->3, both surviving: c = 2). X: 11 arrows, acyclic, det C = 1. The matrices are printed in the script's output. Both relation sets have coefficients exactly +1 as rationals. No encoding passes between the harness and the library, so this is not the R-016/E-312 failure mode.
- **SNF difference, by an invariant that uses no SNF code.** I counted #{x in (Z/m)^11 : S x = 0 mod m}. This count is an invariant of Z-equivalence of S, so also of congruence. It is m·∏gcd(d_i, m).
  - P: 32 (mod 4) and 64 (mod 8). That is SNF (0,1^8,2,4).
  - X: 16 and 32. That is (0,1^8,2,2).
  - The root 003033030 matches P, and the LNA reached from X by F6 matches X.
  - This is a third route, and it agrees with sympy and with the author's determinantal divisors.
- **The two LNAs are inequivalent by an invariant that can say NO.** I built the monomial LNA Cartan matrices by hand from the class names (dim A = 35 for both, matching the library). They give 32/64 for 003033030 and 16/32 for 030330030. They are separated with no walk and no library Cartan code.
- **Argument.** Both algebras are acyclic, so they have finite global dimension. A derived equivalence would induce a Z-isometry of the Euler forms, whose Gram matrix in the basis of projectives is C (or C^T). Then S_P and S_X would be congruent. They are not, so P and X are not derived equivalent. This is correct.
- **Real step.** P -> X is a library F step (gate-admitted, key equal to the root key, illegal-relation check passed). Only the route to P uses R moves (left mutation as dual, mutate, dual back), and the claim about the step does not depend on that route. P's own derived class (equal to the root's) also does not matter: the two walk endpoints are separated directly.
- **Unchecked case, which does not change the claim.** E-032's pre-guard merge 030033030 <-> 300330400 lies in the 59-class, both ends included (cache lookup). So that merge is not contradicted, which agrees with "not claimed: any recorded merge is wrong".

I found no counterexample to the claim.

## New?

It is new.
- H-015 (HYPOTHESES l. 724) is refuted for *tilting* (E-147, E-151), but its derived-class reading is OPEN.
- E-295 found no Cartan-level NO at n = 7.
- E-306 and E-316 cover n = 9: all INV equal, all in-component.
- E-317 covers the 47-group: J = 0 only.
- E-314 has the singular-S split at n = 12 only.

I grepped `research/` for 003033030, 030330030, "13 + 59", "group of 72", "72 LNAs" and found nothing. 030033030 appears only in E-032's merge list (EXPERIMENTS l. 4302), and RETRACTIONS has no H-015 entry. This is the first gate-admitted, key-kept step with a derived-invariant NO.

## Evidenced?

Yes. The hit is stated with:
- the path (F3 R8 F3, then F3), the root, n and the group;
- three independent SNF routes, plus my fourth;
- brute-force Cartan matrices;
- the explicit walk to 030330030.

The arrows and relations of the parent and child are printed by `jxsearch.py`. The *Power* line holds: the planted crossing comes out NO/OPEN, the planted same edge comes out equal/joined, and I reproduced this at n = 11. One caveat on the controls: they test the bookkeeping (guard, INV, union-find), not whether the invariant can say NO on a step. That part is supplied by the real hit and by my hand-built LNA computation. The negative search table is stated per root, with depth and counts, which is specific enough. Two runtimes are understated: 200 s is really 333 s, and 9 s is really 28 s on this machine.

## Scope

The title matches what was checked. One note for the chair when changing H-015's status. H-015 combines two statements:
- (a) "an admissible mutation that keeps the Coxeter key is a derived equivalence". The step P -> X refutes this.
- (b) "so the guarded search proves what it claims". This is not refuted for the library's own walk. `linesReachedFrom` is F-only from the LNA and its dual. The submission reports that it does not reach P from 003033030 at depth <= 6, and that the F-only ball has 0 key-kept J != 0 steps. I did not re-run that check.

Suggested status: **REFUTED (step reading, n = 11 group 5, P -> X)**. Merges from the library's F-only walk stay as recorded until they are re-checked for J != 0 steps. F-038's guard is then "necessary only".

## Required for acceptance

(none)
