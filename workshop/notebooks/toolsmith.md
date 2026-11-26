# Toolsmith notebook

## Rounds 003-130 (short)
- `orbitCensus`/`OrbitsTask` in `batch.py`. Traps: foreground sleep is blocked (use Bash run_in_background with an until-loop, or `timeout 9m bash -c 'until ...'`); `pkill -f` kills your own shell; more than 3 jobs on 4 cores inflates timings; NEVER `git stash`. python-flint is NOT in .venv (use sympy SNF). Create the scratchpad dir before redirecting into it.
- Lessons: grep FINDINGS/EXPERIMENTS before sizing a job; "no hit" needs a positive control; `canonicalKey` is LABELLED (E-164); check encodings first (E-312). LNA names have n-2 digits.
- T10 R122-125 (E-318, E-323): key-kept J != 0 crossings at n = 11, 12; `isTiltingAtVertex` and `refuseNonTilting` are in the library.
- T5 R127 (E-329): the gate and `isTiltingAtVertex` are wrong when a cycle is upstream of i (owed, R-018).
- T10 R130-134 (E-337, E-341, E-344): 0 refusals at depth 3+3; 3 key-kept refusals at E-332's crossing roots, all at the horizon step.

## T10 R137 (E-350): depth-7 option walks, `137/toolsmith_t10walk.py`
- n = 11 13-class: L_in = 0 because no in-class refusal occurs there. The referee's (B) at E-316 n = 9: loses 93-180 in-class nodes per depth-4 walk. Default stays off; the formal decision is at 140.
- `shard`/`merge` mode exists for depth 8 (overnight-sized, not approved).

## T10 R139: Hom check of F-057 (a)'s nine J = 0 steps, `139/toolsmith_homj0.py` (+ `.out`, `_prereg.txt`)
- Three verdicts per step: LIB `isTiltingAtVertex`; HOM `129/theorist_jsilt.test` (imported unchanged, full chain maps over Q); OWN (new: own path enumeration, ideal = span of p r q, Hom(T,T[-1]) = sum ker(e_iAe_k -> sum e_iAe_t) + Hom(Y,Y[-1])).
- Result: 9/9 walk-b J = 0 steps plus 5/5 walk-a steps are tilting by all three, with End0 ~ child Cartan. 0 disagreements in 139 step evaluations. Runs in 4 s.
- Controls said NO ×3: walk J != 0 steps (2), E-316's 19 n = 9 steps (H[-1] = J = 2/3/4), E-080. Also 103 J = 0 path moves to E-316's parents are tilting.
- Honest framing: the outcome is forced by E-351 (acyclic, no illegal relation => LIB = Hom), so this is a record closure, not a discovery. OWN's approximation check is trivially true (no NO control of its own).
- The assignment said "at E-316". I read it as E-332's walk (F-057's open item) and used E-316 as controls. That reading is stated in the submission.
- `theorist_jsilt.py` is a reusable Hom engine: `TH.test(A, k, child)`, `TH.build(arrows, rels)`. Fast at n <= 12.

## Beliefs
- On acyclic parents, J = Hom(T,T[-1]) at every recorded step I have touched. The remaining gap for "J = 0 steps are derived equivalences of the walk's nodes" is H-023 (child = End(T)); it is checked at Cartan level only.
- In the n = 11 13-class, the option removes only cross-class subtrees at depth <= 7 (E-350 A). At n = 9 F-010 it removes in-class nodes (E-350 B), so default-on is not harmless in general.
- d >= 2 alone is not J != 0 (E-338 needs a non-monomial piece).
- Library walks are F-only; mixed F/R joins are not the option's business until a mixed library walk exists.

## Would do next
- Algebra-level End(T) ≅ child on one walk step: build End(T) from the Hom engine's composition and compare presentations (H-023, small n first). This would upgrade "Cartan level".
- T10 leftovers: depth 8 sharded (OVERNIGHT proposal), the 4 unrun duals (003033030, 023033030, 203033030, 223033030 dual; ~1-8 min each), depth 5+ at E-316.
- Measure the option's real cost: an off walk WITHOUT the spy vs an on walk, same root, alone on the machine.
- Owed: `isMutable`/`isTiltingAtVertex` should refuse on a cycle upstream of i (R-018), with a regression; an illegal-relation refusal in `isTiltingAtVertex` (E-351 exception).
- Older pending: E-149's 154 lost edges with `sameIdeal`.
