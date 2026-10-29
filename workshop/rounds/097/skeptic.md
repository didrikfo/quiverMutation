# In a random sample of acyclic monomial/commutativity parents (n = 4..6, search depth 3), every step the search takes satisfies E-235's hypothesis; the degree-1 relation-arrow case of E-248 is never taken (discarded by `isIllegalRelation`), its degree-2 analogue is taken and checks out

author: skeptic · round: 097 · kind: negative (a planned hit did not occur; coverage check passes)
thread: T10 (and T5) · bears on: H-015, E-235, E-248, E-249
scope: sample only (450 search roots, 855 + 514 census roots), n = 4..6 · random acyclic monomial/commutativity parents with a planted relation omega: i ~> t parallel to an arrow i -> t (and a planted Anick chain a b c plus arrow i -> t) · dimension data only (no multiplication) · F steps only (the search walks no R) · J = 0 steps only among the steps taken · random sample, not exhaustive

## Response to referee

1. Title/claim: done. Title and Scope now state the sample (450 search roots; 855 + 514 census roots), n = 4..6, depth 3. The Claim no longer says "every" outside the sample.
2. Roots: the reviewer's "300 per row" is per n, not per row: the census argument 300 is per n, so `4 6 300` is 3 x 300 = 900 attempts and `5 6 300 pat2` is 2 x 300 = 600 attempts. The table's "300 roots each" meant per n (now worded so). Not every attempt gives a root (the generator rejects non-minimal ones): printed with `skeptic_rootcount.py` (same seed, same rng use) 855 roots from 900 attempts and 514 from 600. The search ran exactly 150 roots per n = 450. Claim and table now give these numbers; the totals 956 and 563 PAT pairs come from these 855 and 514 roots.
3. Done. "dim End(T)" for the 271 cyclic-child steps now reads "equals dim End(T) up to arrow length 12 (truncated count, stable between caps 10 and 12; a heuristic, not a proof)".
4. Done: PAT and DEG2PAT are defined in the body (section Definitions), as coded in `skeptic_omega_census.py`.
Not re-run: the three reproduction commands (the referee reproduced them; no code changed). The line-by-line table totals were not re-added by the referee; they are the script's printed counts.

## Definitions (as coded)

- PAT at (parent, vertex v): some relation r of the parent has source v, and its target t equals the target of some arrow v -> t (a relation omega: v ~> t parallel to an arrow). Code: `pat = any(source(r) == v and any(a[1] == target(r) for a in arrowsOutOf(q, v)) for r in rels)`.
- DEG2PAT at v: arrows a, b, c with a: v -> x, b: x -> y, c: y -> t, zero relations ab and bc in the (minimal) relation set, and an arrow v -> t. Code: `deg2pat` in `skeptic_omega_census.py`.
- Taken step: `procedure.isMutable` gate passes, no relation of the (unreduced) child satisfies `arrowPaths.isIllegalRelation`, and the Coxeter key of the reduced child is None or equals the parent's key.
- Census roots: 300 attempts per n; rejected attempts give no root.

## Stated beforehand (what a hit and a miss mean)

Written before the runs. The J = 0 premise has two links: (L1) the repo's steps 1-7 give End(T) (E-193, E-199, E-235); (L2) T is tilting/generating, i.e. J = 0 at that step. This round tests only L1, on the steps the search takes.
- HIT = a step the search takes (gate, illegal-relation test, key guard all passed; parent acyclic) where (H) of E-235 fails, or where the independent dim End(T) differs from the repo child. Meaning: at that step "child = End(T)" is not established, so every "J = 0 path / join" that contains such a step is a statement about the rewrite, not about derived equivalence; the premise would be unsupported there, and the E-157..E-163 style joins would need the step list audited for it. A hit does NOT refute J = 0 itself.
- MISS = no such step. Meaning: L1 is covered for every taken step of the sample. It says NOTHING about L2: dim End(T) = Cartan(child) is how E-181 showed End cannot see J, and all taken steps here are J = 0 anyway. A miss cannot raise the premise's credibility beyond "the rewrite is End(T)"; the generation gap stays open.

## Claim

On random acyclic parents (450 roots for the search, 855 roots (900 attempts, n = 4..6) and 514 roots (600 attempts, n = 5..6) for the two censuses; relations of arrow length 2-3, zero or commutativity, minimal presentations, one omega || alpha or Anick chain planted) the real `mutationSearchDepthFirst` (depth 3) takes 5 136 steps, and the censuses (steps at the roots, gate / illegal test / key guard recomputed) 2 119 more. In every taken step of the sample the parent is acyclic, so E-235's hypothesis (no loop at i of any degree) holds, and the independent dim Hom_K(T_a, T_b) (round-050 code) equals the repo child's Cartan matrix: 0 mismatches in 5 066 search steps and 1 848 census steps with an acyclic child, and for the 271 taken steps with a cyclic child (leaves) a truncated path count (arrow length <= 10 and 12, stable) equals dim End(T) up to arrow length 12 (a truncated count, heuristic, not a proof). The E-248 pattern (a relation out of i whose target t has an arrow i -> t) occurs in 1 519 (parent, vertex) pairs and in none is the step taken: step 5 always writes the relation (alpha*, rbar) along t -> i* -> t, a vertex-repeating path, and `isIllegalRelation` (search.py l. 441-451, run before `reducePathAlgebra`) discards the step. So the degree-1 loop of E-248 never enters a taken step. The degree-2 analogue (Anick chain a b c from i with ab = bc = 0 and an arrow i -> t, i.e. Ext^3(S_i, S_t) != 0) IS taken (238 steps, all with a cyclic child) and passes the dim check.
Not claimed: that J = 0 holds; anything about multiplication in the child; R steps; parents deeper than the sample reaches; n > 6; that discarding the PAT steps is harmless (it loses edges, see Next). Refuted by one taken step with a loop at i, or one taken step whose dims differ from End(T).

## Evidence

Search run (`skeptic_omega.py 4 6 150 3`, 450 roots, seed 97): 5 136 taken steps (361+1 / 1402+7 / 3303+62 acyclic + cyclic child at n = 4/5/6), every parent acyclic and with a minimal relation set (also at depth 2, 3 descendants), 0 PAT steps, End = child on all 5 066 acyclic-child steps (all J = 0; the gate refused the J != 0 ones). Census (all vertices of the roots, gate / illegal / key guard recomputed as the search does; `skeptic_omega_census.py`):

| census | pairs with PAT | PAT taken | taken, acyclic child, End = child | taken, cyclic child, trunc dims = End | DEG2PAT taken | mismatches |
|---|---|---|---|---|---|---|
| standard, n 4-6, 855 roots (300 attempts per n) | 956 | 0 | 1 180 | 23 | 9 | 0 |
| planted Anick chain, n 5-6, 514 roots (300 attempts per n) | 563 | 0 | 668 | 248 | 229 | 0 |

All 1 519 PAT pairs are pre-reduce illegal with a cyclic child (a few are also gate-refused or key-dropped). The non-PAT gate-refused pairs are J != 0 with End = child dims (E-181 again, not new; 0 of them taken). A first run with non-minimal roots (zero relation a b c over an already zero b c) showed 4 PAT steps whose raw child is cyclic but whose reduced child is acyclic with End = child and J = 0; the search drops those too (illegal test is before `reduce`). They were an artefact of non-minimal roots (the generator now rejects them; descendants are reported as found, all minimal), so I do not count them, but they show the discard rule is a code property, not a theorem.

## Reproduction

```
timeout 10m .venv/bin/python -W ignore workshop/rounds/097/skeptic_omega.py 4 6 150 3        # 73 s, workshop/rounds/097/skeptic_omega_search_out.txt
timeout 10m .venv/bin/python -W ignore workshop/rounds/097/skeptic_omega_census.py 4 6 300 97       # 15 s
timeout 10m .venv/bin/python -W ignore workshop/rounds/097/skeptic_omega_census.py 5 6 300 97 pat2  # 10 s (outputs *_census_out_*.txt)
```

## Prior record

E-249 closed "no loop at i" for the 25 failing parents and states the loop classes for one step (quiver level, random graded quivers); E-235 needs only (H); E-223: step 7 runs only at `noCycles` nodes; E-193 / E-199: step-7 kernel = red(Z_1 e_i), acyclic D(i). New here: (1) the E-248 pattern is, for the taken steps of the real search, empty because of the `isIllegalRelation` discard, a fact about code; (2) the degree-2 analogue is taken and its leaves (cyclic children) were never dim-checked before, here 271 are; (3) End(T) is the round-050 independent code on non-line, non-LNA parents with relations (E-181 covered LNAs and 38 parallel-arrow steps). grep of `research/` for "illegal" together with E-248/E-249 found no recorded link between them (E-033 introduced the `continue`).

## Code changed

None in `quivermutation/`. New scripts `workshop/rounds/097/skeptic_omega.py`, `skeptic_omega_census.py` (they exec the round-050 `skeptic_tilt.py` functions; no tests touched).

## Next

- Theorist: the E-249 remark "loops of degree >= 1 at v != i" needs no more work for the search (H_0 each step, E-235 hypotheses hold), but the degree-2 case sits in E-235's (H) only via the cited Lemma A (E-237, E-240); the 238 taken DEG2PAT steps are test material for the Lemma A hypotheses on a non-hereditary parent.
- Toolsmith: the pre-reduce `isIllegalRelation` prunes every omega || alpha step (1 519 pairs in the sample). Are these valid J = 0 edges the walk loses (cyclic child, leaf)? Size a count of PAT steps whose reduced child is acyclic or whose dims equal End(T) with J = 0 (my 4 artefact cases suggest the reduced-child rule matters). This is a NO-power question for T10, not a coverage one.
- Not done: relations of the child compared by multiplication; non-monomial coefficients; R steps; n = 7 (the n = 7 failing parents are covered by E-249).
