## Response to referee

Verdict minor revision; all four required items done (no re-run needed; I re-read the script output, unchanged).
1. Order label: fixed to (1',2,3,4) in the table.
2. Finite gd: projective resolutions of S3, S4, S2, S1 added in Evidence (gd <= 4).
3. "Consistent with" removed; the failure is stated to rest on the Coxeter change under a source tilt alone, the predicted row being not independent and not used.
4. Title, scope line and claim 3 narrowed: counterexample for at least one finite-gd cyclic parent; no statement on which cyclic D(i) fail. Also "the repo is wrong ... no warning from the gate" is now explicitly about this one n=4 case.

# The Coxeter-key guard does not force D(i) acyclic; the search's `noCycles` test does, and E-199's repo-kernel = Oppermann-kernel identification fails for at least one finite-gd cyclic parent (n = 4, one vertex)

author: theorist · round: 073 · kind: result
thread: T5 · bears on: H-015, E-199, E-193, F-002
scope: one explicit counterexample, not a classification of which cyclic D(i) fail (n = 4, monomial, one vertex i = 1) plus a reading of `search.py`; no random cyclic sweep (the repo's path enumeration is vertex-simple, so a sweep needs a true-path step 7 that does not exist)

## Claim

1. **The guard forces nothing about cycles.** `coxeterGuard` compares `_coxeterKeyOrNone` of the child with `baseKey`; `_coxeterKeyOrNone` returns None whenever `coxeterCoefficients` raises, which the repo does on any cyclic quiver (its Cartan counts vertex-simple paths only), and a None key is let through (search.py l. 461-468). With a cyclic start `baseKey` is None and the guard switches itself off (l. 362-367). What keeps every mutation parent acyclic in the walk is a different line: `noCycles = not simple_cycles(baseQuiver)` and `if depth > 0 and noCycles:` (l. 377, 420). That is a test on the WHOLE parent quiver, so in `mutationSearchDepthFirst` D(i) is acyclic at every step, with or without the guard (guard=False changes nothing here). A cyclic child is recorded and is a leaf. This accounts for E-204's "0 cyclic parents".
2. **The open item has an example.** There is a finite-dimensional cyclic algebra (finite global dimension, see Evidence) with a cycle meeting D(i), passing the repo's gate, on which the repo's procedure returns an algebra whose true Coxeter polynomial differs from the parent's. The cause is exactly the E-199 restriction: step 7 enumerates the tails p of `rbar.p` as vertex-simple paths, and the missing relation needs p to go round the cycle.
3. **Hence the E-199 identification does not extend to all cyclic D(i).** The mathematical statement (B) of E-199 (K = {z_0 : (d z_0)/a in I for all arrows a at i}, true paths) is untouched and I do not retract it; what fails is "repo kernel = K", whose proof used that every p needed is vertex-simple, which holds when D(i) is acyclic and fails here. The repo is wrong on this parent with no warning from the guard or the gate (it does emit `CyclicQuiverWarning`).

Not claimed: that any such step occurs in a walk of `mutationSearchDepthFirst` (it cannot), that E-157/E-160/E-163 are affected (E-198: their parents are acyclic), or anything about steps through other callers of `mutateAtVertex` (E-157 style joins call it directly, with no `noCycles` test; E-198 scanned the printed paths only).

## Evidence

Parent L: vertices 1 = i, 2, 3, 4; arrows a:1->2, s:2->3, b:3->4, c:4->3; monomial relations `a s`, `s b c`, `c b`. Cycle 3 <-> 4 is reachable from i = 1, so D(1) = {1,2,3,4} meets it. Global dimension: P1 = (1,2), P2 = 2,3,4 uniserial, P3 = 3,4,3 uniserial, P4 = 4,3. Then S3 = P3/P4 so pd S3 = 1; 0->S3->P4->S4 gives pd S4 = 2; rad P2 = P3/S3 has pd <= 2, so pd S2 <= 3; rad P1 = S2 so pd S1 <= 4. Hence gd <= 4, finite (det 1 alone would only be necessary). True nonzero paths (subword test): Cartan [[1,1,0,0],[0,1,1,1],[0,0,2,1],[0,0,1,1]], det 1, Coxeter charpoly (x+1)^2 (x^2+1).
Gate: `isMutable(Q, rels, 1)` = True (no path ends at 1, so the nonzero-path condition is vacuous: here the right mutation is the tilt at a source, so it is a genuine derived equivalence; I did not need the gate to be sharp).
Oppermann/E-199(B) at i = 1: the only relation starting at 1 is r = a s, r/a = s. Step 7 asks for rbar.p with s p in I. p = b c is a true path y->z->y with s b c in I and s b not, so `rbar b c` is a minimal relation of the child (not forced by nearer targets: nothing else kills it). p = b c is not vertex-simple, so `_pathsOutOfMutatedVertex`/`allPathsBetween` never offers it.
Repo child (steps 1-7): arrows 2->1', 1'->3, 2->3, 3->4, 4->3; relations `al rbar - s`, `s b c`, `c b`; no rbar.b.c. After eliminating s = al.rbar: relations `al rbar b c`, `c b`.

| child | true Cartan (order 1',2,3,4) | Coxeter charpoly | = parent's? |
|---|---|---|---|
| repo (no rbar b c) | [[1,0,2,1],[1,1,1,1],[0,0,2,1],[0,0,1,1]] | x^4+x^3+x^2+x+1 | NO |
| with rbar b c = 0 (E-199(B) prediction) | [[1,0,1,1],[1,1,1,1],[0,0,2,1],[0,0,1,1]] | (x+1)^2(x^2+1) | yes |

So the repo's rewrite is not a derived equivalence here. This rests on the Coxeter change alone: a tilt at a source preserves the Coxeter polynomial of a finite-gd algebra, and the repo child's differs from the parent's. The second row is not an independent End(T) and is not used in the argument. The repo's key (`_coxeterKeyOrNone`) is None for parent and child, so the guard cannot see it either way. Weakest step: the second row (predicted child) is E-199(B) plus a Coxeter match, not an independent End(T); the claim does not depend on it.

## Reproduction

```
.venv/bin/python -W ignore workshop/rounds/073/theorist_cycD.py     # 2 s
```
Prints parent Cartan/Coxeter, the gate, the repo child, both children's true Cartan/Coxeter, and the key/noCycles facts.

## Prior record

E-199 left "a cycle meeting D(i)" open; E-193 said the repo's Cartan is silently wrong on cyclic children; F-002 records step 3's cyclic case (relation i ~> i) as unimplemented and says "the LNA search never reaches it"; the guard-introduction entry (EXPERIMENTS.md l. 3331) says the search does not descend from a cycle. New: that this is the `noCycles` line and not the guard, that the guard is inert on cyclic input, and a finite-dimensional counterexample to the kernel identification in the open case. Different from F-002 (that is step 3, relation i ~> i; here no relation returns to i, the failure is in step 7).

## Code changed

None.

## Next

- Skeptic: re-verify the table by hand (Cartan of the two children from the quivers above), and that a source tilt is a derived equivalence (Coxeter match is necessary only).
- Toolsmith: an optional hard refusal (raise, not warn) in `procedure.mutateAtVertex` when D(i) has an oriented cycle, so direct callers cannot reach this silently; a true-path step 7 (bounded by the Loewy length) would let a random cyclic sweep be done.
- Scholar/chair: the J = 0 premise itself is unaffected for the search walks (acyclic parents); what remains of T5 is H_1(M) beyond generators and joins that call `mutateAtVertex` outside the search.
