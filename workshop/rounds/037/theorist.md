# The gate alone does not force d_i = 2 at J_i != 0: gate-admitted hand examples with (d_i, dim J_i) = (3,1) exist, but none shown to lie on a walk

author: theorist · round: 037 · kind: negative
thread: T5 · bears on: E-128, E-131, E-126

## Claim

No proof of "d_i = 2 whenever J_i != 0 on walks" can come from the gate (L1) plus local structure: the 6-vertex algebra
T1 = kQ/(1-2-5-6 = 1-3-5-6), Q: 1->{2,3,4}->5->6, v = 5, is gate-admitted with d_1 = dim e_1Ae_5 = 3 and dim J_1 = 1
(J_1 = k(p_2 - p_3)'s class, p_k = 1-k-5), and E-128's layered algebra gives (3,2). So the d >= 3, J != 0 pairs are gate-admitted.
What I could not show is that any such algebra is on a walk (derived equivalent to an LNA). T1 is not: its Coxeter key
(1,1,-5,-10,-5,1,1) is not an n = 6 LNA key (dual LNAs included), so it is in no walk. This does NOT claim d >= 3 with J != 0
is impossible on walks; E-131's absence there stays empirical, and is now a statement about the derived class, not the gate.

## Evidence

Hand (`hand`, 2 s): T1 gate True, (d,J) at i=1 = (3,1), key not an LNA key. T2 (two commutation relations) (3,2), key not LNA.
Enumeration (`enum 6`, 3 s): cores = three paths i -> v of lengths in {(2,2,2),(2,2,3),(2,3,3),(3,3,3),(1,2,2),(1,2,3)}, 1 or 2 out-arrows,
every nonempty set of pairwise commutations after the out-arrows (up to 6 vertices, no extra attachments needed at n = 6),
legal, acyclic: gate refuses 7, gate-admitted without (d >= 3, J != 0) 27, gate-admitted with d >= 3 and J != 0: 183, of which 167 have
no LNA key (so are not on a walk, by the Coxeter-polynomial invariant) and 16 have an n = 6 LNA key.
Those 16 have key coincidence only, which is not derived equivalence (E-125: keys do not decide). I mutation-BFSed each (gate-admitted
vertices, reduction, canonicalKey dedup, about 230 algebras each, 54 s for 400): no LNA or dual LNA appears in any of the 16 mutation classes
(all 16 appear to be one class: the counts agree to within levels); inside that class gate-admitted rows with (3,1), (4,1), (5,1), (5,2), (6,2),
(7,1) occur, so if that class were derived equivalent to an LNA, E-131's d = 2 pattern would fail. Not shown; the BFS is a bounded sample
(frontier not exhausted), so "no LNA" is a miss in a bounded search, not a proof of non-equivalence.

Why no proof is available from L1 alone: L1 says J_i is a proper subspace of e_iAe_v; any proper subspace not containing a path
class is allowed, and T1 realises one of dim 1 in a 3-dim space. A proof of d = 2 would need a derived-class invariant that excludes
the fan-with-commutation shape (e.g. a Coxeter-polynomial / Euler-form constraint), which is the T3 invariant question.
The step I am least sure of: reading "key not an LNA key" as "not on a walk" -- it needs the key set at n = 6 to be the full set of
Coxeter polynomials of LNAs (it is `lnaKeys(6)`, built from allOfLength(6) and duals) and gl.dim finite for T1 (it is acyclic).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/037/theorist_d3.py hand          # 2 s
timeout 10m .venv/bin/python workshop/rounds/037/theorist_d3.py enum 6        # 3 s, writes theorist_d3_hits_n6.json
timeout 10m .venv/bin/python workshop/rounds/037/theorist_d3.py bfs file:workshop/rounds/037/theorist_d3_hits_n6.json:0 150   # about 40 s each
```

## Prior record

E-128 (L1, layered (3,2) example), E-131 (285 rows all d = 2, 28 rows d >= 3 all J = 0), E-125 (keys uninformative). New: explicit
gate-admitted (3,1), and an enumeration showing that d >= 3 with J != 0 is common among gate-admitted fans, with the only
obstruction to walks being derived-class membership. Not a rediscovery as far as grep shows.

## Code changed

None (new script `workshop/rounds/037/theorist_d3.py`, reads the round 033 script's helpers). No tests run.

## Next

- toolsmith/experimentalist: walk from the 28 d >= 3 rows' algebras and from LNA classes at n = 6, 7 is complete in principle at n = 6:
  close the BFS of every n = 6 LNA derived class (small) and tabulate (d,J) exhaustively; that settles n = 6 outright.
- theorist: find which invariant separates T1 / the 16 key-coincidences from LNA classes (T3).
