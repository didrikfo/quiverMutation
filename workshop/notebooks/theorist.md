# Theorist notebook (rewritten round 037)

## What I believe now
- J_i = Hom(S_v, e_iA) (E-122); gate = "no single nonzero path in J" (E-114); L1: dim J_i <= d_i - 1 (E-126). L1 is only a restatement of the gate.
- Round 037: gate does NOT force d_i = 2: hand T1 (1->{2,3,4}->5->6, 1-2-5-6 = 1-3-5-6, v=5) is gate-admitted with (d,J) = (3,1); layered m=3 gives (3,2). T1 has non-LNA Coxeter key, so not on a walk. Any proof of d=2 must use derived-class structure, not the gate.
- 16 gate-admitted (d>=3,J!=0) six-vertex fans have an LNA key, but 16 bounded BFSes (~230 algebras each, one class) found no LNA: key coincidence only (E-123). Not a counterexample, not evidence for the pattern either.
- L2 (acyclic: Hom(T,T[-1]) = sum J_i) conditional on AI 2.31, not compared with the text. Cyclic L1 untested.

## What I tried
- 037: `rounds/037/theorist_d3.py` (hand / enum 6 / bfs). 034: `theorist_dimji.py`. 031: `theorist_t3.py`, `theorist_layers.py`. 029: `theorist_circuit.py`, `theorist_keys.py`. 026: W exact match.

## Next
- Close (exhaustively) the n = 6 LNA derived classes and tabulate (d,J): settles n = 6 for d>=3.
- n = 7 enumeration with extra attached vertex (`enum 7`, size it with a dry count first).
- Find an invariant excluding fan-with-commutation shapes from LNA classes (T3).
- Blind spots: a BFS miss is not non-equivalence; a key match is not equivalence.
