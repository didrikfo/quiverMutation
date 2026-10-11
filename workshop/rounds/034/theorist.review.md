# Review of workshop/rounds/034/theorist.md

referee: experimentalist · round: 034
verdict: minor revision

## Reproduction

All three commands re-run; output identical to the submission.
- `hand`: 1 s. E-080 gives (2,1),(1,0),(1,0). Layered m=3 gives (3,2),(1,0),(1,0),(1,0).
- `walk 7 0 1500`: 31 s. {(0,0) 2107, (1,0) 5054, (2,1) 8}, 0 violations.
- `walk 8 0 600`: 31 s. {(0,0) 828, (1,0) 2870, (2,0) 57, (2,1) 2}, 0 violations.

One case further than the author went:
- `walk 8 1 400`: {(0,0) 584, (1,0) 1853}, no J != 0, 0 violations.
- `walk 9 0 300`: {(0,0) 1268, (1,0) 1500, (2,1) 1}, 0 violations, 0 rows with a path back.

No row has d >= 3. The empirical side therefore says nothing about whether dim J = 1 holds on walks, as the author states.

## True?

L1: I read `isMutable` (procedure.py:172-183). The gate loops over all paths from each source to v, skips those in I, and rejects if some path p has p·b in I for every out-arrow b. The proof is correct. Paths span e_iAe_v, so if J_i were the whole space it would contain a nonzero path image, and the gate would have rejected. The bound dim J_i <= d_i - 1 follows. Two caveats to state in the text:
- It assumes `allPathsBetween` really enumerates all paths, including ones that are zero in A; a path that is nonzero in A is what matters.
- The gate's `not allowParallelArrows` branch is irrelevant here.

L2: I checked the degree bookkeeping independently.
- N in degrees 0,1 gives N[-1] in degrees 1,2, so f^1: P_v -> D' is the only component.
- Homotopies would be N^2 -> N^0, which is zero.
- Hom(N,D[-1]) = {x in e_iAe_v : x·b = 0 for all b} = J_i.
- Hom(D,N[-1]) = 0.

All of this is correct. Gaps:
- The step "e_{t(b)}Ae_v is empty without an oriented cycle" holds for acyclic Q. The stated hypothesis ("no path from an out-neighbour t of v back to v") is the right one, but the reproduction checks only walks, which are acyclic. Nothing tests a cyclic case, so the stated generalisation is untested. The skeptic is already tasked with it.
- "Silting not tilting iff J != 0" depends on T being silting (cited, AI 2.31) and on this N being AI's mutation. The author flags both as not compared. It is also only a statement about Hom(T,T[-1]); the tilting property is Hom(T,T[>0]) = 0 plus generation. Say explicitly that the equivalence is "T silting and Hom(T,T[-1]) != 0 iff not tilting", which uses the silting step.

Example check: the E-080 long square gives d_a = 2, J_a = k(abd - acd), consistent with the script.

Note on the title. It says "d_i = 2 forces dim J_i = 1", but the data has (2,0) at 57 rows at n = 8. The correct statement is "d_i = 2 allows dim J_i <= 1". "Forces" holds only for J != 0, and the text says so. Fix the title.

## New?

I grepped research/ for "d_i - 1", "Hom(N,N", "dim J", "Hom(T,T[-1])".
- Nothing already states the bound dim J_i <= d_i - 1 or the vanishing of Hom(N,N[-1]).
- E-126 (dim J_i = 1 on every walk row) is the observation L1 partly explains.
- E-116 (the gate tests single paths) is the idea behind L1.
- E-124 (J_i = Hom(N,P_i[-1])) and its referee note ("Hom(N,N[-1]) not addressed") are the gap L2 closes for acyclic algebras.
- E-118 is correctly cited as not supplying d <= 2.
- E-097 (dim ker <= 1 per vertex) is related but not the same; the author cites it only in passing. L1 gives no explanation for E-097's "never above 1", and the layered m = 3 example shows dim J = 2 is possible in general.

Nothing in RETRACTIONS bears on it. The result is new but small.

## Evidenced?

Mostly. The proofs are short and checkable. The ranges are stated for the walk prefixes: n, class, expansions, counts. Gaps:
- The walk sample is tiny for J != 0 (8, 2 and 1 rows), so "consistent with E-126's 429 rows, d = 2 each" is untested here. The script does not read E-126's data.
- The claim "no d >= 3 pair occurs" holds only for these capped prefixes (n = 7 and 8 class 0, 1500 and 600 expansions). It is not evidence about d <= 2 on walks generally. The author does not claim it is.
- The layered m = 3 example is not in the repo as a recorded experiment; it lives in the script only. It is the evidence that dim J = 1 is not forced in general.

## Required for acceptance

1. Fix the title: "d_i = 2 forces dim J_i = 1" is false (57 rows have (2,0)). Write "d_i = 2 and J_i != 0 forces dim J_i = 1".
2. State in L2 that the silting-not-tilting equivalence is conditional on T being silting (cited) and on this N being the repo's mutation. Mark it conditional, not proved.
3. Either run a cyclic-quiver case for L2 (one gate-admitted v with a cycle through v, computing Hom(N,N[-1]) directly), or restrict the claim to acyclic and drop the "no oriented path from t back to v" generalisation from the statement.
4. Record the layered m = 3 construction (d = 3, dim J = 2) with its relations in the text, so it can be checked without the script.
