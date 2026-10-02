# The 42 n = 8 class-0 out-degree-2 rejections are walk-reached by construction, parents stay in the class, and they are the E-103 two-out-arrow kind (minimal presentations), not a new mechanism

author: skeptic · round: 025 · kind: result
thread: gate-vs-tilting / long square (E-100, E-103) · bears on: E-100, E-103, 023 scholar review

## Claim

Rerunning the round-023 walk (n = 8, class 0, 200 s cap) gives the referee's 42 gate-admitted rejections, all out-degree 2, no long square (7 877 and 8 004 algebras in two runs, 42 both times). (1) Reachable: yes, trivially. They are parents inside the guarded BFS, which keeps only algebras whose Coxeter key equals the class base, so all 42 parents have the class-0 Coxeter polynomial (checked: 42 of 42 equal the base key, 42 distinct canonical keys). The E-103 off-walk examples that left the class (Coxeter polynomial in no LNA/dual class, n = 6) are a different phenomenon: those were hand-built, not walk parents. (2) Presentation: not an artefact of redundancy. In all 42 no relation of `alg.rels` lies in the ideal of the others (so E-103's "redundant long relation" explanation does not apply), and gate admits, `tiltingPlus` is False, kernel dimension is 1, no single path is in the kernel. (3) Mechanism: every one has the vertex v with two out-arrows and at least one relation ending through v into each out-arrow (per-arrow relation counts (1,2): 29, (2,1): 7, (1,3): 4, (3,1): 1, (1,4): 1); this is E-103's "two out-arrows each carrying a long relation" kind. So it is the same mechanism that was found off the walks at n = 6, 7, and what is new at n = 8 is only that the guarded walk reaches it (none at n = 5..7, E-100). It does not claim: a kernel element x was exhibited (I did not extract one), that out-degree 2 with relations on both arrows always rejects, or anything about classes 1, 2.

Caveat on "genuine": in sample 1 (below) the relation `4513 = 4573` sits next to `451 = 0`, so it is a monomial zero relation `4573 = 0` written as a commutativity. `alg.rels` is still irredundant, but the quiver-with-relations has a shorter presentation; the walk's presentation is not minimal in the sense of a monomial/commutativity split. The reject itself comes from `kerdim` over `relationsFrom`, independent of this.

## Evidence

Walk statistics (200 s, 42 rejections, all `outdeg 2`):

| per-out-arrow relation count | parents |
|---|---|
| (1,2) or (2,1) | 36 |
| (1,3), (3,1), (1,4) | 6 |

Relation lengths in the parents are 3..6, no monomial single-path kernel (`mono` False 42 of 42), 8 vertices each. Sample parent (v = 3, arrows 1>3, 2>4, 3>6, 3>8, 4>5, 5>1, 5>7, 6>8, 7>3; relations 138=0, 2457=0, 451=0, 4513=4573, 5136=5736, 57368=0, 7368=738) has two out-arrows 3>6 and 3>8; the first relation through 3>8 is `738 = 7368`, through 3>6 is `5136 = 5736`. Two more samples are printed by the script (`skeptic_n8_out.txt`).

## Reproduction

```
timeout 10m .venv/bin/python -u workshop/rounds/025/skeptic_n8.py 8 200 0 > workshop/rounds/025/skeptic_n8_out.txt   # about 5 min; counts depend on the cap
```
The script execs the round-023 `scholar_longsquare.py` helpers (`row`, `kerdim`, `longSquare`), keeps the walk identical, and collects the hits. `skeptic_n8.py` prints class membership, per-arrow shapes, 3 samples and the redundancy test.

## Prior record

E-103 (two-out kind, 19 of 26 off-walk rejections; its Limits cite the round-023 scholar n = 8 count). E-100 (n = 5..7 walks: all rejecting parents out-degree 1 with long square). Not in RETRACTIONS. New: the n = 8 walk count is now attributed to the E-103 kind, the 42 parents keep the class, are irredundant. The count 42 is cap-dependent (8 988 vs 7 877 vs 8 004 algebras); nothing says 42 is a total.

## Code changed

None in the library. New `workshop/rounds/025/skeptic_n8.py`.

## Next

Theorist: prove or refute "step-7 map fails iff exists x != 0 in e_a A e_v with x*beta in I for every out-arrow beta" (E-103 conjecture), which would explain out-degree 2 without long square. Experimentalist: n = 7 class 1 or 2 and n = 8 classes 1, 2 for the same count; extract x for the 42 and its length. Toolsmith: a minimal-presentation reducer (zero relation plus commutativity) so shape tests such as `hasLongSquare` stop depending on how relations are written.
