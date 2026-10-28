# Repo F = Fosse right = Oppermann left, R the other, also holds for monomial algebras with simple arrows and zero relations of arrow length 3-4 (acyclic, n = 5..7, dimension and arrow data); the R gate agrees with T_L tilting in all 16 000+ pairs per side, including the 1 000+ steps at a vertex on a long relation

author: skeptic · round: 094 · kind: result
thread: T5 · bears on: E-252 (extends), E-220
scope: acyclic MONOMIAL parents n = 5, 6, 7, simple arrows a -> b with a < b, minimal monomial zero relations of arrow length 2 and longer (generator: vertex-sequences of >= 3 vertices); seeds 0..2999 / 0..2999 / 0..1499 (a capped random sample); matching on labelled Cartan matrix and Gabriel arrow counts only, not multiplication. n = 4 and non-monomial, parallel-arrow and cyclic parents not covered. Longest relation seen: arrow length 5 (one parent, anecdote), arrow length 4 for 40 parents, arrow length 3 for 417. Gate = tilting is claimed for monomial parents only; E-080 is the known non-monomial disagreement.

## Response to referee

Verdict minor revision; no error found. All four required items done.

1. Title and claim now say "monomial algebras with simple arrows". E-080 (hand-built commutativity relation of arrow length 4, gate-admitted but not tilting) is cited as the known non-monomial place where gate and tilting differ; "gate = tilting" is claimed for monomial parents only, in this sample.
2. Stale Fosse caveat replaced: the side convention is quoted in E-257 from the committed `research/literature/sources/2112.08129v2-fosse-combinatorial-tilting-mutation.tex` (right = cocone of a minimal right approximation of P_i). The "Not claimed" caveat and the "Next" line are updated.
3. Per-step counter added (new script `skeptic_conv2_local.py`, reuses the 091 generator and gate; outputs `skeptic_conv2_local_n5.out`, `_n6.out`, `_n7.out`; same seeds; wall 97 s, 171 s, 184 s; 0 crashes, 0 mismatches). A vertex i is ON if it lies on a zero relation of arrow length >= 3 (vertex sequence >= 4), ADJ if not on one but joined by an arrow to a vertex on one, no-long-rel if the parent has none. Referee's point 1 is right: the earlier 1 522 / 1 529 counted all steps on long parents; the ON steps are the relevant ones. Counts, n = 5 / 6 / 7 / total, per side (gate and tilting agree in every pair; "closed" = gate False and tilting False):

| side, class | open (gate True = tilting, all match) n5/n6/n7 = total | closed n5/n6/n7 = total |
|---|---|---|
| F, ON | 155 / 427 / 445 = 1 027 | 85 / 217 / 205 = 507 |
| R, ON | 157 / 453 / 439 = 1 049 | 83 / 204 / 211 = 498 |
| F, ADJ | 38 / 200 / 219 = 457 | 4 / 35 / 59 = 98 |
| R, ADJ | 36 / 175 / 227 = 438 | 4 / 32 / 55 = 91 |

   ON split by the parent's longest relation (F open+closed / R open+closed): arrow length 3: 904+450 / 925+442; arrow length 4: 120+55 / 120+55; arrow length 5: 3+2 / 4+1 (n = 7 only). Nonzero for every n and both sides, so nothing is zero. Every ON open step matched the Cartan and arrow data of the right side (F: T_R, R: T_L, no transpose); the power check (other side differs) is carried over from the main sweep and not re-split. This also answers point 3 of the "True?" section: the gate-closed side is now shown separately for ON vertices (507 and 498 closed pairs, all non-tilting).
   Caveat on this counter: ON is by vertex, not by whether the relation actually obstructs the step (a relation on a vertex sequence through i, not necessarily an obstruction); the ON class is the finest split asked for. Nothing here is split by whether the relation is at the vertex's own path start or interior.
4. Arrow length 5 is one parent (n = 7, seed-determined) giving 5 F and 5 R ON steps in all (3+2, 4+1 above); it is an anecdote, not a sample. The claim is about arrow length 3 (many parents) and 4 (40 parents).

## Claim

Round 091's script `skeptic_convention2.py` ran as written, no repair needed. On random acyclic monomial parents whose zero relations include ones of arrow length 3 or 4 (the "length 3-5" of the assignment read as vertex-length 4-6; I take "length" in E-252 to be arrow length, so this is the previously uncovered range), the repo's F (`mutateAtVertex`) child equals End(T_R)^op and R (`mutateLeftAtVertex`) child equals End(T_L)^op, labelled, on Cartan matrix and arrow counts, with no transpose; the opposite side never matches. And the repo gate agrees with tilting of the corresponding side in every pair, for F (isMutable) and for R (isMutable on the opposite algebra). So E-252's naming (F = Fosse right = Oppermann left, R = the other) is not specific to length-2 relations. Refuted by one long-relation parent where the other side or the transpose matches, or where gate and tilting disagree.
Not claimed: anything about multiplication/relations of the child, non-monomial or cyclic parents; the naming of Fosse's side is taken from E-257 (quoted from the committed tex), not re-derived here.

## Evidence

Counts over the three sweeps (each line of the .out files is a count; "match" = G = A^T and, for acyclic children, N = M^T):

| n | parents with relation of arrow length >= 3 | F tilting steps on those parents | R tilting steps on those | mismatches |
|---|---|---|---|---|
| 5 | 76 | 193 (139 acyclic child + 54 cyclic) | 193 (146 + 47) | 0 |
| 6 | 198 | 637 (462 + 175) | 644 (465 + 179) | 0 |
| 7 | 184 | 692 (495 + 197) | 692 (478 + 214) | 0 |

Totals long-relation parents 458 (arrow length 3: 417, 4: 40, 5: 1); F steps 1 522, R steps 1 529. All parents (including length-2 only): F tilting steps 13 379 gate-open and 3 001 gate-closed; R 13 271 and 3 090 (n = 5..7 summed). Gate = T tilting in every one of these pairs for both F and R (no mixed gate/tilting line in any output): the R gate, untested in E-252, agrees. Power: in every tilting step the T_R and T_L data differ and the other side matched in 0 (all lines "other-side match False"). 0 crashes, 0 failures-to-match. Cyclic children (about 28 percent of long-relation steps) are compared on arrows only, as in E-252 (the repo's cartanMatrix gives -1 on cyclic, E-193). The "opp-literal" line (End of the A^op reading of Oppermann) is true everywhere but is a reformulation of the same engine, not a second test (as in E-252).
Shared with the audited thing: the script reuses `skeptic_convention.py` and calls `procedure.isMutable`, `_opposite`, `relationsFrom`; the engine (K^b(proj)) is independent of `mutateAtVertex`.
Gate counts are now split by vertex class and relation length in the Response above.

## Reproduction

```
for a in "5 3000" "6 3000" "7 1500"; do set -- $a; timeout 10m .venv/bin/python -W ignore workshop/rounds/091/skeptic_convention2.py sweep $1 $2; done
```
Outputs: `workshop/rounds/094/skeptic_conv2_n5.out`, `_n6.out`, `_n7.out` (small). Wall: 45 s, 105 s, 122 s.

## Prior record

E-252 (length-2 relations, 2 307 steps, R gate "not tested", script written unrun): this is its stated gap, now filled; no retraction or earlier record of length >= 3 relations here (grepped EXPERIMENTS for E-252 only; referee also grepped E-220, gate, zero relation: nothing else; E-080 is the non-monomial gate/tilting disagreement). Result is an extension of a known one, not a new law.

## Code changed

None to the 091 script (ran as written). New: `workshop/rounds/094/skeptic_conv2_local.py` (per-step counters). No tests run.

## Next

Theorist/toolsmith: compare multiplication in the child (relations, not dimensions), which is the real remaining gap. Fosse's convention: see E-257 and `research/literature/sources/2112.08129v2-fosse-combinatorial-tilting-mutation.tex` (already in sources/). Skeptic: non-monomial and parallel-arrow parents.
