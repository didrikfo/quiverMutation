# Every step the guard admits at n = 6, 7 is a genuine tilting mutation, and the exact criterion (Ladkani 2.3(c)) catches the one step the guard exists for

author: scholar · round: 001 · kind: result (with a proposal)
thread: T5 · bears on: H-015, F-016, F-038, F-047, E-032, E-034

## Claim

For every mutation step that the gate (`procedure.isMutable`) admits, from every LNA of length 6 (depth 6) and
length 7 (depth 4) and their opposites, walking distinct algebras: (a) Ladkani's iff-criterion for
`T^+_k` to be a tilting complex (arXiv:1001.4765 Prop. 2.3(c), with linear combinations of paths)
holds, and (b) the Cartan matrix of the rewritten, reduced algebra equals `r^+_k C (r^+_k)^T`
exactly (Prop. 3.6), i.e. is Z-congruent to the parent's. 61,718 steps (29,822 at n=6, 31,896 at n=7),
0 failures. The guard never fires in this range (the unguarded walk admits no step the guarded walk
refuses), which is F-038's "changes no answer at n <= 7" measured on a second invariant.
The one place the guard fires, the E-032 ALARM path, is reproduced: at step 7 the gate admits, and
2.3(c) says NOT tilting, the Cartan congruence fails, and the key moves -- three independent signals agree.
It does NOT claim H-015 is proved: no guard-admitted non-tilting step was found, but none was
expected at these sizes (see Evidence), and the deep region where the guard matters is not covered.

## Evidence

Second-invariant candidates surveyed (literature in `research/literature/`, R-008):

| candidate | computable here? | verdict |
|---|---|---|
| Avella-Alaminos-Geiss | needs gentle | closed, R-008, F-019 |
| Hochschild cohomology | yes (Bardzell) | `HH*(A)=k` for every LNA (2312.14699); blind on the start, not computed on mutated non-monomial nodes |
| tau-periodicity / Coxeter periodicity | yes | only rules out, silent below n=10 (F-046 area, math/0611201) |
| hereditary form (H-015's suggestion) | only at relation-free nodes | rare on a path; weak |
| Z-congruence / SNF profile (F-047) | yes | sound, but Cartan-only |
| **Ladkani 1001.4765 Prop 2.3(c)** | yes, a rank computation | **an iff for the step itself**; nothing in `research/` FINDINGS/EXPERIMENTS/HYPOTHESES uses it (grepped) |

Prop 2.3(c) is stronger than any invariant here: if `T^+_k` is a tilting complex the step is a derived
equivalence (Rickard), whatever the polynomial says. So it tests the guard's *claim* directly, per step,
instead of hoping a second invariant separates.

| walk | algebras expanded | admitted steps | guard passes | tilting (2.3c) | Cartan = r C r^T | not tilting |
|---|---|---|---|---|---|---|
| n=6, depth 6, unguarded | 9,476 | 29,822 | 29,822 | 29,822 | 29,822 | 0 |
| n=7, depth 4, unguarded | 8,988 | 31,896 | 31,896 | 31,896 | 31,896 | 0 |
| n=6, depth 5 guarded / unguarded | 4,582 | 14,142 | 14,142 | 14,142 | 14,142 | 0 (identical both ways) |

The gate's own refusals are not counted. Illegal-relation and cyclic children: none occurred.
Why this is weak evidence about H-015: the guard is inert here, so guard-passing = gate-passing and the
table says "at n <= 7, depth <= 6 the gate is already sufficient", not "the guard is". Depth is small: F-038
says corruption needs seven clean steps at n=10 to come back to a line.

The ALARM path (E-032, F-038; the relation dual of `03033030`, steps `[4,6,4,6,9,4,4,6]`, using the
`relationDual()` start), per step:

| step | 1-6 | 7 (vertex 4) | 8 (vertex 6) |
|---|---|---|---|
| gate | admits | **admits** | admits |
| 2.3(c) tilting | yes | **NO** | (yes, on an algebra already outside the class) |
| Cartan = r C r^T | yes | **NO** | yes |
| Coxeter key held | yes | **moves** | held relative to the wrong start |

So the exact criterion refuses precisely the step the guard refuses, and it does so from the parent
alone, before the rewrite. It also shows the wrong-key step is genuinely a *non-tilting T*, not a bug
in the rewriting that the polynomial happens to expose.

Unfinished: n=8 depth 2 (858 algebras at depth 1 took 51 s; killed when the chair called wrap-up).

## Reproduction

```
timeout 10m .venv/bin/python workshop/rounds/001/scholar_h015.py 6 --depth 6 --unguarded   # 3 min 40 s
timeout 10m .venv/bin/python workshop/rounds/001/scholar_h015.py 7 --depth 4 --unguarded   # 6 min 10 s
RD=1 timeout 10m .venv/bin/python workshop/rounds/001/scholar_h015_f038.py                  # seconds; the ALARM path
```
Output is a cross-tabulation `(guard, tilt, cong)`; anything but `('guard','tilt','cong')` is listed
under "suspicious". Without `RD=1` the f038 script uses the plain opposite and the path does not apply.

## Prior record

H-015 asks for a second invariant along guarded paths; E-034 tested one cospectral pair (n=9, depth 6);
F-016 checked 280 newly allowed steps by polynomial only; F-047 gave SNF profile (Cartan-only, hence
automatically preserved by any genuine tilt, by Ladkani 3.6). Not in `research/`: applying Prop 2.3(c)
(the literature file 1001.4765 recommends replacing `isMutable` with it, "How it lands on our problem", but
no finding records doing so) and the observation that the E-032 step-7 failure is a failure of 2.3(c).
Nothing here contradicts RETRACTIONS; no rediscovery found, but the n=6,7 null is expected from F-038.

## Code changed

New, mine, not in `quivermutation/`: `workshop/rounds/001/scholar_h015.py` (audit) and `workshop/rounds/001/scholar_h015_f038.py`
(ALARM path). No library code touched, no tests run or needed. Caveat: `tiltingPlus` in the script is
the 2.3(c) rank test; it is not yet in the library or under a test.

## Next

- Chair/toolsmith: promote `tiltingPlus` to `procedure.py` as `isTilting(quiver, relations, vertex)`, with a
  test pinning the ALARM step 7 (False) and every vertex of every n=4..6 LNA (True where the gate admits).
  Then the search could use it as the exact gate and keep the Coxeter guard as a cross-check, which is a
  genuine sharpening of F-016 (necessary -> iff).
- Overnight proposal (not run): the same audit on the region where the guard fires: n=9 and n=10 from the
  leftover orbits of H-013, depth 6-8 with `fingerprint.Visited`; expected cost about 2x per level from
  the n=7 depth-4 figure (6 min), so n=10 depth 8 is hours; size with `--plan` (depth 1 only) first.
  The decisive outcome to look for: a step with guard=passes and tilt=False.
- Skeptic: is 2.3(c) as coded right for non-monomial parents (I row-reduce with `idealBasis` and reduce
  the products `p*beta` by pivots of `(i, head(beta))`)? The ALARM reproduction is the only positive test.
