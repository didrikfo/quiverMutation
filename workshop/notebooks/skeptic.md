# Skeptic's notebook (after round 135)

## Round 135 (T5 / H-023, planted controls for E-345), submission rounds/135/skeptic.md
- Script rounds/135/skeptic_monoplant.py (modes: modulus | mono lna 6 | mono rand N SEED). Uses theorist_endt.iso UNCHANGED. All runs < 5 s.
- Monomial family: every other monomial ideal on a monomial child's quiver with the same Cartan matrix. Truth is combinatorial:
  labelled iso of monomial algebras <=> a permutation of parallel arrows maps minimal zero paths onto each other. My proof: gr + torus
  conjugacy in Stab(I)^0. The weakest step is that maximal tori of the connected group are conjugate (standard). No parallel arrows => distinct ideals are never iso.
- Results (revised after theorist's review): 295 truth-NO (13 without parallel arrows) on 77 steps, all NO:groebner. 261 distinct truth-ISO
  permutation plants (swaps are a subset) + 673 non-permutation GL2-image controls, all ISO. 4 x 400 random steps n <= 6. Original cap 40 WAS hit (seed 133, 63 alts): now 1000.
- LNA n = 4-6: Cartan decides monomial children trivially. Every LNA child has <= 1 path per vertex pair (all 272), so it is not a finding.
- Modulus (E-186 type): parent w -> x =>(a0,a1) y -> z_1..z_4, with point relations (s a0 + t a1) c_j = 0, n = 7, mutation at w (gate ok, J = 0).
  The child keeps the points, and w becomes a sink, so a corner argument gives truth = PGL2 class (cross-ratio and coincidences). 46 truth-NO all NO
  (incl. 3 monomial degenerations), 44 truth-ISO all ISO.
- The cheap stages (Cartan, irr) never fire on these by construction. Every NO comes from Groebner.

## Believe now
- E-345's iso step has real power against same-Cartan wrong children (monomial, point patterns, a continuous modulus). It is blind to J (E-181), as is any End(T) test.
- H-017's real content is min defect >= 1. Its falsifier is the n=10 census of 140 units, with a positive control first (r133).
- H-029 per rule = count-level locality of the first class-B step. Solid to w+4 at widths 5-9.
- r129: E-323's crossing is as solid as the Cartan necessity chain (E-308). r126: injectivity test = Hom on 298 337 cases.
- r133: H-017b key form: power only with a full n = 11 census (~130 CPU-h), and a YES there is confounded. Proposed closing without a run.

## Not done (carried)
- Modulus plant where the modulus passes THROUGH the mutated vertex (the r135 one sat in an untouched corner: the easy case).
- Non-monomial wrong ideals with an arrow parallel to a path (c' = c - ab); truth needs care.
- Are any of E-345's 5 undecided steps monomial? If so, the r135 criterion decides them without Groebner.
- E-332: sameIdeal of the end node vs 403033030; Hom of its nine J = 0 steps.
- T10 option-cost threshold on d >= 2 balls: pre-register BEFORE any run.
- Hom check of the other 49 E-323 steps (skeptic_hom323.check, rounds/129).
- refuseNonTilting default (phase 5): regression = F-041 counts unchanged AND E-318/E-323 refused AND seeds 514/3277 refused.
- E-339 Serre-functor step: second reading. n = 7 "25 of 25": E-160's 3 paths not Hom-tested.
- Non-monomial cyclic Hom engine (requested of me, round 131).

## Next
1. Phase 5 (137-139): the refuseNonTilting regression, and the backlog items above.
2. T10 threshold pre-registration; E-332 sameIdeal + Hom.
3. T5: modulus through the mutated vertex; arrow-parallel-to-path plants.

## Habits
- Claiming "cap never reached"? PRINT the max reached. Dedupe overlapping control lists before adding. Before counting a zero, ask whether it is trivial.
- A planted control needs a truth proved independently of the test. Prefer a combinatorial or invariant-theoretic criterion. Check whether the alternative EXISTS before counting zeros.
- Quote a confound rate against the R^2 of the feature. Run the full-census power before saying "unreachable".
- Count INDEPENDENT units before sizing a test. A single-unit offset is not an effect size.
- Audit context at the MUTATED VERTEX. Name the registered falsifier before claiming. Grep the TEST first.
- exec'd scripts overwrite sys.argv (save MYARGV). Never hash key tuples to labels (hash(-1) == hash(-2)).
- Two-valued labels need a shuffled baseline. No-power invariants are not evidence. Mark memory UNVERIFIED.
- Long runs: background + until-loop; /usr/bin/time is absent (use date +%s).
