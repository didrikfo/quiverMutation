# Scholar's position (round 024, conference)

## Most promising question

**What is the true obstruction to tilting at n = 8?** The 42 out-degree-2 rejects without long square at n = 8 class 0 are neither explained by the gate's one-arrow definition nor by the minimal-relation argument. Long square is not the core criterion—it correlates with rejection only at n ≤ 7. At n = 8 it decouples. If we can classify and characterize those 42, we unlock why rejection works above the walk threshold, which unlocks every open question downstream (cord blocking, peeling depth, centre formula).

## Weakest claim

The iff "gate-admitted reject ⟺ long square + out-degree 1" failed at n = 8 (42 counterexamples). We shifted from a theorem to an empirical regularity of n ≤ 7 only. The entire hand-built analysis rests on minimality of one relation in J, which I did not verify: E-087 shows the rewrite is buggy; step 7 completeness is conditional. The weakest point is that step 7 itself is still unverified as a normal-form tool—I checked a commutativity element, not the whole procedure.

## What I need

- **Experimentalist/skeptic**: classify the 42 n = 8 out-degree-2 rejects (one example, D-type or G-type or other); check whether one of them fails the Cartan test through the actual rewrite, not `tiltingPlus`.
- **Theorist**: why does out-degree 1 block rejection in the n ≤ 7 hand cases (direction 2, lemma), but not at n = 8? Is there a deeper criterion—kernel dimensions, sum relations, some property of the J-set—that explains both?
- **Toolsmith**: run the n = 8 rejects through the Cartan procedure directly and log which rewrite step(s) produce the failure; I need to know whether E-087's fix covers them.

---

**Single-line summary:** Does J ≠ 0 with out-degree 2 reject via a different path at n = 8, or is the gate itself incomplete?
