# Skeptic, round 024 (conference)

## Most promising question

**Do the 42 n=8 out-degree-2 rejections without long squares fail to be reachable by key-preserving walks?**

If yes, the walk-restricted version of E-102's iff ("reject iff long square on key-preserving walks") survives. The E-105 reachability check (n=7 reduced-relation kind) gates this: if those examples are unreachable by the walk, the walk restriction is the right one and the step-7 criterion holds as an iff at all n, not just n <= 7. If no—they ARE reachable—then the criterion fundamentally breaks at n=8 and we need a different characterization of the step-7 reject.

## Weakest claim the workshop relies on

Scholar's "reject ⟺ long square" holds universally (E-105). Revised: it holds on key-preserving walks and at n <= 7 off walks, but the 42 gate-admitted n=8 class-0 rejects with out-degree ≥ 2 and no square break the iff. The referee scoped it correctly; the claim is now conditional and untested for its condition. Until the walk reachability of those 42 is checked (or their Coxeter class found to be non-LNA), the iff is fragile.

## What I need from another persona

- **Experimentalist:** Run a batch check of the 42 n=8 class-0 rejects: for each, is the Coxeter polynomial reachable by the key-preserving walk from some n=8 class-0 LNA? (Same method as E-105's unreachable examples: compare Coxeter key to the orbits of class-0 at n=8.)
- **Theorist:** If those 42 ARE reachable, what does their existence tell us about step 7? Do they belong to a different x-value class that has x*beta in I for some but not all beta?
