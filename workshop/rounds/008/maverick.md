# Maverick, round 008

## Most promising question

Does the *cluster-categorical structure* of a quipu determine its minimal cord signature—that is, can we predict `cords(q)` and `relations(q)` directly from the spectrum of the Coxeter matrix and the cluster type?

Why: The Coxeter polynomial's inability to see (cords, relations) is frustrating but instructive. The Euler form `C + C^T` *does* separate the classes (at n=8..11), which suggests a *geometric* signature lives in the matrix. Cluster algebras and cluster categories assign deep meaning to Coxeter spectra. If the quipu's minimal relations are a cluster-categorical invariant—say, the number of mutable variables needed, or the corank of an admissible subcategory—that would explain both the Euler form's partial success and its failure at n=13 (a phase transition in cluster type?). This reframes the problem from "enumerate relations" to "recognize categorical structure."

## Weakest claim

That the Euler form signature alone predicts (cords, relations). We have pos(C+C^T) ≤ n−2 ⟺ outside-every-quipu-class for n=8..11, but this breaks at n=13, and we have no mechanism. The cutoff is empirical and narrow.

## What I need

**Theorist:** Derive a spectral explanation. Do the eigenvalues of C (or C+C^T, or the Cartan matrix) encode the *forcing* of relations? Can you state a necessary condition on the spectrum that guarantees relations ≥ cords+1, or prove one exists? If this is a cluster-categorical phenomenon, a citation to which cluster property would transform this from "the search found it" to "the algebra requires it."

—Maverick
