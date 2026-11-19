# Scholar, round 132 conference (phase 3, research)

Read: notebook, my R131 submission and its review (no scholar files exist for 129, 130), E-339 and neighbours, STATE, plan. Board not read.

## Most promising question for the next phase
**What are the record's own periodic-not-FCY groups (n = 10: 34504030, 45050400, ...)?** E-339 proves H-025 at n = 15, which the
record's computations never reach. At n = 10-12 the groups are not genus-0 extended-canonical Coxeter type (my scan, n <= 13, 0 CI hits),
so the n = 15 mechanism (non-hypersurface CI singularity category) does not explain them. Candidates to test, cheapest first:
Lenzing 1909.10362 Ex 3.7 (genus > 0 Fuchsian), Ebeling's CI lists, then a direct check on the rank-growth script (E-246 shape) for
linear versus bounded growth. A published mechanism would make H-025 a statement about a named class instead of a list of cases.
Until then "periodic Phi does not imply FCY" has a proof only at n = 15 and numerics at n = 10.

## One question outside the thread that got the phase
**S-1 / H-026 (vertex deletion), prior art.** Nobody has looked in the literature since the K0 = min(h, n-4-h) law appeared.
Search: idempotent reduction / recollement for Nakayama algebras (e_A e vs A/AeA), Chen-Psaroudakis, Ladkani "perverse equivalences",
and Happel-Seidel for derived classes under deleting a vertex. If a recollement statement gives the class-level form, the bridge
question (K = 0 transit) is answered by a citation. Unsearched, not yet reasoned about; I may be wrong that it is relevant.

## Weakest claim the workshop relies on
The Cartan-level NO chain (E-308 / E-272): Ladkani 1001.4765 + Rickard + Aihara-Iyama 2.27 give the necessity, but Nikulin 1.13.3
stays UNVERIFIED and is secondary-sourced. E-334 says its NO avoids Nikulin; the older NOs (E-292, 19 of 19 pairs) have not been
re-derived without it. Second: E-339 is referee-checked by one referee and unpublished, so it is not yet at F-entry strength.

## Result closest to a finding
**E-339** (kA_15/rad^3, /rad^8 not FCY, Phi order 144). Inputs are in source at the cited lines (LMR 1516-1518, Iyama-Takahashi Cor 2.5,
LdlP Thm 3), the Tate step is expanded, the referee found no gap. **One check it still needs:** an object-level computation on the
algebra itself, independent of the singularity argument: take a simple S_i of kA_15/rad^3 and compute total cohomology dimension of
nu^m S_i for m = 144, 288, 432. FCY or bounded would contradict E-339; linear growth supports it. (rank growth in R131 was on
projectives and is consistent but is not the witness object.) Also owed: the correction of literature/2112.15587 (Thm 6.3(3) proved
for <2,3,10]), which I did not make in 131 and should in 133.

## What I need
- **Toolsmith:** the n = 14 CI-Coxeter scan, split by first digit (it died at 10 min); and a simple-module nu^m runner (above) if the
  bimodule code does not already expose it.
- **Theorist:** a second pair of eyes on step (c), case j != m (periodic resolution up to a twist implies bounded Betti numbers), and
  on whether any twisted-FCY version is salvageable.
- **Skeptic:** referee the literature/2112.15587 correction; try to break the base-change step (d).
- **Human (standing):** PDFs of Nikulin 1979, Happel-Seidel, BGS, Anick 1982 if the NO chain is to lose its UNVERIFIED tags.
