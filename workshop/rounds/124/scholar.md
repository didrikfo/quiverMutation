# Scholar, round 124 (conference, phase 1, cycle 002)

Read: my notebook, rounds 121-122 submissions and reviews (no scholar file in 123), E-315, E-319, E-320, E-322, STATE, plan. Nothing run, nothing else edited. Board not read.

## Most promising question for phase 2
Does ENS's (T) (letter digraph of d acyclic) hold on the real Construction-2.2 parents, and does the weight (C) of E-322 imply it? If yes, H-031 is "ENS cancellation (cited, k^n transfer mine) + (H3)", and the (LF)/(H1w) proofs of E-319/E-322 become a second, independent route to the same step. If (T) fails where (LF) holds (or the reverse), the two routes differ and that is itself a result. Cheap: feed the 82 LNA instances of `theorist_lnaparent.py` to the `scholar_ens.py` predicates. Open since R122 ("(T) for general parents: guess E-284 weight"); E-322 now supplies the weight.

## One promising question outside the thread
H-025 (Phi periodic => FCY, novelty retry): dormant since the backlog, never searched properly. One reading round: Herschend-Iyama (n-representation-finite), Lenzing-de la Pena, Kussin-Lenzing-Meltzer, Bocian-Skowronski; goal is a prior-art line or a theorem number, before any F-entry. Second choice: H-027 silting-versus-tilting via Aihara-Iyama Thm 2.27 (E-308 left it open).

## Weakest claim the workshop relies on
Mine: every Cartan-level NO leans on necessity (derived eq. => P C P^T = C'), with a published source (BHL 0906.3422 Prop 2.1, Ladkani 1001.4765 Lemma 3.5 + AI Thm 2.27) but Nikulin Cor 1.13.3 (the lattice step behind q_S) is still secondary-only (29 secondaries, no primary). Also BHL Prop 2.2/Phi needs C invertible over Q. Workshop-wide: the J = 0 premise, now reduced to (H3) on real parents plus Hom/End(T) (E-322).

## Closest to a finding this phase
E-322: (LF) proved for Construction-2.2 parents under (C); (H1w), (H3), H(I) = 0 machine-checked on 82 vertex instances. Caveats it states itself: the 82-instance (LF) check is vacuous (acyclic letter quivers), (C) is sufficient not necessary, no loop parents with degree-1 loops. One check it needs: a parent with a genuine degree-1 loop (a real f, not a formal M1-M6 model) that satisfies (C) and passes (H3) in the engine, i.e. one non-vacuous instance, plus whether every gate-admitted J = 0 parent can be reweighted to meet (C) (count on the 600 E-256 parents).

## What I need
- theorist: a dump of the E-322 instances as letters, degrees and d (any format; I adapt `scholar_ens.py`), and one non-vacuous loop parent.
- human/dispatcher: PDFs of Chekanov 2002 (Invent. Math. 150), and the primaries Nikulin 1979, Happel-Seidel, BGS, Anick 1982 (not on arXiv); otherwise they stay UNVERIFIED and phase 5 retires them.
- referee: a verdict on the k^n transfer of ENS (word-level, argued by me, checked for commutativity by the theorist only).

## What the literature changed
H-031 no longer needs completion: KY Lemma 2.10 and DWZ Thm 4.6 do not apply (complete, no-loop), Oppermann has no general lemma, but ENS (proofs of Cor 3.11 and Thm 3.10, Move III) gives a non-completed substitution plus trivial-part step under (T)+(P); BHL gives a second source for H-024's congruence step. R121's "M4 all Theta in C" was wrong (3 Theta stay).

## What phase 2 should not do (published or settled)
- Do not search again for a completed-setting splitting lemma (KY, DWZ, Van den Bergh App. A, Oppermann): read, negative.
- Do not re-prove the H homotopy or its sign (ENS Cor 3.11 proof; sign machine-checked R122, mirrored in the engine).
- Do not try to re-derive necessity of P C P^T = C' from scratch: BHL Prop 2.1 states it for finite-dimensional algebras.
- Do not hunt Nikulin 1.13.3 on arXiv again (null in R122); cite it as secondary or ask for the PDF.
- Do not describe E-309 and ENS remainders as isomorphic: only letter counts match.
