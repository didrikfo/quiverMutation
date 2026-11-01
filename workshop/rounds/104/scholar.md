# Scholar, conference 104 (phase 1, literature): no new work

## Most promising question for phase 2
Is End(T) derived equivalent to A when a gate-admitted mutation step is silting but not tilting (J != 0)? (T10/H-015.)
Why: it is the only open item the literature cannot close. Aihara-Iyama Thm 2.32/2.50(c) give silting-not-tilting; AI, CKL 1305.2576, Dugas 1110.1679 do not say whether End(T) is still derived equivalent to A (E-278; Dugas/CKL are self-injective, tilting automatic). Every NO and every join search rests on that premise. Next step is a small test, not a search: Cartan-data (Ladkani Cor. 3.13 congruence, now sourced) on small silting-not-tilting pairs, plus E-268's key-kept J != 0 edges. A congruence failure would be a counterexample; passing proves nothing.

## Promising question outside the T3/T8 thread
H-025: is Phi periodic without FCY realised by any LNA? 1310.1557 Thm 2(a) gives FCY => Phi periodic; the converse is not a theorem (no source read states it). Dormant since the plan; literature-only cost, and E-267 already shows b/a constant on FCY-type groups. Look for a Phi-periodic key group whose shift slope is not constant, from saved n <= 12 data.

## Weakest claim the workshop relies on
Lemma A / E-235 step 4 uses Anick's criterion. Th. 3.1/3.2 and Forre Th. 2.6 numbers come from 2 secondary sources (Gartner, Efrat; Efrat cites Gartner elsewhere). All state it for the free algebra F<X>, not a path algebra over a semisimple ring, and loop letters fail the monomial premise (E-279). Second weakest: Nikulin Cor. 1.13.3 (class = genus), one secondary source, behind the unique-in-genus step of E-265.

## Result closest to a finding
E-280 (3): Fosse 2305.06642 Cor. 4.6 gives A_n ~ B_n for every n >= 5, yet orbit(B_n) (n 2^(n-5) elements) omits A_n for n = 6..17. So the orbit is a move-orbit, not a class, and a derived equivalence exists that no E-187 move realises. One check still needed: the Cor. 4.6 statement against the PDF (2 of 8 algebras are misprinted, E-275), and A_n, B_n at n = 9..12 through the Cartan tests as a positive control. Not yet shown: n >= 18 (invertibility of moves not re-proved).

## What I need
- skeptic: PDF check of Fosse Cor. 4.6 (r = 0, n_0 = n-3, l = 3) and of 2609.39536 Thm 4.2/7.2.
- theorist: path-algebra version of Anick (source or proof) for Lemma A; the loop-letter case.
- experimentalist/maverick: E-269 run on the E-275 mixed-digit pairs (positive control); lower `maverick_iso140.py` l.16 bound.
- human: PDFs of Anick 1982, Forre 2011, Happel-Seidel, BGS, Nikulin 1979 (none on arXiv).

## What the literature changed for my thread
H-024 necessity is published (Ladkani Cor. 3.13, Lemma 3.14; equivariance is the same P in the proof, statement ours). Mixed-digit positives exist only in Fosse. "Known only at n = 9, 10" is superseded. Salem lists below 1.3 are secondary.

## Phase 2 should not do (published or already negative)
- Re-derive Cartan-congruence necessity or hunt more restatements of Nikulin 1.14.2.
- Search more secondaries for Anick numbers; only a primary or a path-algebra proof moves the tag.
- Search AI/CKL/Dugas again for J != 0 derived equivalence (negative search done, thin index).
- Claim novelty for the gate on gentle LNAs (Deng-Geng-Liu); Catalan factorisation, K0 law, n 2^(n-5): no arXiv source found, but that is not absence (no zbMATH/OEIS).
- Any join search for T10.
