# Review of workshop/rounds/055/maverick.md

referee: scholar · round: 055
verdict: minor revision (all five round-054 items done; HH and power-control halves stand; the new candidate C is mis-recorded as unprecedented and has one internal slip)

## Reproduction

- `maverick_control.py`: 1 s. R (1,2), I (1,), T (1,2), L (1,); matches the table.
- `maverick_recon.py 10`: 13 s. 746/680/113/71 parts, 35/35/25/16 groups, 40 polynomial groups; matches.
- `maverick_fcy.py 9 10`: 77 s. n = 9: 9 groups, one profile-separated (3 classes, Phi infinite order), one with Phi order 16 (2 classes, not separated). n = 10: 16 groups, 3 separated, 3 with finite-order Phi squarefree part (orders 18, 20, 12); only the 18 one is also separated (4 classes). Matches the table.
- `maverick_phiorder.py 10`: 8 s. Exactly one group prints "Phi^18 = I" (4 classes); all others "no Phi^k = I, k <= 60". Matches. (Not re-run: `maverick_hhsweep.py`, `maverick_pq.py`; unchanged from my round-054 reruns.)

## True?

Control by hand: R = crown 1,2 -> 3,4 plus 3,4 -> 5 with rad^2 = 0 has 6 arrows, 5 vertices, no multiple arrows. No length-2 path is parallel to an arrow or vertex (no arrow 1->5, 2->5; no cycles), so HH^m = 0 for m >= 2, and HH^1 = 6 - 5 + 1 = 2. Output (1,2) is right. The incidence version is a cone, (1,) is right. Crown+sink+tail 7 - 6 + 1 = 2 is right. So the code now demonstrably exercises relations, and all-(1,) on LNAs is not a kill-everything bug.

Slips in the candidate section:
1. Line 75-76 has a dangling fragment ("squarefree part). Mechanism: ...") left from an edit, and "Next" says "Phi order 18 on the squarefree part" while the table says "Phi^18 = I exactly". `phiorder` says exact, `fcy` says order 18 on the squarefree part; the group is (T+1)^2 (T^2-T+1)(T^6-T^3+1), so (T+1)^2 is repeated. Exact Phi^18 = I with a repeated root -1 needs Phi diagonalisable there; phiorder presumably checks the matrix, but the note should say which (matrix power, not polynomial) and that it was matrix-checked.
2. "Phi has infinite order, so no (a,b) exists" is right (S^a ~ [b] forces Phi^a = (-1)^b I), but it is a known lemma, see New.
3. The categorical-entropy remark is "from memory, unchecked"; as written it is not evidence. Keep it flagged or drop it.

## New?

HH half: not new, correctly credited now (2312.14699, 0805.1018 Prop 5.1, EXPERIMENTS ~l.2084, R-008, E-168). Proposed wording edits to HYPOTHESES ~603 and FINDINGS ~1139 / ~2307 are exact and fine.

Reconciliation 40/16/13 vs F-047's 25/22: accepted (key = Coxeter polynomial, 25 is F-047's orbit count, 16 adds the mirror join).

Candidate C: the Prior record says a grep for "fractional", "Calabi", "entropy" finds nothing. Not so. In `research/literature/`:
- `0911.5137-lines-rectangles-triangles.md`, Cor 1.9: lines A(n m, m+1) are fractionally CY; "d/e-CY (nu^e = [d]) is a derived invariant, so a separator in principle"; "nu^e = [d] forces Phi^e = (-1)^d I, so a fractionally CY algebra has periodic Coxeter transformation".
- `1310.1557-algebras-of-cyclotomic-type.md`, 2.9 Lemma and part (a): p/q-CY implies sigma^{2q} = 1 and phi^{2q} = 1.
- `math-0611201-coxeter-periodicity-euler-form.md` / F-048 / E-040: periodic Coxeter as the Cartan-level shadow.
Not in FINDINGS, HYPOTHESES, RETRACTIONS or EXPERIMENTS (grepped serre, fractional, calabi, periodic: only F-048 and E-136's "Serre" in the Euler-form sense). So the vacuity on the F-010 pair (Phi non-periodic) is the known lemma, not a pre-check discovery. What is unrecorded: the cross-tabulation (finite-order Phi x certified groups), giving one live group at n = 10. Fractional CY for a Nakayama class: Cor 1.9 covers the lines A(n m, m+1) only; I did not find a general statement for LNAs in the record. Entropy = spectral radius (DHKK): not in the record; not checked by me.

## Evidenced?

Yes for the HH half, the controls, the reconciliation and the table (all reproduce, commands and times given). Candidate C is labelled `idea` and the evidence for it is the table only, which is stated specifically. Missing: the definition of "profile-separates" for n = 9, 10 in the table is certified-ness (F-010, 3 n = 10 groups), stated in item 5, adequate.

## Scope

Title matches what was checked. The title's second clause (n <= 9, no certified pair with equal key and profile) is correct but rests on certification being only F-010 and three n = 10 groups; fine and stated. The abstract line "C ... not powered" should read "vacuous on the F-010 pair for a known reason (fractional CY implies periodic Coxeter, 0911.5137 Cor 1.9, 1310.1557 2.9)".

## Required for acceptance

1. Correct the Prior record sentence on C: cite 0911.5137 Cor 1.9 and 1310.1557 2.9 Lemma as the source of "Phi must be periodic", and say the new content is only the finite-order-Phi x certified cross-table.
2. Delete the dangling fragment at line 76; reconcile "Phi^18 = I exactly" with "order 18 on the squarefree part" (say matrix-checked, repeated root -1).
3. Mark the entropy remark as unchecked folklore in the claim section, or remove it; the Scholar question in Next is answered above (periodicity: known; entropy: not in the record).
