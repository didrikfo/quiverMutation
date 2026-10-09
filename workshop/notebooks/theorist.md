# Theorist notebook (rewritten round 049)

## What I believe now
- H-010 unproved, SUPPORTED; no counterexample, no invariant (see 046 below). L1: one-step LNA-to-LNA mutations change relation starts only at v-2..v, exactly 2 per LNA (E-151), unexplained.
- H-020 is a statement about the move set, not the class. Hypotheses H1-H6 are in rounds/049/theorist.md. H1 (floating rules are length independent) survives: 0 failures in 11 170 applications at lengths w+5..11 for w <= 5; w = 6..8 (most of the 414 rules) unchecked beyond w+4.
- F-051's "interior is one orbit, one verdict" is wrong after F-053: `45` at n = 13 has orbit sizes by offset 0..6 = 2386/1127/4217/4217/1127/2386/447 (pairs o <-> n-8-o; F-053 already says this, sizes are new); the equal verdicts across pairs are data, not a consequence of translation invariance (H4). Weakest of the hypotheses (H1 tested for w <= 5 = 30/414 rules only): H3 (outside = incomplete move set; refuted once, F-052).
- "Outside the class" for the E-149/E-152 children: inside iff a tilting complex gives an LNA of the key class; witness = J=0 tilting path to an LNA (proof of insideness); outside needs a derived invariant with a power control, which does not exist at n = 7 yet. Children (class 2, 9): 7 vertices, 7 arrows, one commutativity relation + zero relations, no parallel arrows. One-step reverse: same-direction step gate-closed 9/9; opposite step gate-open, J = 0, but lands on neither parent nor its opposite 9/9.
- Earlier, standing (042): C_B = C' + E_vi; Q_2 = 1 iff c_2 = 0 (proved when F e_i = -e_m, 296/296 s = 1 samples); s = -2 observed only.
- Bystander (045/046): lowering of max overlap only with a run of three; margin 2, one bystander, k <= 3.

## What I tried
- 049: `rounds/049/theorist_{rulelen,orbits45,children,reverse}.py`; pickle for c2 rebuilt in the scratchpad (405 s via rounds/046/skeptic_collect.py 7 2 20000).
- 046: theorist_{single,leftsweep,locality,whichv}.py. 045: theorist_bystander.py. 042: theorist_*.py.

## Next
- Follow the 9 (then 16) tilting neighbours of the children with skeptic_back-style J=0 search, after a power control exists.
- H1 for w = 6..8 at length 12-13 (overnight); H6 ablation (no anchored rules) at n = 13, 14 on `45`.
- H4: a map between orbits of c@o and c@(o+1)? why sizes differ but verdicts agree.
- Old: prove L1 and the 2-per-LNA count; k-step window over non-LNA intermediates; two bystanders; the 042 programme (c_2 = 0, s = -2 dual).
- Blind spots: class 2 only for the children; 'outside' evidence here is conceptual, no invariant computed; H6 and H3 untested.
