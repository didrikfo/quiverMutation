# Scholar's notebook (after round 006)

## Believe
- Ladkani 1001.4765 Prop 2.3(c) = Aihara-Iyama 1009.3370 Thm 2.32(b) = `tiltingPlus`: one linear map
  (p |-> (p beta)_beta over arrows beta out of k). Two papers, one test; not independent of the code.
- E-032 ALARM step 7 is a correct rejection: parent has c = [8,6,4]+[8,10,4] nonzero, c*(4>9)=0
  (commutativity relation through the vertex). Silting not tilting; key moves.
- The repo performs the right mutation; the mirror test disagrees with the key on steps 1,3,6.
- CHZ 2509.12983 Cor 3.6 (path-wise) as summarised passes step 7; Prop 3.5 (socle) fails it. So Cor 3.6
  is monomial-only. UNVERIFIED against the PDF: arxiv.org is blocked by the proxy (WebFetch EGRESS_BLOCKED).
- Literature does not bear on H-015 (guard sufficiency) beyond per-step exactness.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same (script sense). Say which.
- Earlier: starts n=5..7 gate <=> 2.3(c); no second gate-admitted rejection at n<=7 (E-055, E-057).

## Did
- R001: scholar_h015*.py. R002: nonmono control. R006: rounds/006/scholar_step7.py (witness),
  scholar_sides.py (left vs right test along ALARM path). Seconds each.

## Next
- Read the CHZ PDF (need network) for the hypothesis on I in Prop 3.5/Cor 3.6; then fix literature note.
- Independent check: build End(mu^-_{P_4}(A)) from the cone and compare Cartan with repo child.
- Overnight audit n=9/10 (Menu 4); recommend `isTilting` promotion only after a second gate-admitted rejection.
- Lessons: a paper's iff can be the same computation as the code's test; say so. Check hypotheses (monomial).
