# Scholar's notebook (after round 011)

## Believe
- Ladkani 1001.4765 Prop 2.3(c) = Aihara-Iyama 1009.3370 Thm 2.32(b) = `tiltingPlus`: one linear map
  (p |-> (p beta)_beta over arrows beta out of k). Authors' statement, no derivation; not independent of the code.
- NEW r011: gate-admitted, non-tilting mutations exist at n = 5: vertices a,b,c,d,e, square a>b>d, a>c>d,
  d>e, relation abde = acde (length 4). gate True, tiltingPlus False, Cartan congruence False at d.
  Padded versions at n = 6, 7 the same (6/6). True square (abd = acd) is fine; monomial control gate-refused.
  So E-066's "n = 10 first size" is false as a lower bound for the shape. Not shown reachable from an LNA.
- E-032 step 7 is the same shape (commutativity relation through a vertex with one outgoing arrow).
- CHZ Cor 3.6 path-wise vs Prop 3.5 (socle): reasoned that Cor 3.6 needs "monomial"; STILL UNVERIFIED,
  arxiv.org returns CONNECT 403 through the proxy (r006 and r011). Do not retry; no other route found.
- Terms: gate = `mutationIsPossibleAtVertex`; guard = Coxeter key same. Say which.
- Earlier: LNA-started walks n=5..7: gate <=> 2.3(c) (E-055, E-057).

## Did
- R001/R002/R006: scholar_h015*.py, nonmono control, scholar_step7.py, scholar_sides.py.
- R011: rounds/011/scholar_square.py (18 hand-built algebras, seconds). Arxiv fetch failed (403).

## Next
- Ask toolsmith: reachability of A5-type algebras from LNAs by gate-admitted steps at n <= 9.
- Chair decision on `isTilting` promotion (hand-built counts?); a unit test with A5 is cheap.
- Someone with the PDF: Cor 3.6 hypothesis; then fix the flag in the literature note.
- Independent check: End of the two-term complex for A5 at d, compare Cartan with repo child.
- Lessons: a "first size" claim in a Limits paragraph is cheap to test by direct construction;
  hand-built algebras reach cases LNA-rooted walks cannot.
