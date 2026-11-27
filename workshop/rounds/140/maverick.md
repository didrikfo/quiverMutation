# Maverick, round 140 (closing conference, cycle 2)

speculation level: idea unless stated. Read: notebook, r138 submission and review, E-350..E-358, STATE, plan.

**Most promising question.** Why exactly two periodic tau-orbits, and can (p, r) be predicted? E-355/r138: six certificate algebras (n = 10-15) give two orbits with (8,2)/(6,1), (9,1)/(4,1), (7,2)/(6,1), (11,2)/(6,1) twice, (8,1)/(9,1). Prediction target: a weighted-projective-line / Fuchsian model (2112.15587) giving (p, r) from the weights. Test: N_16(8) <2,3,11] (nobody has run it), with both nu codes. A formula that fits 6 and then predicts a 7th would turn a table into a theorem-shaped statement. Reason: it is the only open point that could explain F-058 rather than re-verify it.

**Outside that thread (dormant line).** T5/H-023 cyclic children: the gate's cycle failure is real for random cyclic monomial algebras but unreachable from LNAs (r131). Does any guarded walk from an LNA ever leave LNAs into a cyclic parent? One bounded search over n <= 8 with the cyclic Hom engine (`AlgG`) would settle whether T5's open cyclic half is empty in practice.

**Weakest claim.** "Two families = two orbits" is about uniserials only. Orbits are seen through the stalk test, so non-uniserial members and non-uniserial periodic orbits are invisible, and the nu2 dimension-vector caveat is still open. A coherent wrong nu passes every per-step check (identity plant, both primes); only the two-code agreement defends it, and the two codes share the interval convention (E-347).

**Closest to a finding.** E-356/E-355 pair: F-055 rebuilt independently (19/19 NO, 19/19 YES, equals E-292 on 38/38). Of mine, the GF(2)/GF(3) reproduction in both codes (E-355, 16 runs, 0 disagreements) is supported. One check it still needs: a q_S-independent NO control (an inequivalent pair by construction, not by q_S separation) for F-055; for F-058 a control at m = 132 for the FCY side (only 2 S,P and 3 uniserials of N_14(3) are controlled).

**What I need.** Theorist: the (p, r) formula or a weighted-line model of the two orbits (open since r134). Toolsmith: nu_bimod over non-uniserial objects (indecomposable summands of a projective resolution of an orbit member) so the orbit claim can leave uniserials; N_16(8) run time on both codes. Skeptic: an inequivalent-by-construction pair for F-055.

**Open items, one line each.**
- H-024/F-055 scoped NO test (q_S-separated): settled at n = 9..12, E-356 rebuild agrees.
- H-024 float bound / unsourced necessity step (Nikulin 1.13.3): still open (UNVERIFIED; the 0.1 slack is counts-only).
- H-024 singular S (2 of 6 n = 12 low-lam groups), planted power on a second group: still open.
- H-025/F-058 second prime GF(2), GF(3), both codes: settled.
- H-025 FCY size control at m = 144, 60: settled; at m = 132: still open.
- H-025 module AR-sequence route to the two families: retired (empty at n = 14-15, powerless at n = 10-11).
- H-025 why two orbits: still open (my lead question); non-uniserial orbits and the 3 one-shift groups: still open.
- H-025 two-nu agreement: settled for what it checks; it cannot detect a shared convention error, so not a proof.
- H-032 n = 11 walk is tilting at every J = 0 step (E-357, 0 disagreements / 140): settled; sameIdeal end-node (E-354): settled.
- H-032 default for `refuseNonTilting`: keep off, since it costs 93-180 in-class nodes at E-316 (E-350); decision rests on the chair.
- H-032 depth 8, 4 duals, depth 5+ at E-316: still open (overnight-sized); the generic "guard cannot leave a class through J != 0" question stays open.

Most promising question in one line: why exactly two periodic tau-orbits, and does a weighted-projective-line model predict (p, r), tested on N_16(8) <2,3,11]?
