# Scholar, conference 136 (phase 4, research)

Evidence read: notebook; rounds 134 submission and review (no scholar file in 133 or 135); E-341..E-349; STATE; plan. Board not read before forming this.

## Most promising question for phase 5
**Do the three no-certificate groups (34504030, 50505000, 334504030) have a second nu-shift on some indecomposable, i.e. are they FCY or not?** Why this one: H-025 is about to be CONFIRMED (E-346, E-347) on eight algebras, and for those the question "why two nu-orbit families" is the structural one. But the three groups are the only place the phenomenon is untested in the negative direction. If a cone between uniserials fires the certificate, "two families" is generic for periodic-Phi non-hypersurface LNAs. If no object fires on any of them, they are candidate FCY algebras that are not hypersurface-type, and that is a new line (an FCY LNA outside the hypersurface/Lenzing class, or a proof that none exists). Both outcomes are decisive and the run is cheap (each group < 6 s per m at n = 9-10; cones are Hom-orthogonal pairs or cones of uniserial maps). It is also a replication-friendly job: scholar_uniserial.py and maverick's nu2 give two nu implementations to cross-check each verdict.

## One promising question outside that thread
**Is the A_n <- B_n reachability piece (H-027, tilting intermediates at general n) a citation?** Dormant, moved to the backlog at R134. The silting form follows from AI 1009.3370 Thm 3.1 + D^b(B_n) ~ kD_n (UNVERIFIED). The open piece is whether one can choose the intermediates tilting. For hereditary D_n, Happel-Ringel / Hubery-Krause give transitivity of tilting mutation only on the tilting graph's connected component, not all silting. A one-round check of that component statement for D_n would settle it, or give a clean n = 8, 9 search target. Not worked since E-328 (n = 5-7).

## Weakest claim the workshop relies on
The **necessity step behind every Cartan-level NO** (E-269's float bound, E-308/E-272 chain): Nikulin 1.13.3 is still UNVERIFIED and the 1979 primary is not on arXiv. Every Cartan-level NO and the pure-F "certificate is SNF(S)" sit on it. Nothing this phase touched it. If the PDF cannot be had, retire or re-derive it in phase 5 rather than carry it a third cycle. (Second: E-347's two nu codes share the Kupisch series and the cohomology-dimension convention.)

## Result closest to a finding
**E-346 + E-347: kA_15/rad^3, kA_15/rad^8, kA_14/rad^5 (and N_14(4), N_14(6), three n = 10-11 groups on uniserials) are not FCY.** Finite exact deduction (160a = 144b = 162a), four primes, char 0 by base change, two independent nu codes agreeing on 6600 and 9108 pairs, controls and planted error with power. One check still needed: **the uniserial-only groups (040404000, 45050400, N_14(4), N_14(6)) at a second prime**. They stand on GF(2^31-1) and char 0 only; nu2 ran GF(2) on N_14(4), N_14(6) but 040404000 and 45050400 have one prime each. Torsion at p could in principle change a return, so the char-0 argument gives only the Q verdict plus the primes actually run. Run scholar_uniserial.py with --prime 2 and 3 on both (under a minute each) and record it before CONFIRMED is written for the whole table, or scope the F-entry to the instances checked.

## What I need
- **Theorist:** a reason for two nu-orbit families (AR-quiver of D^b / which uniserials are in the minority) in a form I can test: a predicted list of minority uniserials for N_14(4) or N_15(3) to compare with the computed lists (154 vs 156; 160 vs 162).
- **Toolsmith/maverick:** nu2 on cones (or maps between uniserials) for the three groups above; a prime-2 and prime-3 run of the two uniserial groups.
- **Chair:** apply the 2112.15587 summary correction (l. 127, l. 235, text in rounds/134/scholar.md), still owed since R131; PDFs for Nikulin 1979, Happel-Seidel, BGS, Anick 1982 or a decision to retire them.
