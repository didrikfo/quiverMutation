# Conference 116 - experimentalist (no new runs)

**Most promising question.** Is a key class in S-1 a single move orbit (or derived class)? H-026's K0 = min(h, n-4-h) now has 10 fits, 3 pre-registered and held, but every number is a count of Coxeter-key classes of the image after deleting a vertex. "Class" has never been checked to be an orbit. That is the F-candidate blocker, and it is testable: at n = 13 and 15 (class sizes 504-4597 rows) take the 2 image keys at K = K0 and run the move-orbit / discover step on a sample from each. Either they are different orbits (the split is real structure) or they join (the split is a key-label artefact and the law is about the guard, not the algebras).

**Outside the main thread (T10 got most of 113-115).** T7: a gap-2 shape that CAN lower overlap at k >= 6. Today T7 has no control (E-288: no pair at k <= 4 can lower), so its "stays at overlap 2" is unfalsifiable. Finding one shape that lowers, then asking whether the E-230/241/260 cells do, is a cheap job with `--plan` first; it has been dormant since 107.

**Weakest claim.** The K-threshold law (S-1). It is a fit to key-level counts at head ends, K >= K1 only, n <= 21, with no derivation and no K < K0 data at n = 19 (K1 = K0 there, so "no split above" is all we see). Runner-up: "25 of 25 joined" rests on the J = 0 premise and assumed generation.

**Closest to a finding.** E-304/E-301: K0 = min(h, n-4-h) on 10 lone-3 classes, 3 predictions made before running (5, 6, 3) and all held. Not a finding yet because: (i) key class not shown to be an orbit (above); (ii) tail ends are only tested by a mirror equivariance; (iii) non-lone classes (1, 2, 3 keys, none) show it is not a general law. One check it needs: the orbit check above at n = 15 h = 3 (K = 3), plus ideally n = 19 h = 5 at K1 = 3 (~50 min overnight) to see K < K0 behave.

**Need.**
- toolsmith: a cheap orbit-membership call on rows (move-orbit canonical form, or confirmation that `canonicalKey` is not orbit-complete) so I can test key class vs orbit without writing my own.
- theorist: a reason why min(h, n-4-h): is it the distance of the 3 to the nearer end of the tail after deletion? A derivation would let me predict the non-lone depths instead of tabulating them.
- skeptic: name the planted error you would accept as a power test for the threshold itself (none exists now).
