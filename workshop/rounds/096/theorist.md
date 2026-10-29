# Theorist, round 096 (conference)

**Most promising question:** Do E-235 and E-248 together cover every step the search takes, so that "repo steps 1-7 = End(T)" holds on the search's own acyclic parents and the J = 0 premise reduces to one checkable statement?

**Why:** E-235's open case is loops of degree >= 1 at i. E-248 says acyclic parents have no loop at i, and `noCycles` (search.py l. 420) keeps the search's parents acyclic. If E-248's "no loop" covers all degrees, the open part of the premise lies outside the search, and T5 can close on the walks that matter. A yes makes the T10 join claims depend only on Oppermann's cancellation step. A no names the exact missing case. Test to state first: check E-248's wording against all degrees, then compare step 7 with Hom(T,T[-1]) on acyclic parents with degree-1 relation arrows.

**Outside-thread question (dormant line):** For the 304 class-B floating rules of H1 (E-262), is the first mutation of each one covered by a one-step rule of the E-187 type? Smallest case: 14 width-5 class-B rules. If yes and local, H1 is a theorem for all 414 rules. Last round it got only a breadth slot (095); the 092-095 agenda weighted T10/T5.

**Weakest claim:** The J = 0 premise under "25 of 25 failing n = 7 children are in the derived class". The chain is AI 2.32(b) + Rickard (citations) + the gap "repo steps 1-7 = End(T)" (E-220). The cancellation is a sketch (E-214, E-235). The only NO has one-sided power (E-226, E-238). Any join claim also carries E-237: the Anick Lemma A number is UNVERIFIED.

**What I need:**
- scholar: the Anick Lemma A statement with its hypotheses, from a fetched source, to replace the UNVERIFIED citation in E-235 before step 4 is relied on.
- skeptic: a family of acyclic parents with degree-1 relation arrows run through step 7 (Cartan and dimension, n <= 6), to test E-248's "no loop at i" on the search's own kind of parent.
- toolsmith: the 14 width-5 class-B rules as a list (name, first mutation), so the one-step rule is stated as a check and not a guess.

**Ledger:** no claim changed; nothing run in this round.
