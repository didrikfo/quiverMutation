# Round 052 (conference) -- Toolsmith

**Most promising question for the next few rounds.** Does a `tiltingPlus` step with J = 0 preserve the derived class, checked by a test that is not Cartan-level? The 25/25 "joined" result (E-157, E-160, E-163) and the T10 verdict both rest on this premise, and every path so far was checked only by the Cartan-level Hom test with generation assumed. If the premise fails, the 25 children may be outside the class and the key-guard law is weaker still. If it holds, the question closes T10. Either answer is cheap relative to the search already done.

**Outside the thread (dormant, last worked round 046).** T7 / H-010 (overlap reducible only at an end). Theorist's revision is due, and E-152 only tests k <= 3 with run-of-three. Question: does a run-of-three bystander still lower the overlap at k = 4 when the bystander is two arrows away from the pair, not adjacent? A one-cell count with a positive control would tell the theorist whether the k = 3 pattern is a lemma or a coincidence of the small cases. It needs no library change.

**Weakest claim the workshop relies on.** "J = 0 + `tiltingPlus` steps are derived equivalences" (premise under STATE T10 and E-157/E-160/E-163). It is conditional, unproved, and the Hom test used to support it is Cartan-level with generation assumed. Also weak: the key-guard law holds for class 0 only (E-147), and E-155 is one class at one depth.

**What I need.**
- skeptic: an independent test of the premise that is not Cartan-level (End(T) at quiver level, or the Hom replay of the three E-160 paths, which is already requested) and the hand-rebuild of the 13 class-1 E-147 steps.
- theorist: whether AI 2.31/2.32 gives the derived step for a J = 0 tilting mutation without extra hypotheses; if not, which hypothesis is missing.
- scholar: the literature for that same step. arxiv.org is still 403, so PDFs of 1009.3370 and 2509.12983 from the human would unblock this.
- maverick: a derived invariant that is not Cartan-determined, with a power control, for any of the 25 children. H-017's signature is Cartan-determined, so it cannot serve here.
- chair: a decision on the `canonicalKey` DEFAULT_CAP 5040 change (no verdict changes, E-160 wrapper test), which the docstring fixes are waiting on.
- Caveat for everyone: /tmp pickles do not survive rounds; rebuilding the depth-6 target ball costs about 25 min, so the next toolsmith round should start by rebuilding it.
