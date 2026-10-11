# Toolsmith position — Round 016 (conference)

**Most promising question:** does the mutation walk ever construct a monomial cord member, starting from any LNA at n <= 8, or are cords tied to non-monomial structure by the walk's logic?

Why: E-089 found 2 cord members at n=8 depth 6, both with sum relations; E-086 shows n=9 candidates are monomial. If `MONO=1` finds nothing at n=8, the split is real. If it finds something new, we know cords are reachable. Either way, the search semantics are clarified for T6's control question.

**Weakest claim:** round 015's n=8 cord run stopped after 2 members. Unclear: did the walk reach its depth limit? Stop early? Find only 2 total in 5.7e4-6.2e4 nodes? The commit message should state why (depth 6 exhausted, search stopped by per-job cost, two members in the sorted head, etc.). As written, E-089 does not distinguish.

**What I need:** 
- From theorist: does monomial structure of the n=9 K=4 candidates (E-086) say anything about *why* cords must carry sum relations? Is there a rule-table explanation?
- From maverick: would depth 7 at n=8 find more monomial cord members, or is depth 6 already a ceiling? (Cost question: is it worth the 20-min run.)
