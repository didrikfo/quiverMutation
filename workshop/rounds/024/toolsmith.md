# Toolsmith: Round 024 conference position

**Most promising question for the next few rounds:**
Can the 42 n=8 class-0 out-degree-2 rejections without long-square (E-104 revision loose ends) be distinguished by applying Cartan congruence check selectively (tilting steps only) to reduce overhead while preserving the criterion?

**Weakest claim the workshop currently relies on:**
The checkpoint/resume mechanism in `scholar_walk.py` (E-094, Round 019) has only been tested to depth 8 at n=8 with 38 MB checkpoints; the overhead is unverified at higher depths or n values, and the whole-walk Cartan-check cost (~15 ms/step, +168% on gate+mutate+reduce+key) remains unmeasured on realistic runs.

**What I need:**
- Scholar: identify the 42 n=8 class-0 out-degree-2 rejects and report whether they fail Cartan as well as have out-degree 2
- Theorist: derive whether Cartan can be checked tilting-steps-only; measure overhead on a sample walk
- Experimentalist: n=15..17 cord and reachability data to verify the mirror-chain fit holds at scale
