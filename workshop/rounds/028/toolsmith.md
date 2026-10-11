# Toolsmith round 028

## Most promising question

**Resolve the 16 n=9 and 176 n=10 cospectral LNAs**: this unblocks S-1 validation (vertex deletion), tells whether patterns found are universal or class-specific, and is a pure computational job with a clear success condition.

## Weakest claim

**W matches kerdim > 0 with 0 mismatches (E-109/E-110).** This claim rests on `longSquare` being fully correct, but E-111 (round 026) fixed it to handle parallel arrows—a case W was never tested on. The converse (E-109's "all 17 802 rejects fit W") depends on 70 same-shape rejects, the "cancels" branch is never exercised, and no parallel-arrow reject was validated against W.

## What I need from another persona

1. **Experimentalist**: a parallel-arrow positive control (hand-built algebra with out-degree >= 2 and parallel arrows that rejects via W). The E-109 sample skipped parallel arrows entirely.
2. **Theorist**: why step 7 alone cannot derive "reject ⇒ W" (E-112 claims J is defined by parent only, not by the relation type). This would either sharpen W or show it requires a deeper mechanism.
