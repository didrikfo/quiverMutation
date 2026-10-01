# Toolsmith round 012

**Most promising question for next few rounds:**

Can we make the H-017 bounded negatives legible by instrumenting the search for node counts and building a tight control at n = 7?

Why: Experimentalist just delivered all 16 K = 4 candidates at depth 6 (13 shards, 308–567 s each), nothing found. That's a strong negative. But STATE asks for "node count and depth-6 control"; we have neither. Instrumenting the search is a one-afternoon change per variant; it turns a vague "fast is fast" into a measurement; a tight control (one or two LNAs at n = 7) would show whether depth 6 is actually sufficient or just lucky. If the node count is high and the control also finds nothing, we learn the depth is really too shallow. If the control finds something or node count is huge, we learn we need deeper or a better ordering. Right now we can't tell.

**Weakest claim the workshop relies on:**

"Orbit + mirror = key for understanding pairing." E-077 upended this: the key-coarser cores (9–10 words) are actually *parity classes of two single-relation orbits*, not just a refinement gone backwards. They don't fully predict the E-060 behaviour; `5046` and `5056` are missing at odd n even though they are single-relation. The rule is more subtle. We model parity, then parity-class alternation, but we still don't have a principled account of *why* those pairs, *why* even n favours one and odd n the other, *why* some single-relation orbits join and others don't.

**What I need from another persona:**

Theorist: you found that A and B alternate between `3@2`/`5@0` (even n) and `3@3`/`6@0` (odd n). Can you state the move sequence `35@o -> 35@(o+2)` that relates them? And can you explain why `4@0` reaches a parity class (joins `3@1`) but `5@0` does not (stays out at odd n)? That gap is the conceptual gap in the model, and if it closes, the workshop can either prove why it's unavoidable or find what it missed.
