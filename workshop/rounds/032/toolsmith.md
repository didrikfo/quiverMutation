# Toolsmith position statement, Round 032

**Most promising question for next few rounds:**
Can the two-term kernel structure J_i (specifically, why Gamma_i has at most 2 edges on LNA-derived algebras) explain why long circuits are impossible, and therefore why W (length-2 ground path) is necessary but insufficient?

Why: This question ties together all the weak points. The converse of W rests on limited data (70 same-shape rejects, no "cancels" or parallel-arrow test). But the deeper issue is that W is only a symptom: the doubled/tripled-arrow controls are gate-admitted but fail the real test (Cartan). If J_i structure forces acyclicity, then we have a mechanism that explains why W works on walks (all J != 0 are length-2 paths) without needing to validate W exhaustively. The resolve-the-cospectral-LNAs task becomes a necessary base-rate check on this mechanism.

**Weakest claim the workshop relies on:**
The code's W is sound. E-121 shows it is True on "cancels" (where J = 0, so it cannot distinguish nn from nz) and misses the tripled-arrow H chain (E-112). On hand-built controls, gate admission is satisfied but Cartan fails. We have no positive control for parallel-arrow rejects, and the 61 loose D' rejects may be gate-admitted for a different reason than the 70 same-shape ones.

**What I need from another persona:**
- From **theorist**: a proof or invariant showing why Gamma_i components have at most 2 edges on LNA-derived algebras (not just on the n = 8 c0/c1 sample). This is the crux: if it holds, W becomes a corollary, not a primary claim.
- From **experimentalist**: walk data showing whether the hand-built "cancels" and parallel-arrow controls are gate-admitted (and if so, why); and dim J_i on all walk rows at n = 9/10 to test the 2-edge claim in data.
