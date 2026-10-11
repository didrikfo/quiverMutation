# Scholar's statement (round 012, conference)

## Most promising question

**Can we find an LNA-reachable gate-admitted non-tilting mutation at n ≤ 9?**

This question bridges the hand-built anomaly (E-080, n = 5) with the field's open cases. If the answer is no, `isTilting` characterizes exactly the gate-admitted tilting mutations reachable by walks; if yes, the hand-built A5 is not special and `isTilting` requires rethinking. Either path narrows H-015 toward a theorem.

## Weakest claim

The workshop relies on the caveat: "hand-built gate-admitted rejections do not count for `isTilting` promotion; only LNA-reachable cases do." This is a reasonable hypothesis (walks are how mutations arise in practice), but it is not proven—it is a policy decision. If the answer to the question above is no, we can upgrade it to a finding; if yes, we must choose between promoting `isTilting` anyway or admitting the boundary is subtler than the gate alone.

## What I need

1. **Toolsmith:** check if any LNA walk at n ≤ 9 reaches an algebra that fails `tiltingPlus` at some vertex (either gate-admitted or gate-refused; either answers the question). Scope: n = 5..9, all placed cores, a shallow search (same reachability as E-076/E-077 if possible).

2. **Theorist:** derive the one-map identity independently (Ladkani 2.3(c) = Aihara-Iyama 2.32(b) = `tiltingPlus`) from first principles, for the small commutative instance at n = 5..7. This closes the gap between the code and the literature and makes CHZ's dependence on it explicit.

3. **Anyone with PDF access:** verify CHZ Cor 3.6, especially the "monomial" hypothesis. ArXiv 403 blocks the route; ask the chair if the paper is cached or has a DOI link.
