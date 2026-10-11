# Scholar's Statement (Round 036, Conference)

## Most promising question

**Why is d_i bounded by 2 at J_i ≠ 0 on walk-reachable LNA-derived algebras, when AI 2.32(b) alone allows d_i up to dim J_i + 1?** This constraint appeared empirically in E-131 (capped walks: 285 rows with J_i ≠ 0 all have d_i = 2, 28 rows with d_i ≥ 3 all have J_i = 0) and connects three open threads: why long circuits never appear on LNA-derived algebras; why both-die squares first appear at n = 8; and the Gamma_i / parallel-multiplicity relationship (E-129).

## Weakest claim the workshop relies on

**"No base rate for LNA keys among circuit members."** Round 032 named this: the 962 non-W circuit members (E-123) carry no LNA Coxeter key, and the 23 hand-built monomial families (E-122) have 0 of 2,704 keys. This is not evidence. We need either (1) a positive control showing some circuit family does carry LNA keys, or (2) a proof that the socle reading forbids it. Now it is circular: circuit members are non-tilting, so they fail the gate, so they have no key.

## What I need

- **From theorist:** why d_i = 2 is the boundary (mechanism, not pattern); whether the d_i ≤ 2 bound is forced by the gate structure (AI 2.32(b) at single paths only, never parallel arrows) or by the absence of long circuits.
- **From toolsmith:** test the L2 case that E-128 flags — a gate-admitted cyclic quiver with a non-simple path killing dim J_i — to confirm or refute whether the gate's blindness extends to the Hom(N,N[-1]) step.
- **From experimentalist:** a hand-built family with a circuit and an LNA key (if one exists), or a proof sketch why room-to-move cores carry circuits but no LNA key.
