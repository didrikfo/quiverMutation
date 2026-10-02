# Review of workshop/rounds/023/scholar.md

referee: experimentalist · round: 023
verdict: minor revision

## Reproduction

- `--hand`: 1 s. All seven rows A-G match the table exactly (gate, out(v), dim J, tiltingPlus, longsq).
- `5 --class 0 --budget-sec 300`: 24 s, 6 240 algebras, only J=0 rows (13 320 out-degree 1, 3 300 out-degree 2), no J != 0. Matches the n = 5 row (16 620 steps).
- `6 --class 0 --budget-sec 300`: 300 s, 46 004 algebras (note says 46 858; cap-dependent). J != 0: 1 112 steps, all out-degree 1 and longsq (note: 1 139). J = 0 with out-degree 2: 32 699. Same shape, counts differ by cap/timing as the note says.
- n = 7 not re-run. I ran one case further: n = 8, `--class 0 --budget-sec 200` (8 988 algebras).

## True?

The one case further falsifies the walk statement. n = 8 output:

```
('J!=0','tiltingPlus',False,'outdeg',1,'longsq',True ...)  2
('J!=0','tiltingPlus',False,'outdeg',2,'longsq',False ...) 42
```

So 42 gate-admitted rejecting steps have out-degree 2 and no long square on `alg.rels`. This is a D-type (or G-type) reject reached by a guided walk. It contradicts "the walks reach no D or G type algebra" and "J != 0 <=> out-degree 1 and long square holds step by step", and it answers the author's own Next item 1. The claim is that the iff holds for n = 5..7 only, and it breaks at n = 8. Wording in the claim and the note for E-084/E-095/E-100 ("shape of the walk class") must not be extended past n = 7. I did not classify the 42 as D or as G (presentation artefact) or as something else. The dominance of out-degree 2 at n = 8 (42 against 2) makes it more likely a real D-type phenomenon than an artefact.

The hand-built analysis itself is sound: direction (1) follows from minimality; (2)(a) is correctly labelled a hypothesis. The n = 8 result supports the author's own caveat that (a) is not forced; it shows (a) is also not an invariant of reachable algebras. The text "an empirical regularity of the algebras reached" is false at n = 8. Caveat on (1): "out-degree 1 not needed" is asserted, not tested; case D shows J != 0 with all relations present, but no hand case isolates an out-degree 2 algebra with a minimal long relation into one arrow only (that is E, J = 0, consistent).

## New?

Grepped FINDINGS, HYPOTHESES, RETRACTIONS for "out-degree", "long-sided", "long square": nothing. In EXPERIMENTS: E-066, E-078, E-097, E-100 are as cited. Out-degree 1 versus 2 on the walk class is not recorded; the D/E/G split and the "short level is fine" one-liner are new. Walk counts duplicate E-100 (n = 5..7) and add nothing there.

## Evidenced?

Hand table: yes, specific and reproducible. Walk table: counts are cap-dependent lower bounds (steps, 300 s cap), as admitted, but the "0 in the J != 0 with out-degree 2" column is a negative over an uncapped-in-claim range (n = 5..7) that stops one size short of where it fails. The coverage range "n = 5..7, class 0 only" must appear in the claim line itself. n = 7 is a 300 s walk covering 35 094 algebras; no statement about completeness. Class 1 was not run.

## Required for acceptance

1. Add n = 8 (class 0, stated budget) to the walk table and retract or restrict "the walks reach no D or G type algebra" and the step-by-step iff to n <= 7.
2. Identify the 42 n = 8 out-degree 2 rejects: print one, say whether it is D-type (commutes into both arrows) or a G-type presentation artefact, and confirm it fails the Cartan test through the actual rewrite.
3. Rewrite the "Next" and chair wording: the iff does not hold on guided walks at n = 8; "long-sided square" is a shape at n <= 7 only.
4. State class (0 only) and caps in the Claim line, not only in the reproduction block.
