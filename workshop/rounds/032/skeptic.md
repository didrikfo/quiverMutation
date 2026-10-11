# Skeptic round 032

## Most promising question

**Can we hand-build a both-die mutation sequence at n=6 under a capped walk, and what does the outcome tell us?** If it exists, E-122's claim fails and the "no circuit" argument loses its wall. If it doesn't, why not: walk-reachability, algebra, or definition?

## Weakest claim

E-122 states capped walks at n=6,7 contain no both-die row (in contrast to loose walks); this is presented as an observed fact without mechanistic explanation. Loose rows have it (both-die at n=6,7; half-W elsewhere); capped rows don't. But capped and loose are shape definitions only—no algebraic property is given for why the shape forbids the both-die pattern. Risk: the distinction is vacuous or the query missed cases.

## What I need

- **experimentalist:** Build an n=6 both-die square (two mutations killing both kernel terms, verified by code). Does the walk reach it when paths go under a capped core (height >= 2 wall)? If yes, E-122 is false. If no: why—Cartan fails, the path gets stuck, or algebra forbids it?

- **theorist:** Why do capped cores structurally prevent both-die rows? If it follows from the kernel structure (two-term alone cannot forbid it per E-115), then it's walk geometry or the algebra's height. State it precisely.
