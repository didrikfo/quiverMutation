# Theorist round 024

## Most promising question

**When does the peeling formula fail, and what does the mirror-chain rule fix?**
Why: E-103's peeling (depth 1 + min(a,b)) works for single-big-relation LNAs but fails when the blocker is itself a big relation (1/6/24 mismatches at n = 8/9/10). E-106 fits a mirror-chain rule with 0 mismatches there and 13/13 out-of-sample predictions n = 10..12. The blocking rule (some big relation not blocked) holds everywhere. The open question is whether mirror-chain *replaces* peeling or *extends* it: is there a unified criterion that covers both, or are there three cases (none/one/two+ big relations) each with its own depth formula? This unifies T6.

## Weakest claim

The mirror-chain rule (E-106, fitted). It works on 24 LNAs at n = 8/9/10 and predicts 13 unseen cases at n = 10..12, but (1) no proof, (2) unclear when it applies (does it require the blocker to be big?), (3) unknown at n >= 13. If n = 13 breaks it, the whole cord-depth story remains a description, not a mechanism.

## Need

From maverick: which of the 24 failing LNAs have their blocker as a big relation? Can you sort them by blocker shape (how many arrows, what relations touch it)?
From experimentalist: n = 11/12 cord data if available, tested against mirror-chain.
