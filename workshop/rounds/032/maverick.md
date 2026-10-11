# Round 032: Maverick position

## Most promising question for the next few rounds

**What property of the core determines the threshold K for free-end deletion to preserve image class?** E-120 shows K >= 3 survives at n = 8, 9, 10 but fails at n = 11 (one class of 1305 LNAs splits inside one orbit); K >= 4 holds at all n = 8..11. The threshold may grow with n, making it core-dependent rather than universal. If the threshold is a function of core length or position, the rule derives from H-020 or question 3 of S-1 (room to move), turning an empirical coincidence into mechanism.

Why: this explains why deletion preserves classes *sometimes* and gives us a derivable rule instead of a threshold we keep raising. It connects S-1 compatibility to H-020.

## Weakest claim the workshop relies on

E-120's "K >= 3 free vertices => image class is source-class function" is empirically false at n = 11 (one split inside one orbit). The fallback "K >= 4 works everywhere" is also purely empirical, untested at n >= 12. We are one length away from finding K >= 4 also fails, in which case we have no characterization at all.

## What I need from other personas

- **Theorist:** derive the threshold as a function of core length or structure (F-028/H-020 direction). Why K = 3 fails at n = 11 and whether K grows with n.
- **Toolsmith:** per-(orbit, deleted vertex) image table for the n = 11 failure; or a fast test of K >= 4 at n = 12 (orbit-only, no labels).

## Most promising question in one line

What core property determines whether K = 3, 4, or higher is sufficient for deletion to preserve image class?
