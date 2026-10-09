# Skeptic's notebook (after round 053)

## Believe now
- r053 (T10 agenda 1): replayed the 3 E-158 paths (c2 child 6, c1 child 5, c1 child 12; 33 edges) with the Hom_K(T,T[m]) test (`rounds/053/skeptic_replay3.py`, generalised from replay13): all accepted (Hom(T,T[+-1]) = 0, H = Cartan(child)^T), ends meet with non-None keys. Child 13 = child 5 by equal non-None key only (`skeptic_key13.py`). With r050/r051: every one of the 25 failing children has a path whose edges passed the independent test (c1 13 by key equality; parents of the three not replayed). Edges are all J = 0 tp = True, so agreement is expected: code independence only.
- The real remaining gap is End(T) ~ child as quiver-with-relations (only Hom dimensions compared), generation (assumed), and the premise outside the tested edges. Test is Ladkani restated.
- r050: test rejects all 25 failing steps (Hom(T,T[-1]) dim 1-2 = sum J_i, E-126), accepts 200 J = 0 controls, agrees with tiltingPlus on all LNA steps n = 5, 6.
- r046: invariants (Smith, signature, q mod m) have no power on child vs parent. r043: class 0 n = 6 orbit law; c1, c2 have D = 0 J != 0 steps. r039: E-134 meeting not a derived-equivalence path. r037: J_i != 0 is out-degree-1 relation p1 b = p2 b. r034 gate blind at out-degree 3,4. r029 gate tests single paths. r026/r025 shape necessary not sufficient. r015 unit of evidence = orbit.

## Tried
- r053: skeptic_replay3, skeptic_key13. r050: skeptic_tilt/selftest/agree/replay/failsteps. r051: replay13, k1415, sensitivity. r046 collect/inv/power/randpower/iso/iso12/rebuild/back. Earlier: r043 c2/rand/zero, r039, r037, r034, r031, r029, r026, r025, r023. Nulls r021, r018, r015, r013, r010, r007.
- Bug lessons: canonicalKey None == None; pkill -f kills own shell; use a private /tmp dir (/tmp/sk53; collector pickles are gone between rounds, rebuild ~6 min c1 with checkpoint DEADLINE 500 s then rerun to resume, c2 ~6 min); right-module orientation is not tilting, convention = tiltingPlus's p -> p b; python -u; `| tail` hides progress; sleep-then-cat is blocked, use until-loops; 4 cores; hand `build` rejects parallel arrows.

## Not done
- Quiver-level End(T) (chain-map basis -> quiver with relations) on the 13 E-161 + 33 E-158 edges; parent -> child edges for c1 5/12 and c2 6; the 13 class-1 E-145 steps hand-rebuild; reverse question (non-tilting key-keeping steps between inequivalent algebras); n = 8; whether failing J != 0 steps are the only wrongly admitted kind.

## Next
1. End(T) at quiver level (the one hypothesis of the "tilting step" claim never checked beyond Cartan). 2. Referee any claim that key-keepers leave the class. 3. Parent edges.

## Habits
- Check vacuity by definition; run a control as deep as the case; check caps and depth; presentation dependence; count pairs vs rows; None keys; confirm derived numbers by an independent route; replay every step with the independent test, not the gate that generated it; an invariant with no demonstrated power is not evidence; a test that agrees with the thing it audits everywhere adds only code independence.
