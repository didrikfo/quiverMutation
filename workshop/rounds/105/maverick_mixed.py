"""Round 105 maverick: H-024 control.  E-269 exhaustive Phi-commutant Z-congruence test (rounds/099/maverick_iso140.py test(), threshold l.16 fixed to > 1+1e-6,
closed-form bound n-2+lam+1/lam+slack) on the 6 E-275 mixed-digit derived-equivalent pairs (YES expected), then the NO control (E-269 cross pairs, key group of 140, n=12).
Digit strings are relLengths (n-2 entries, 'a' not needed: all digits < 10 except none).  usage (repo root): timeout 10m .venv/bin/python -W ignore workshop/rounds/105/maverick_mixed.py [mixed|ctl]"""
import sys; sys.path.insert(0, "workshop/rounds/099"); sys.path.insert(0, "workshop/rounds/097")
import maverick_iso140 as M
PAIRS = [(9, "3030400", "4003030"), (10, "30305000", "50003030"), (11, "303030400", "400303030"), (11, "304005000", "003040030"),
         (12, "3030305000", "5000303030"), (12, "3040050000", "0040040030")]
dg = lambda s: tuple(int(c) for c in s)
def run(n, a, b, slack):
    M.N = n
    try:
        r = M.test(a, b, slack); return {k: v for k, v in r.items() if k != "P"}
    except AssertionError as e: return "NOT APPLICABLE (assert: signature/det/bound) " + repr(e)[:120]
mode = sys.argv[1] if len(sys.argv) > 1 else "mixed"
if mode == "mixed":
    for n, a, b in PAIRS:
        for slack in (0.5, 30):
            print(n, a, b, "slack", slack, run(n, dg(a), dg(b), slack), flush=True)
else:
    import itertools, collections
    from quivermutation import coxeterTables as ct
    from quivermutation.quipuForms import quipuForAlmostSeparateLNA
    M.N = 12; keys = ct.lnaKeys(12); groups = collections.defaultdict(list)
    for rl, kk in keys.items(): groups[kk].append(rl)
    a0 = (2, 2, 2, 3, 0, 2, 3, 0, 3, 0); grp = [v for v in groups.values() if a0 in map(tuple, v)][0]
    q = collections.defaultdict(list)
    for rl in grp:
        r = quipuForAlmostSeparateLNA(12, list(rl))
        if r is not None: q[r].append(tuple(rl))
    reps = [(r, rl) for r in q for rl in q[r][:2]]
    for (r1, a), (r2, b) in itertools.combinations(reps, 2):
        print("within" if r1 == r2 else "ACROSS", a, b, run(12, a, b, 0.5), flush=True)
