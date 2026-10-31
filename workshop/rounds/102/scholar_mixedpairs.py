"""Round 102 scholar: mixed-digit derived-equivalent LNA pairs at n = 9..12 from Fosse arXiv:2305.06642 (Cor 4.3 = corollary:2Rels, Cor 4.6 = cor:EquivNakayamaAlgebras).
The 8 algebras of Cor 4.6 are re-implemented here DIRECTLY FROM THE PRINTED FORMULAS (1-based start n_i, arrow count l_i), not from LinearNakayamaAlgebra.classPreservingOrbit.
Checks per source algebra: (a) printed members admissible?, (b) same Coxeter key (invariants.coxeterKey), (c) agreement with the repo orbit,
(d) Cartan-level tests of E-274 (Phi charpoly, q_S, Z[Phi]-module search) on one pair per source algebra. MIXED = >= 2 distinct nonzero digits.
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/102/scholar_mixedpairs.py [nmin nmax]"""
import sys, itertools
sys.path.insert(0, "workshop/rounds/097")
from quivermutation.nakayama import LinearNakayamaAlgebra as L
from quivermutation import invariants
import sympy as sp


def rels(a): return [(s, l) for s, l in a.relations()]
def build(n, R):
    """R = list of (start, arrows); None unless a valid LNA with the paper's standing assumptions."""
    R = sorted(R)
    if any(s < 1 or l < 2 or s + l > n for s, l in R): return None
    if len({s for s, _ in R}) != len(R): return None
    if any(s + l >= s2 + l2 for (s, l), (s2, l2) in zip(R, R[1:])): return None
    d = [0] * (n - 2)
    for s, l in R: d[s - 1] = l
    return L(n, d)

def cor46(n, R, typo=False):
    """The 8 algebras of Cor 4.6 as printed. R = [(n_i, l_i)]. typo=True uses the printed 'n_0+l_0' in the last length of the third/fourth right-hand entries."""
    ns = [s for s, _ in R]; ls = [l for _, l in R]; r = len(R) - 1
    e = lambda i: n - ns[i] - ls[i] + 1            # n - n_i - l_i + 1
    rev_starts = [e(i) for i in range(r, -1, -1)]  # (n-n_r-l_r+1, ..., n-n_0-l_0+1)
    rev_len = ls[::-1]
    lastL = n - ns[r] - ls[r] + 3
    out = {}
    out["1L"] = [(ns[i], ls[i]) for i in range(r + 1)]
    out["1R"] = list(zip(rev_starts, rev_len))
    out["2L"] = [(ns[i], ls[i]) for i in range(r)] + [(ns[r], lastL)]
    out["2R"] = list(zip([ls[r] - 2] + rev_starts[1:], [lastL] + rev_len[1:]))
    out["3L"] = [(ls[0] - 2, ns[0] + 2)] + [(ns[i], ls[i]) for i in range(1, r + 1)]
    lastR3 = ns[0] + ls[0] if typo else ns[0] + 2
    out["3R"] = list(zip(rev_starts, rev_len[:-1] + [lastR3]))
    out["4L"] = [(ls[0] - 2, ns[0] + 2)] + [(ns[i], ls[i]) for i in range(1, r)] + [(ns[r], lastL)]
    out["4R"] = list(zip([ls[r] - 2] + rev_starts[1:], [lastL] + rev_len[1:-1] + [lastR3]))
    return out

def key(a): return tuple(invariants.coxeterKey(a))
nmin, nmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (9, 12)
tot = bad_printed = bad_fixed = repo_dis = 0; examples = []
for n in range(nmin, nmax + 1):
    for a in L.allOfLength(n):
        R = rels(a)
        if len(R) < 2 or any(l < 3 for _, l in R) or not a.hasAlmostSeparateRelations(): continue
        if len({l for _, l in R}) < 2: continue          # mixed digits
        if a != min(a.classPreservingOrbit(), key=lambda x: x.relLengths) and False: pass
        k0 = key(a); orb = set()
        for typo, bucket in ((True, "printed"), (False, "fixed")):
            for name, RR in cor46(n, R, typo).items():
                b = build(n, RR)
                if b is None or key(b) != k0:
                    if typo: bad_printed += 1
                    else: bad_fixed += 1
                    if typo and len(examples) < 6: examples.append((n, a.className(), name, RR, None if b is None else b.className(), b is not None))
                elif not typo: orb.add(b)
        tot += 1
        if orb != set(a.classPreservingOrbit()): repo_dis += 1
print("n=%d..%d mixed sources %d; members failing as printed %d; failing with n_0+2 %d; sources where the 8 fixed members != repo orbit %d" % (nmin, nmax, tot, bad_printed, bad_fixed, repo_dis))
for ex in examples: print("printed-failure", ex)

# --- breakdown of the printed failures, then Cartan-level controls on a sample ---
import types, random, collections
nn = mm = 0
for n in range(nmin, nmax + 1):
    for a in L.allOfLength(n):
        R = rels(a)
        if len(R) < 2 or any(l < 3 for _, l in R) or not a.hasAlmostSeparateRelations() or len({l for _, l in R}) < 2: continue
        for name, RR in cor46(n, R, True).items():
            b = build(n, RR)
            if b is None: nn += 1
            elif key(b) != key(a): mm += 1
print("printed members that are not admissible LNAs: %d; admissible but different Coxeter key: %d" % (nn, mm))
sys.modules["cypari"] = types.SimpleNamespace(pari=None)       # cypari absent in .venv; module_iso_search unusable, the other tests need no pari
sys.path.insert(0, "workshop/rounds/097")
import maverick_cartanZ as mc
random.seed(1); pool = []
for n in range(9, 13):
    cands = [a for a in L.allOfLength(n) if len(rels(a)) >= 2 and all(l >= 3 for _, l in rels(a)) and a.hasAlmostSeparateRelations() and len({l for _, l in rels(a)}) >= 2]
    random.shuffle(cands); pool += [(n, a) for a in cands[:6]]
x = sp.symbols("x"); stat = collections.Counter()
for n, a in pool:
    orb = sorted(a.classPreservingOrbit(), key=lambda t: t.relLengths)
    # partner: another member of the orbit, or a itself plus a length-2 relation (Cor 4.3) if a free start exists
    partners = [b for b in orb if b != a]
    for s in range(1, n - 1):
        if a.relLengths[s - 1] == 0:
            d = list(a.relLengths); d[s - 1] = 2
            b = build(n, [(i + 1, v) for i, v in enumerate(d) if v])
            if b is not None and b.hasAlmostSeparateRelations(): partners.append(b); break
    for b in partners[:2]:
        Ca, Cb = mc.cartan(n, a.relLengths), mc.cartan(n, b.relLengths); Sa, Sb = Ca + Ca.T, Cb + Cb.T
        cp = sp.expand(mc.phi(Ca).charpoly(x).as_expr() - mc.phi(Cb).charpoly(x).as_expr()) == 0
        da, db = int(Sa.det()), int(Sb.det())
        try: q = mc.disc_iso(Sa, Sb)[0] if abs(da) > 1 else "vacuous"
        except Exception as e: q = "err"
        stat[(cp, mc.sig(Sa) == mc.sig(Sb), da == db, q)] += 1
        print(n, a.className(), b.className(), "charpoly_eq", cp, "sig_eq", mc.sig(Sa) == mc.sig(Sb), "det", da, db, "qS", q, flush=True)
print(stat)
