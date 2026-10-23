"""Round 082 (idea): per-EIGENVALUE signature of G = C^-1 + C^-T under M = Phi^T (a finer invariant than the per-rational-factor signature
of E-209), then (stage 2) a full-lattice equivariant-isometry search via the M-invariant positive-definite majorant H = G (I - 2 pi_N).
usage (repo root): timeout 10m .venv/bin/python workshop/rounds/082/maverick_eigsig.py [--plan] [stage]   stage = sig (default) | iso
"""
import sys, time, collections
sys.path.insert(0, "workshop/rounds/074")
import numpy as np
from maverick_isometry import data, shortVectors
from quivermutation import coxeterTables as ct, freeMoves as fm

def eigsig(G, M, tol=1e-6):
    """dict: rounded eigenvalue angle/2pi -> sign of w* G w  (M diagonalizable, finite order => G-orthogonal eigenspaces for lam != conj-inverse)."""
    lam, W = np.linalg.eig(M.astype(float)); out = {}
    for k in range(len(lam)):
        w = W[:, k]; s = (w.conj() @ G.astype(float) @ w).real
        ang = round((np.angle(lam[k]) / (2 * np.pi)) % 1.0, 4)
        out[(ang, k)] = (np.sign(s) if abs(s) > tol else 0)
    return out
def sigProfile(G, M):
    """sorted tuple of (angle, #pos, #neg) over eigenvalue angles (multiplicity resolved)."""
    d = collections.defaultdict(lambda: [0, 0])
    for (ang, k), s in eigsig(G, M).items():
        d[ang][0 if s > 0 else 1] += 1
    return tuple(sorted((a, p, q) for a, (p, q) in d.items()))
def isoClasses15(N=11):
    lnas, orbits = fm.derivedOrbits(N, free=True); parent = {}
    for k, mem in orbits.items():
        for m in mem: parent[tuple(m)] = tuple(mem[0])
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for m in list(parent):
        a, b = find(m), find(fm.mirrorRow(N, m))
        if a != b: parent[a] = b
    keys = ct.lnaKeys(N); groups = collections.defaultdict(list)
    for rl in lnas:
        C = np.array(ct.lnaCartanMatrix(N, rl), float); Ci = np.linalg.inv(C)
        if np.linalg.eigvalsh(Ci + Ci.T).min() > -1e-9: continue
        Phi = np.rint(-Ci.T @ C).astype(np.int64); I = np.eye(N, dtype=np.int64); P = Phi.copy(); ok = False
        for m in range(1, 201):
            if (P == I).all(): ok = True; break
            P = P @ Phi
        if ok: groups[keys[rl]].append(tuple(rl))
    mem = next(v for v in groups.values() if len(v) == 15)
    reps = {}
    for m in mem: reps.setdefault(find(m), m)
    return mem, list(reps.values()), find
if __name__ == "__main__":
    plan = "--plan" in sys.argv; args = [a for a in sys.argv[1:] if not a.startswith("--")]
    t0 = time.time(); rd = lambda s: tuple(int(c) for c in s)
    if plan:
        print("plan: sig stage = load derivedOrbits n=11 (~50 s), 15 eigendecompositions 11x11 (ms); n=10 controls 90000000/00000030/D_10/50505000/34504030 (ms)"); sys.exit()
    for s in ["90000000", "00000030", "22222230", "00000000", "50505000", "34504030"]:
        G, M, _ = data(rd(s), 10); print("n=10", s, sigProfile(G, M))
    mem, reps, find = isoClasses15()
    prof = collections.Counter()
    for rl in mem:
        G, M, _ = data(rl, 11); p = sigProfile(G, M); prof[p] += 1
    print("15-group: distinct eigenvalue-signature profiles:", len(prof), f"{time.time()-t0:.0f}s")
    for p, c in prof.items(): print(c, p)
    print("9 reps:", ["".join(map(str, r)) for r in reps])
