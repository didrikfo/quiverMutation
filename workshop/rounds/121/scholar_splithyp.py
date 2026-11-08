"""Round 121 (scholar), H-031 / H-023 literature check: test, on the mutated dg quivers kM of the loop models
(workshop/rounds/106/theorist_loopfilt.py; Def. 3.2 of arXiv:1504.02617v1), the hypotheses of the published
trivial-part / splitting statements:
  DWZ 0704.0649v4 Def. 4.1 (no loops, eq. no-loops), Thm 4.6 (Splitting Theorem; complete path algebra, potential);
  Keller-Yang 0906.0761v3 Thm 2.3 ("all quivers are quivers without loops"), Lemma 2.10 (complete Ginzburg dg algebra of a
  DIRECT SUM (Q',W') + trivial (Q'',W''): projection Gamma -> Gamma' is a q.i.).
Mechanical checks per model (H = homological degree as in Oppermann):
  (a) loops of kM by degree; degree-0 loops; (b) Ginzburg shape: a degree-2 loop at every vertex and #deg1 = #deg0 letters;
  (c) linear part of d: which letters have a linear term (KY: 'trivial' iff that part pairs letters bijectively);
  (d) KY Lemma 2.10 direct-sum hypothesis on the RAW generators: C = letters with a linear term + their linear targets;
      d(C) in T(C) and no d(rest) uses a C letter?  (e) the 'reduced' remainder: does the rest still contain Theta letters?
Control (power): the same test (d) on the theorist's substituted M5 generators (workshop/rounds/118/theorist_m45.py, E-309)
must say YES (E = trivial part is a direct summand) while raw kM says NO.
Usage (repo root): .venv/bin/python -W ignore workshop/rounds/121/scholar_splithyp.py"""
import importlib.util, io, contextlib, sys
sp = importlib.util.spec_from_file_location('lf', 'workshop/rounds/106/theorist_loopfilt.py')
lf = importlib.util.module_from_spec(sp); sp.loader.exec_module(lf)

def name(L):
    if L[0] == 's': return L[1] + '*'
    _, a, p, b = L
    return (a or '') + p + ('/' + b if b else '')

def split_test(gens, dval, theta):
    """gens: letters; dval: letter -> {word: coef}; theta: set of Theta letters (to be cancelled).
    Returns (C, mixing list, trivial_is_linear, rest_thetas)."""
    lin = {x: {w[0]: k for w, k in dval[x].items() if len(w) == 1 and k} for x in gens}
    C = set()
    for x in gens:
        if lin[x]: C.add(x); C.update(lin[x])
    mix = []
    for x in gens:
        for w in dval[x]:
            uses = set(w) & C
            if x in C and set(w) - C: mix.append((name(x) if not isinstance(x, str) else x, 'd uses non-C', w))
            if x not in C and uses: mix.append((name(x) if not isinstance(x, str) else x, 'd uses C', w))
    # trivial part linear and bijective: every C letter with a linear term has exactly one, and targets distinct
    tgts = [t for x in C for t in lin[x]]
    bij = all(len(lin[x]) <= 1 for x in C) and len(tgts) == len(set(tgts)) and len(tgts) * 2 == len(C)
    rest_th = sorted((name(x) if not isinstance(x, str) else x) for x in theta if x not in C)
    return C, mix, bij, rest_th

def ideal_test(gens, dval, C):
    """Referee item 4 (KY's PROOF, not its statement): the projection killing C is a dg map iff the two-sided ideal (C)
    is d-stable, i.e. every word of d(c), c in C, contains a C letter. Returns the failing (c, word) pairs."""
    return [(c, w) for c in gens if c in C for w, k in dval[c].items() if k and not set(w) & C]

def run(model):
    q = lf.MODELS[model](); m, letters, wt = lf.setup(q)
    dval = {L: {w: k for w, k in m.partial({(L,): 1}).items() if k} for L in letters}
    deg = {L: m.deg(L) for L in letters}
    verts = sorted({m.src(L) for L in letters} | {m.tgt(L) for L in letters})
    loops = sorted((deg[L], name(L)) for L in letters if m.src(L) == m.tgt(L))
    n0 = sum(1 for L in letters if deg[L] == 0); n1 = sum(1 for L in letters if deg[L] == 1)
    ginz = all(any(m.src(L) == v == m.tgt(L) and deg[L] == 2 for L in letters) for v in verts) and n1 == n0
    theta = {L for L in letters if L[0] == 'm' and L[3] is not None}
    C, mix, bij, rest = split_test(letters, dval, theta)
    print(f'== {model}: {len(letters)} letters, degrees {sorted(set(deg.values()))}, vertices {verts}')
    print(f'  (a) loops (deg, letter): {loops}; degree-0 loops: {[l for l in loops if l[0] == 0]}')
    print(f'  (b) Ginzburg shape (deg-2 loop at every vertex and #deg1 = #deg0 = {n0}/{n1}): {ginz}')
    print(f'  (c) letters with a linear term in d: {[(name(x), [name(w[0]) for w in dval[x] if len(w) == 1]) for x in letters if any(len(w) == 1 for w in dval[x])]}')
    print(f'  (d) KY L.2.10 direct-sum on raw generators: C = {sorted(name(x) for x in C)}; trivial part linear+bijective: {bij}; mixing terms: {len(mix)}  -> {"YES" if C and not mix and bij else "NO"}')
    for t in mix[:4]: print('      e.g.', t[0], t[1], [name(y) for y in t[2]])
    bad_id = ideal_test(letters, dval, C)
    print(f'  (d\') dg-ideal version: words of d(C) outside the ideal (C): {len(bad_id)}  -> {"YES" if C and not bad_id else "NO"}', [(name(c), [name(y) for y in w]) for c, w in bad_id])
    print(f'  (e) Theta letters outside the trivial part (reduced part, quadratic d, Lemma A territory): {rest}')
    return m, letters, dval

for mod in ['M1_loop1_dphi=ba', 'M2_loop1_dphi=0', 'M4_two_alpha', 'M5_loop_two_paths', 'M6_noloop_control']:
    run(mod)

# control: theorist's substituted M5 generators (E-309): E = trivial part, B = rest
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    sys.path.insert(0, 'workshop/rounds/118'); import theorist_m45 as T
gens = list(T.claim5); theta = {'ta', 'saa', 'sma', 'tmp', 'smmp', 'samp'}
C, mix, bij, rest = split_test(gens, T.claim5, theta)
print(f'== CONTROL M5 after the E-309 substitution: C = {sorted(C)}; linear+bijective {bij}; mixing {len(mix)} -> {"YES" if C and not mix and bij else "NO"}; Theta left in reduced part: {rest}; dg-ideal failures {len(ideal_test(gens, T.claim5, C))}')
# planted error: re-introduce one mixing term (d ta = pa + w fb) -> must say NO
bad = dict(T.claim5); bad['ta'] = {('pa',): 1, ('w', 'fb'): 1}
C2, mix2, bij2, _ = split_test(gens, bad, theta)
print(f'== PLANTED (d ta = pa + w fb): mixing {len(mix2)} -> {"YES" if C2 and not mix2 and bij2 else "NO"}')
# planted for (d'): a word of d(ta) avoiding C (w fb is in (C) only if fb is; use two non-C letters)
nonC = sorted(x for x in gens if x not in C)[:2]
bad2 = dict(T.claim5); bad2['ta'] = {('pa',): 1, tuple(nonC): 1}
print(f'== PLANTED for (d\') (d ta = pa + {" ".join(nonC)}): dg-ideal failures {len(ideal_test(gens, bad2, C))}')
