"""Round 085 (scholar): what F-047 names 2400230 / 3304000 / 3060000 by.
Prints the almost-separate flag, the quipu (thm:QuipuToAn of arXiv:2305.06642 inverted
in quipuForms) and whether classPreservingOrbit (corollary EquivNakayamaAlgebras) joins them."""
from quivermutation import nakayama as nk
names = ["2400230", "3304000", "3060000"]
algs = {s: nk.LinearNakayamaAlgebra.fromClassName(s) for s in names}
for s, a in algs.items():
    print(s, "n =", a.length if hasattr(a, "length") else "?", "relLengths", a.relationString(),
          "almost separate:", a.hasAlmostSeparateRelations(), "quipu:", a.quipuName())
    orb = {x.relationString() for x in a.classPreservingOrbit()}
    print("   orbit size", len(orb), "contains:", [t for t in names if t != s and t in orb])

# the project's wider move orbit (lnaMoves.orbitOf): which members are almost separate, and what they are named
from quivermutation import lnaMoves, lines
for s in names:
    a = algs[s]
    orb = lnaMoves.orbitOf(9, a.relLengths)
    out = {}
    for r in orb:
        m = nk.LinearNakayamaAlgebra.fromClassName(r)
        out[m.relationString()] = (m.hasAlmostSeparateRelations(), m.quipuName() if m.hasAlmostSeparateRelations() else "")
    print(s, "lnaMoves orbit size", len(orb))
    named = {}
    for k, (ok, q) in out.items():
        if ok: named.setdefault(q, []).append(k)
    print("   almost-separate members by quipu:", {q: len(v) for q, v in named.items()}, "| 2400230 in orbit:", any(k == algs["2400230"].relationString() for k in out), "| 3304000 in orbit:", any(k == algs["3304000"].relationString() for k in out))

# are 3004000 / 3024000 (almost separate, same quipu as 2400230) joined to 2400230 by Fosse's class-preserving ops?
print("--- Fosse class-preserving closure (cor:EquivNakayamaAlgebras ops) ---")
for s in ["2400230", "3004000", "3024000"]:
    a = nk.LinearNakayamaAlgebra.fromClassName(s)
    seen = {a.relationString()}; fr = [a]
    while fr:
        x = fr.pop()
        for y in x.classPreservingOrbit():
            k = y.relationString()
            if k not in seen: seen.add(k); fr.append(y)
    print(s, a.quipuName() if a.hasAlmostSeparateRelations() else "not almost separate", "closure size", len(seen),
          "contains 2400230:", algs["2400230"].relationString() in seen)
