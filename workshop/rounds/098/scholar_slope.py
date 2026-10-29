"""Tabulate (a, b, b/a) per key group from the recorded E-253 shift logs.
Reads workshop/rounds/091/maverick_shift*.out only; no computation of algebras.
Run from repo root: .venv/bin/python workshop/rounds/098/scholar_slope.py"""
import re, glob
from fractions import Fraction
for f in sorted(glob.glob("workshop/rounds/091/maverick_shift*.out")):
    mem, summ = [], None
    for line in open(f):
        m = re.match(r"member (\d+) first k=(\d+) .*shifts=\[(-?\d+)\]", line)
        if m: mem.append((int(m[2]), int(m[3])))
        if line.startswith("SUMMARY"): summ = line.strip()
    pairs = set(mem)
    print(f.split("/")[-1], "members", len(mem), "distinct (a,b):", sorted(pairs),
          "slopes:", sorted({Fraction(-b, a) for a, b in pairs}))
    print("  ", summ)
    sm = re.search(r"order/sign \((\d+), (-?\d+)\)", summ)
    order, sign = int(sm[1]), int(sm[2])
    for a, b in pairs:
        ok = (-1) ** (a + b) == sign and a == order
        print("   a+b parity:", (a + b) % 2, "(-1)^(a+b) =", (-1) ** (a + b),
              "logged sign", sign, "order", order, "PASS" if ok else "FAIL")
