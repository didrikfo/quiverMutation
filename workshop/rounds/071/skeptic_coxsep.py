"""Separation power of the Coxeter polynomial at length 10: fraction of random distinct LNA pairs with different polynomials.
usage: skeptic_coxsep.py [seed]   (300 pairs, seed default 0)"""
import sys, random, importlib.util, os
sys.argv = [sys.argv[0], '10', '7'] + sys.argv[1:2]
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'skeptic_verify.py')).read().split("rules = [x for x")[0]
exec(src)
rows = enum(10); random.seed(seed); n = d = 0
while n < 300:
    a, b = random.sample(rows, 2); n += 1; d += cox(10, a) != cox(10, b)
print("seed", seed, "pairs", n, "different Coxeter polynomial", d)
