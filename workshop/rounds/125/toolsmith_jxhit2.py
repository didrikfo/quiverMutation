"""Round 125 (toolsmith): round-122 `toolsmith_jxhit.py` (brute-force Cartan, det C, Coxeter polynomial, SNF(S) by sympy and by determinantal
divisors, J = 0 balls of the children, the explicit walk) run unchanged from another root of n = 11 group 5 with the same four steps
(root F3 R8 F3 then F3 / R5; dual F8 R3 R3 then F5 / R3).  usage (repo root):
  ROOT=203033030 R=2 WALK=1 timeout 10m .venv/bin/python -W ignore workshop/rounds/125/toolsmith_jxhit2.py"""
import os
src = open('workshop/rounds/122/toolsmith_jxhit.py').read()
src = src.replace("ROOT = '003033030'", "ROOT = os.environ.get('ROOT', '203033030')", 1)
assert "os.environ.get('ROOT'" in src
exec(compile(src, 'jxhit', 'exec'))
