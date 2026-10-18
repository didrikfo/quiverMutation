"""Overnight item: n = 9 F-010 NO control (3060000 vs 3304000) at depth 5+5 and beyond, iso meets.
usage: .venv/bin/python -W ignore workshop/rounds/073/toolsmith_f010_overnight.py PKL --budget-hours H [--depth D]
Runs toolsmith_certificate.py `f010 D 2` (needs c2.pkl from 049/toolsmith_collect.py 7 2 20000 PKL 100). Exit 0 finished, 2 budget spent, 1 error."""
import argparse, subprocess, sys, os, time
ap = argparse.ArgumentParser(); ap.add_argument('pkl'); ap.add_argument('--budget-hours', type=float, required=True); ap.add_argument('--depth', type=int, default=5)
a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__)); t = time.time()
try:
    r = subprocess.run([sys.executable, '-W', 'ignore', os.path.join(here, 'toolsmith_certificate.py'), a.pkl, '2', '4', 'f010', str(a.depth), '2'], timeout=a.budget_hours * 3600)
except subprocess.TimeoutExpired:
    print('BUDGET SPENT after %.2f h at depth %d+%d: no verdict (a bound, not a proof)' % ((time.time() - t) / 3600, a.depth, a.depth)); sys.exit(2)
sys.exit(0 if r.returncode == 0 else 1)
