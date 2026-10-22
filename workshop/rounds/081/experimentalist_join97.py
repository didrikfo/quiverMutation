"""Round 081 (experimentalist): budgeted driver for the n = 9 F-010 3060000 vs 3304000 J = 0 join at 7 + 7 (E-219 open item).
Runs rounds/079/skeptic_f010join.py `side` in resumable 540 s chunks (each < 10 min), A and B in parallel, until both balls reach DEPTH or --budget-hours is spent.
  experimentalist_join97.py plan                       print the sized plan (no work)
  experimentalist_join97.py run --budget-hours H [--depth 7] [--dir DIR]    exit 0 = both balls done and meet printed (result in DIR/meet.txt); exit 2 = budget spent, rerun to resume
Sizing (round 081, 4 cores shared, ~5-9 expansions/s per process): 3060000 ball needs 26 353-6 354 = ~20 000 further expansions to depth 7 (about 1-1.2 h), 3304000 ~11 400 (about 0.6 h); meet needs WL hashes of ~110 000 + 65 000 level-7 algebras (est. 0.5-1 h, unsized).
Run from the repository root."""
import sys, os, time, subprocess, argparse
ap = argparse.ArgumentParser(); ap.add_argument('mode'); ap.add_argument('--budget-hours', type=float, default=0.0); ap.add_argument('--depth', type=int, default=7); ap.add_argument('--dir', default='/tmp/join97')
a = ap.parse_args(); S = 'workshop/rounds/079/skeptic_f010join.py'
if a.mode == 'plan':
    print(__doc__); sys.exit(0)
os.makedirs(a.dir, exist_ok=True); T0 = time.time(); cmd = lambda root, tag: ['timeout', '10m', sys.executable, '-u', '-W', 'ignore', S, 'side', 'x', root, str(a.depth), '540', f'{a.dir}/ck_{tag}.pkl']
def done(tag, log):
    try: return ('completed depth %d ' % a.depth) in open(f'{a.dir}/{log}').read().splitlines()[-1]
    except Exception: return False
pending = {'A': '3060000', 'B': '3304000'}
while pending:
    if time.time() - T0 > a.budget_hours * 3600: print('budget spent, resume by rerunning'); sys.exit(2)
    ps = {t: subprocess.Popen(cmd(r, t), stdout=open(f'{a.dir}/log_{t}.txt', 'a'), stderr=subprocess.STDOUT) for t, r in pending.items()}
    for t, p in ps.items(): p.wait()
    for t in list(pending):
        if done(t, f'log_{t}.txt'): del pending[t]
out = subprocess.run(['timeout', '10m', sys.executable, '-u', '-W', 'ignore', S, 'meet', f'{a.dir}/ck_A.pkl', f'{a.dir}/ck_B.pkl'], capture_output=True, text=True)
open(f'{a.dir}/meet.txt', 'w').write(out.stdout + out.stderr); print(out.stdout); sys.exit(0 if out.returncode == 0 else 2)
