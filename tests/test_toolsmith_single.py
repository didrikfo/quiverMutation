import subprocess, sys
def test_label_block_runs():
    r = subprocess.run([sys.executable, 'workshop/rounds/038/toolsmith_single.py'], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-500:]
    ls = [l for l in r.stdout.splitlines() if l.startswith('n=10 single 3')]
    assert len(ls) == 3
    v = [l.split(') (', 1)[1] for l in ls]
    assert v[0] == v[2] != v[1]   # I1 = {(h,K)=(4,2),(2,4)}, I2 = {(3,3)}
