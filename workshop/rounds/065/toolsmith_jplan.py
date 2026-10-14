"""Round 065 (toolsmith), T10 item 1: size a ball search for a J != 0 path at n = 8, 9.
Wraps rounds/047/experimentalist_deepreplay.py (key-guarded BFS, per-edge tiltingPlus tally) and prints a sizing row.
  toolsmith_jplan.py run N CLS BUDGET_SEC [DEPTH]   run one slice (checkpoint in /tmp/tsm065), parse, print a row
  toolsmith_jplan.py parse FILE...                  parse saved deepreplay outputs
Row: levels (cumulative expansions per depth, seconds), key-kept edges, key-kept tiltingPlus-failing edges (J != 0), s per expansion,
and, when the walk reached a full level, the projected cost of the next level (x2.4 expansions per depth, as measured)."""
import sys, re, subprocess, os
def parse(txt, label):
    lv = [(int(m[0]), int(m[1]), int(m[2]), int(m[3])) for m in re.findall(r'^depth (\d+) expanded (\d+) next (\d+) .*? (\d+)s$', txt, re.M)]
    g = re.search(r"\('guard', 'tilt'\) (\d+)", txt); gb = re.search(r"\('guard', 'NOTtilt'\) (\d+)", txt)
    ed = int(g.group(1)) if g else 0; bad = int(gb.group(1)) if gb else 0
    nk = re.search(r'^n (\d+) class (\d+) .* starts (\d+) distinct algebras (\d+)', txt, re.M)
    print(label, 'n,cls,starts,distinct =', nk.groups() if nk else None)
    print('  levels (depth, cum. expansions, next frontier, s):', lv)
    print('  key-kept edges %d, of which tiltingPlus-false (J != 0) %d' % (ed + bad, bad))
    if lv:
        d, e, nx, s = lv[-1]; print('  s/expansion %.3f; next level %d expansions ~ %.0f s (frontier as printed)' % (s / max(e, 1), nx, nx * s / max(e, 1)))
    return ed + bad, bad
if __name__ == '__main__':
    if sys.argv[1] == 'run':
        n, c, b = sys.argv[2:5]; dep = sys.argv[5] if len(sys.argv) > 5 else '9'
        os.makedirs('/tmp/tsm065', exist_ok=True)
        p = subprocess.run(['.venv/bin/python', '-u', 'workshop/rounds/047/experimentalist_deepreplay.py', n, '--class', c, '--depth', dep, '--budget-sec', b, '--ckpt', '/tmp/tsm065/n%sc%s.ckpt' % (n, c)], capture_output=True, text=True)
        parse(p.stdout, 'n%s class %s' % (n, c))
    else:
        tot = [0, 0]
        for f in sys.argv[2:]:
            a, b = parse(open(f).read(), f); tot[0] += a; tot[1] += b
        print('TOTAL key-kept edges', tot[0], 'J != 0', tot[1])
