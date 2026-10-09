"""Per-n table for the 7 survivors: how many pair (fit d<=6) and hold strict mirror.
Usage: .venv/bin/python workshop/rounds/003/experimentalist_table.py workshop/rounds/003/experimentalist_census_7cores_n8_18.jsonl
Reuses the fit rule of workshop/rounds/002/experimentalist_fit.py (run it per n for the full lines)."""
import sys, json, subprocess, collections, tempfile, re
rows = [json.loads(l) for l in open(sys.argv[1])]
by = collections.defaultdict(list)
for r in rows: by[r['n']].append(r)
for n in sorted(by):
    with tempfile.NamedTemporaryFile('w', suffix='.jsonl', delete=False) as f:
        f.write(''.join(json.dumps(r) + '\n' for r in by[n]))
    out = subprocess.run([sys.executable, 'workshop/rounds/002/experimentalist_fit.py', f.name],
                         capture_output=True, text=True).stdout
    lines = [l for l in out.splitlines() if re.match(r'^\d+ ', l)]
    fit = [l.split()[0] for l in lines if 'd=None' not in l]
    strict = [l.split()[0] for l in lines if 'strict=1' in l]
    print(f"n={n:2d} cores={len(lines)} pair={len(fit)} strict={len(strict)} nofit={[l.split()[0] for l in lines if 'd=None' in l]}")
