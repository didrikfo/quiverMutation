"""List the record entries (F-/H-/E-/R-) whose text mentions a primary source.
usage: .venv/bin/python -I workshop/rounds/137/scholar_citedeps.py [--plant]
Run from the repository root. For each source pattern, prints entry id, file and
the matching line (first 160 chars). An entry is the text from one '## X-NNN'
header to the next. --plant adds a fake source pattern that must match nothing
(control: the scan reports zero entries when the truth is zero)."""
import re, sys, glob

SOURCES = {
    'Nikulin 1979': r'Nikulin',
    'Happel-Seidel': r'Happel.{0,3}Seidel|Happel and Seidel|HS.symmetry',
    'BGS': r'\bBGS\b|Bobi.ski.{0,3}Geiss|Geiss.{0,3}Skowro',
    'Anick 1982': r'Anick',
}
if '--plant' in sys.argv:
    SOURCES['planted (must be empty)'] = r'Zzyzx-Qwertyuiop'

files = ['research/FINDINGS.md', 'research/HYPOTHESES.md', 'research/EXPERIMENTS.md',
         'research/RETRACTIONS.md']
hdr = re.compile(r'^##\s+([FHER]-\d+)')
for name, pat in SOURCES.items():
    rx = re.compile(pat, re.I)
    hits = {}
    for f in files:
        cur = None
        for ln, line in enumerate(open(f, encoding='utf-8'), 1):
            m = hdr.match(line)
            if m:
                cur = m.group(1)
            if rx.search(line):
                hits.setdefault((cur or 'preamble', f), []).append((ln, line.strip()[:160]))
    print('####', name, '--', len(hits), 'entries')
    for (eid, f), ls in sorted(hits.items(), key=lambda kv: kv[0]):
        for ln, s in ls[:2]:
            print(f'  {eid:8s} {f.split("/")[-1]}:{ln}  {s}')
