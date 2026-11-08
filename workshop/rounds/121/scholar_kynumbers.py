"""Round 121 (scholar), referee item 1: number theorem-like environments when they share the SUBSECTION counter
(\\newtheorem{lemma}[subsection]{Lemma}, as in Keller-Yang 0906.0761v3). The round-110 counter ignores \\subsection and
so prints 2.2 / 2.7 / 7.2 for T:DWZreduction / L:reduction / lemma:Mittag-Leffler; here \\subsection (unstarred) and
every environment sharing [subsection] step one counter, reset by \\section (unstarred). Comment text (%) ignored.
usage: .venv/bin/python -I workshop/rounds/121/scholar_kynumbers.py FILE.tex [LABEL_OR_NUMBER ...]
With extra args, prints only rows whose label or number equals an arg (a planted wrong number prints NOTHING)."""
import re, sys
src = open(sys.argv[1], encoding='utf-8', errors='replace').read()
src = re.sub(r'(?m)(?<!\\)%.*$', '', src)
shared = {m.group(1) for m in re.finditer(r'\\newtheorem\{(\w+)\}\[subsection\]', src)}
pat = re.compile(r'\\(section|subsection)(\*?)\{|\\begin\{(\w+)\}')
sec = 0; cnt = 0; rows = []
for m in pat.finditer(src):
    if m.group(1):
        if m.group(2): continue
        if m.group(1) == 'section': sec += 1; cnt = 0
        else: cnt += 1
        continue
    if m.group(3) in shared:
        cnt += 1
        lab = re.search(r'\\label\{([^}]*)\}', src[m.end():m.end() + 200])
        body = ' '.join(src[m.end():m.end() + 300].split())[:90]
        rows.append((f'{sec}.{cnt}', m.group(3), lab.group(1) if lab else '-', body))
want = sys.argv[2:]
for r in rows:
    if not want or r[0] in want or r[2] in want:
        print(*r, sep=' | ')
