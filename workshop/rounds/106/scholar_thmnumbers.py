"""Number theorem-like environments of a LaTeX file by counter order.
usage: python -I scholar_thmnumbers.py FILE.tex [LABEL_OR_NUMBER ...]
Handles \\newtheorem{env}{Name}[section] (own base), \\newtheorem{env}[base]{Name}
(shared with base) and starred (unnumbered). Comment lines (%) are ignored.
Prints: number, env, label, first 90 chars of body. With extra args, prints only
rows whose label or number equals an arg (planted wrong numbers print NOTHING)."""
import re, sys
src = open(sys.argv[1], encoding='utf-8', errors='replace').read()
src = re.sub(r'(?m)(?<!\\)%.*$', '', src)
envs = {}   # env -> counter name
for m in re.finditer(r'\\newtheorem(\*?)\{(\w+)\}(?:\[(\w+)\])?\{[^}]*\}(?:\[(\w+)\])?', src):
    star, env, shared, within = m.groups()
    if star: continue
    envs[env] = shared if shared else env
# resolve chains
def root(e):
    seen = set()
    while e in envs and envs[e] != e and e not in seen:
        seen.add(e); e = envs[e]
    return e
pat = re.compile(r'\\(section)\*?\{|\\begin\{(\w+)\}(\[[^\]]*\])?(\\label\{[^}]*\})?')
sec = 0; cnt = 0; rows = []
for m in pat.finditer(src):
    if m.group(1):
        if src[m.end(0)-2:m.end(0)-1] != '*' and not src[m.start():m.end()].startswith('\\section*'):
            sec += 1; cnt = 0
        continue
    env = m.group(2)
    if env in envs:
        cnt += 1
        body = src[m.end():m.end()+400]
        lab = re.search(r'\\label\{([^}]*)\}', src[m.start():m.start()+200])
        rows.append((f'{sec}.{cnt}', env, lab.group(1) if lab else '-', ' '.join(body.split())[:90]))
want = sys.argv[2:]
for r in rows:
    if not want or r[0] in want or r[2] in want:
        print(*r, sep=' | ')
