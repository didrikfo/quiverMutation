"""Round 114 reply: are the malformed lines of E-263's h=4 hit files rows with an entry >= 10 written as two digits, and
are the 44 'lost' class rows exactly those? Compares the tail-K>=2 class (3800 rows, own enumerator) with the stored files.
usage: experimentalist_hitsvs.py DIR"""
import sys, glob, os
os.environ['TAIL'] = '1'
sys.argv = [sys.argv[0], 'xx'] + sys.argv[1:]
D = sys.argv[2]
src = open('workshop/rounds/114/experimentalist_tail.py').read().split("mode = sys.argv[1]")[0]
exec(src)
n, K1 = 15, 2
tg = {key(n, lone(n, 4)), key(n, lone(n, 7))}
cls = [s + (0,)*K1 for s in rows(n-K1) if key(n, s + (0,)*K1) in tg]
print('class rows (tail K>=2)', len(cls))
lines = []
for f in glob.glob(D + '/experimentalist_scan15_hits_4_*.txt'):
    lines += [l.strip() for l in open(f) if l.strip()]
S = set(lines); bad = {l for l in lines if len(l) != 13}
wide = [r for r in cls if max(r) >= 10]
inS = [r for r in cls if ''.join(map(str, r)) in S]
print('class rows with an entry >= 10:', len(wide), '; class rows whose joined string is in the files:', len(inS), '; absent:', len(cls)-len(inS))
print('absent rows == rows with an entry >= 10 :', {r for r in cls if ''.join(map(str, r)) not in S} == set(wide))
print('malformed lines that are the joined string of a class row:', len({''.join(map(str, r)) for r in wide} & bad), 'of', len(bad))
