import time, sys, urllib.parse, urllib.request, re, xml.etree.ElementTree as ET, ssl
Q = [
 'all:Nakayama AND all:"derived equivalence" AND all:mutation',
 'all:"linear Nakayama" AND all:"derived equivalence"',
 'all:Nakayama AND all:"derived equivalence classes"',
 'all:"tilting mutation" AND all:Nakayama',
 'all:"silting mutation" AND all:Nakayama',
 'all:"derived equivalence classification" AND all:gentle AND all:"number of"',
 'all:"derived equivalence" AND all:Nakayama AND all:Catalan',
 'all:"derived equivalence" AND all:"Nakayama algebras" AND all:"vertex"',
 'all:"mutation" AND all:"quivers with relations" AND all:"derived equivalence"',
 'all:"derived equivalence class" AND all:"gentle algebras" AND all:"counting"',
 'all:"Nakayama algebras" AND all:"Catalan"',
 'all:"tilting" AND all:"Nakayama algebras" AND all:"orbit"',
 'all:"Coxeter polynomial" AND all:Nakayama',
 'all:"tilting graph" AND all:Nakayama',
]
ns={'a':'http://www.w3.org/2005/Atom'}
seen={}
for q in Q:
    url='https://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':q,'max_results':12,'sortby':'relevance'})
    d=None
    for tr in range(4):
        time.sleep(6)
        try: d=urllib.request.urlopen(url,timeout=60).read(); break
        except Exception as e: print('retry',e); time.sleep(20)
    if d is None: print('ERR',q); continue
    r=ET.fromstring(d)
    print('\nQ:',q)
    for e in r.findall('a:entry',ns):
        i=e.find('a:id',ns).text.split('/abs/')[1]
        t=' '.join(e.find('a:title',ns).text.split())
        new = i not in seen; seen[i]=t
        print(('  * ' if new else '    '),i,t[:100])
