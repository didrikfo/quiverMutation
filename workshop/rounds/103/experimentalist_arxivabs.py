import sys,time,urllib.request,urllib.parse,xml.etree.ElementTree as ET
ns={'a':'http://www.w3.org/2005/Atom'}
args=sys.argv[1:]
for a in args:
    if a.startswith('q='): url='https://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':a[2:],'max_results':10})
    else: url='https://export.arxiv.org/api/query?id_list='+a
    time.sleep(6)
    for tr in range(4):
        try: d=urllib.request.urlopen(url,timeout=60).read(); break
        except Exception as e: print('retry',e); time.sleep(20); d=None
    if d is None: continue
    print('\n##',a)
    for e in ET.fromstring(d).findall('a:entry',ns):
        i=e.find('a:id',ns).text.split('/abs/')[1]
        t=' '.join(e.find('a:title',ns).text.split())
        s=' '.join(e.find('a:summary',ns).text.split())
        print(i,t); 
        if not a.startswith('q='): print('  ',s[:900])
