# stdlib-only xlsx -> csv (first sheet). Usage: python3 -I xlsx2csv.py in.xlsx out.csv
import sys, zipfile, re, csv
import xml.etree.ElementTree as ET
ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
z = zipfile.ZipFile(sys.argv[1])
ss = []
if 'xl/sharedStrings.xml' in z.namelist():
    for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns):
        ss.append(''.join(t.text or '' for t in si.iter('{%s}t' % ns['m'])))
sh = sorted(n for n in z.namelist() if re.match(r'xl/worksheets/sheet\d+\.xml', n))[0]
rows = []
def col(ref):
    s = re.match(r'[A-Z]+', ref).group(); n = 0
    for c in s: n = n * 26 + ord(c) - 64
    return n - 1
for r in ET.fromstring(z.read(sh)).iter('{%s}row' % ns['m']):
    row = {}
    for c in r.findall('m:c', ns):
        v = c.find('m:v', ns)
        if v is None: continue
        row[col(c.get('r'))] = ss[int(v.text)] if c.get('t') == 's' else v.text
    rows.append(row)
w = max(max(r) for r in rows if r) + 1
with open(sys.argv[2], 'w', newline='') as f:
    cw = csv.writer(f)
    for r in rows: cw.writerow([r.get(i, '') for i in range(w)])
