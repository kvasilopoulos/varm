"""Minimal pure-Python MATLAB v5 .mat reader (no numpy/scipy). Usage:
   python3 matread.py file.mat            -> list top-level variables with class/dims
   python3 matread.py file.mat var [var]  -> write /<outdir>/<var>.tsv (flat column-major numeric) + print dims
   Library use: load(path) -> {name: node}; node = ('num', dims, values) | ('char', dims, str) | ('cell', dims, [nodes]) | ('struct', dims, {field: [nodes]}) | ('other', cls)
"""
import struct, zlib, sys, os
T = {1: ('b', 1), 2: ('B', 1), 3: ('h', 2), 4: ('H', 2), 5: ('i', 4), 6: ('I', 4), 7: ('f', 4), 9: ('d', 8), 12: ('q', 8), 13: ('Q', 8)}

def _elem(buf, pos):
    tag = struct.unpack_from('<I', buf, pos)[0]
    if tag >> 16:  # small data element
        typ, n = tag & 0xFFFF, tag >> 16
        return typ, buf[pos + 4:pos + 4 + n], pos + 8
    typ, n = struct.unpack_from('<II', buf, pos)
    data = buf[pos + 8:pos + 8 + n]
    return typ, data, pos + 8 + ((n + 7) // 8) * 8

def _nums(typ, data):
    if typ == 16 or typ == 17: return []
    c, sz = T[typ]
    return list(struct.unpack('<%d%s' % (len(data) // sz, c), data))

def _matrix(data):
    pos = 0
    if len(data) == 0: return '', ('num', [0, 0], [])
    _, flags, pos = _elem(data, pos)
    cls = struct.unpack_from('<I', flags, 0)[0] & 0xFF
    _, dimsb, pos = _elem(data, pos)
    dims = list(struct.unpack('<%di' % (len(dimsb) // 4), dimsb))
    _, nameb, pos = _elem(data, pos)
    name = nameb.decode('latin1')
    n = 1
    for d in dims: n *= d
    if cls == 1:
        cells = []
        for _ in range(n):
            t, d, pos = _elem(data, pos); cells.append(_matrix(d)[1])
        return name, ('cell', dims, cells)
    if cls in (2, 3):
        if cls == 3: _, _, pos = _elem(data, pos)  # class name
        _, lb, pos = _elem(data, pos); L = struct.unpack('<i', lb)[0]
        _, fb, pos = _elem(data, pos)
        fields = [fb[i:i + L].split(b'\0')[0].decode('latin1') for i in range(0, len(fb), L)]
        out = {f: [] for f in fields}
        for _ in range(n):
            for f in fields:
                t, d, pos = _elem(data, pos); out[f].append(_matrix(d)[1])
        return name, ('struct', dims, out)
    if cls == 4:
        t, d, pos = _elem(data, pos)
        s = d.decode('utf-16le' if t == 17 else 'latin1', 'ignore') if t in (16, 17) else ''.join(map(chr, _nums(t, d)))
        return name, ('char', dims, s)
    if cls in (6, 7, 8, 9, 10, 11, 12, 13, 14, 15):
        t, d, pos = _elem(data, pos)
        return name, ('num', dims, _nums(t, d))
    return name, ('other', cls)

def load(path):
    buf = open(path, 'rb').read(); pos = 128; out = {}
    while pos < len(buf):
        typ, n = struct.unpack_from('<II', buf, pos)
        if typ == 15:  # compressed elements are NOT padded to 8 bytes
            data = zlib.decompress(buf[pos + 8:pos + 8 + n]); pos += 8 + n
            typ, data, _ = _elem(data, 0)
        else:
            typ, data, pos = _elem(buf, pos)
        if typ == 14:
            n, v = _matrix(data); out[n] = v
    return out

def describe(node, depth=0, maxd=2):
    k = node[0]
    if k == 'struct' and depth < maxd: return 'struct%s{%s}' % (node[1], ','.join(node[2].keys()))
    if k == 'cell': return 'cell%s' % node[1]
    if k == 'num': return 'num%s' % node[1]
    if k == 'char': return 'char:%r' % node[2][:60]
    return str(k)

def tojson(n):
    k = n[0]
    if k == 'num': return {'dim': n[1], 'x': n[2]}
    if k == 'char': return {'str': n[2]}
    if k == 'cell': return {'dim': n[1], 'cell': [tojson(c) for c in n[2]]}
    if k == 'struct': return {'dim': n[1], 'struct': {f: [tojson(c) for c in v] for f, v in n[2].items()}}
    return {'other': n[1]}

if __name__ == '__main__' and len(sys.argv) > 3 and sys.argv[2] == '--json':
    import json
    m = load(sys.argv[1]); json.dump({v: tojson(m[v]) for v in sys.argv[4:]}, open(sys.argv[3], 'w'))
    sys.exit(0)

if __name__ == '__main__':
    m = load(sys.argv[1])
    if len(sys.argv) == 2:
        for k, v in m.items(): print(k, describe(v))
    else:
        for k in sys.argv[2:]:
            v = m[k]; print(k, describe(v))
            if v[0] == 'num':
                print(' '.join('%.12g' % x for x in v[2][:60]))
