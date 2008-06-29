import sys, re, dis, hashlib, cStringIO


def disassemble(co):
    """ Given a code object, return output from dis.disassemble as a string.

    (dis.disassemble writes its output to stdout.)
    """
    tmp = sys.stdout
    sys.stdout = cStringIO.StringIO()
    dis.disassemble(co)
    res = sys.stdout.getvalue()
    sys.stdout = tmp
    return res

pat_norep = re.compile(r'<[^<>]*>')
pat_white = re.compile(r'\s+')

def disassemble2(co):
    acc = cStringIO.StringIO()
    for line in disassemble(co).splitlines():
        line = line[16:].strip()
        if line:
            acc.write(line+'\n')
    
    acc.write('co_names:\n')
    for name in co.co_names:
        if name in '''math exp log sqrt
                      cos sin tan acos asin atan atan2'''.split():
            continue
        acc.write('%8s: %s\n' % (name, eval(name)))
    
    res = acc.getvalue()
    
    while True:
        tmp = pat_norep.sub('NO_REPRESENTATION', res)
        if tmp == res:
            break
        res = tmp
        
    return res

def func_hash(f):
    txt = disassemble2(f.func_code)
    #print txt
    txt = pat_white.sub(' ', txt)
    return hashlib.md5(txt).hexdigest()


if __name__ == '__main__':
    import math
    from math import *
    
    md5sums = []
    
    b = 3.14159
    g = lambda x: x
    def h(n):
        return n + 3
    
    for a in xrange(2):
        def f(x):
            inner1 = lambda t: t/3.0
            def inner2(): return
            t = b + g(42) + h(4)
            return sin(pi * x) + a + t
        md5sums.append(func_hash(f))
   

    def f(x):
        return math.sin(x) + math.cos(x)
    md5sums.append(func_hash(f))

    assert md5sums == ['b821514915e98426c49d93f58e400025',
                       '2bf13d8983c80c8fd773db4534a2c1b6',
                       '8d2ce5ab9152dabc3d49d0732fb84666']
