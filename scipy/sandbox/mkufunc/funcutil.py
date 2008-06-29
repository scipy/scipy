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
        try:
            tmp = str(eval(name))
        except NameError:
            tmp = 'EVAL_FAILED'
        acc.write('%8s: %s\n' % (name, tmp))
    
    res = acc.getvalue()
    
    while True:
        tmp = pat_norep.sub('NO_REPRESENTATION', res)
        if tmp == res:
            break
        res = tmp
        
    return res


def func_hash(f, extra=None):
    txt = disassemble2(f.func_code) + repr(extra)
    #print txt
    txt = pat_white.sub(' ', txt)
    return hashlib.md5(txt).hexdigest()

if __name__ == '__main__':
#    import math
#    from math import *
    
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
        return math.sin(x)
    md5sums.append(func_hash(f))

    def f(x):
        return sin(x)
    md5sums.append(func_hash(f, float))
    
    print md5sums
    #assert md5sums == ['91d13599d610a554dccd6b44cb5ef1f0',
    #                   'be0c54b477180f897cbf7604fc565d18',
    #                   '732d1ef6c1ce8cc92a7f28917496d292']
