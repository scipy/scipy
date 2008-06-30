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

def dis2(co):
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


def func_hash(f, salt=None, verbose=0):
    """ Return the MD5 hash for a function object as string.

    'salt' can be any object that has a representation
    """
    txt = dis2(f.func_code) + repr(salt)
    if verbose:
        print txt
    
    txt = pat_white.sub(' ', txt)
    return hashlib.md5(txt).hexdigest()


if __name__ == '__main__':
    pass
