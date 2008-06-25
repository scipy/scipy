import sys
import re
import cStringIO

import numpy
import scipy.weave as weave

from interactive import Translation

verbose = False
_cnt = 0

typedict = {
    int:    ['NPY_LONG',   'long'  ],
    long:   ['NPY_LONG',   'long'  ],
    float:  ['NPY_DOUBLE', 'double'],
}

class Cfunc(object):
    
    def __init__(self, f, signature):
        global _cnt
        _cnt += 1
        self.n = _cnt
        self.sig = signature
        self.nin = f.func_code.co_argcount     # input args
        self.nout = len(self.sig) - self.nin
        assert self.nout == 1                  # for now
        
        if not verbose:
            rem = sys.stderr
            sys.stderr = cStringIO.StringIO()
            
        t = Translation(f, backend='c')
        t.annotate(signature[:self.nin])
        t.source()
    
        if not verbose:
            sys.stderr = rem
    
        c_source_filename = t.driver.c_source_filename
        assert c_source_filename.endswith('.c')
        src = open(c_source_filename, 'r').read()
        
        self.prefix = 'f%i_' % self.n
        self.allCsrc = src.replace('pypy_', self.prefix + 'pypy_')
        self.cname = self.prefix + 'pypy_g_' + f.__name__
        
    def cfunc(self):
        p = re.compile(r'^\w+[*\s\w]+' + self.cname +
                       r'\s*\([^)]*\)\s*\{.*?[\n\r]\}[\n\r]',
                       re.DOTALL | re.MULTILINE | re.VERBOSE)
        
        found = p.findall(self.allCsrc)
        assert len(found) == 1
        res = found[0]
        res = res.replace(self.prefix + 'pypy_g_ll_math_ll_math_', '')
        return res + '\n'
    
    def decl(self):
        p = re.compile(r'^\w+[*\s\w]+' + self.cname +
                       r'\s*\([^)]*\);',
                       re.DOTALL | re.MULTILINE | re.VERBOSE)
        
        found = p.findall(self.allCsrc)
        assert len(found) == 1
        return found[0]


    def support_code(self):
        arg0type = typedict[self.sig[0]][1]
        rettype = typedict[self.sig[-1]][1]
        n = self.n
        cname = self.cname
        return '''
static %(rettype)s foo_%(n)i(%(arg0type)s x)
{
	return %(cname)s(x);
}

typedef %(rettype)s Func_%(n)i(%(arg0type)s);

static void
PyUFunc_%(n)i(char **args, npy_intp *dimensions, npy_intp *steps, void *func)
{
	/* printf("PyUFunc_%(n)i\\n"); */
	
	npy_intp n = dimensions[0];
	npy_intp is0 = steps[0];
	npy_intp os = steps[1];
	char *ip0 = args[0];
	char *op = args[1];
	Func_%(n)i *f = (Func_%(n)i *) func;
	npy_intp i;
	
	for(i = 0; i < n; i++, ip0 += is0, op += os) {
		%(arg0type)s *in1 = (%(arg0type)s *)ip0;
		%(rettype)s *out = (%(rettype)s *)op;
		
		*out = f(*in1);
	}
}
''' % locals()


def test1():    
    def sqr(x):
        return x * x
    #verbose = True
    for argtypes in ([int, int], [float, float]):
        x = Cfunc(sqr, argtypes)
        print x.cname, x.nin, x.nout, x.sig
        print x.cfunc()
        print '{{{%s}}}' % x.decl()
        print x.support_code()


def write_pypyc(cfuncs):
    fo = open('pypy.c', 'w');
    fo.write('#include "head.c"\n\n')
    for cf in cfuncs:
        fo.write(cf.cfunc())
    fo.close()


def genufunc(f, signatures):
    
    signatures.sort(key=lambda sig: [numpy.dtype(typ).num for typ in sig])
    
    cfuncs = [Cfunc(f, sig) for sig in signatures]
    
    write_pypyc(cfuncs)
    
    declarations = ''.join('\t%s\n' % cf.decl() for cf in cfuncs)

    func_support = ''.join(cf.support_code() for cf in cfuncs)

    pyufuncs = ''.join('\tPyUFunc_%i,\n' % cf.n for cf in cfuncs)
    
    data = ''.join('\t(void *) foo_%i,\n' % cf.n for cf in cfuncs)

    foo_signatures = ''.join('\t%s  /* %i */\n' %
                         (''.join(typedict[t][0] + ', ' for t in cf.sig), cf.n)
                         for cf in cfuncs)
    
    support_code = '''
extern "C" {
%(declarations)s}

%(func_support)s

static PyUFuncGenericFunction foo_functions[] = {
%(pyufuncs)s};

static void *foo_data[] = {
%(data)s};

static char foo_signatures[] = {
%(foo_signatures)s};
''' % locals()

    ntypes = len(signatures)
    nin = cfuncs[0].nin
    
    code = '''
import_ufunc();

return_val = PyUFunc_FromFuncAndData(
    foo_functions,
    foo_data,
    foo_signatures,
    %(ntypes)i,         /* ntypes */
    %(nin)i,            /* nin */
    1,                  /* nout */
    PyUFunc_None,       /* identity */
    "foo",              /* name */
    "",                 /* doc */
    0);
''' % locals()
    
    if 1:
        fo = open(f.__name__ + '_code.cc', 'w');
        fo.write(code);
        fo.close()
        
        fo = open(f.__name__ + '_support_code.cc', 'w');
        fo.write(support_code);
        fo.close()
    
    ufunc_info = weave.base_info.custom_info()
    ufunc_info.add_header('"numpy/ufuncobject.h"')
    ufunc_info.add_include_dir('"."')
    
    return weave.inline(code,
                        verbose=0, force=1, # XXX
                        support_code=support_code,
                        customize=ufunc_info,
                        sources=['pypy.c'])


def test2():

    def sqr(x):
        return x * x
    
    ufunc = genufunc(sqr, [
        (float, float),
        (int, int),
        ])
    
    x = array([0.0, 1.0, 2.5, 12.0])
    print "x =", x, x.dtype
    y = ufunc(x)
    print "y =", y, y.dtype
    
    x = array([0, 1, 2, 15])
    print "x =", x, x.dtype
    y = ufunc(x)
    print "y =", y, y.dtype


def mkufunc(signatures):
    print 'signatures', signatures
    
    class Compile(object):
        
        def __init__(self, f):
            self.ufunc = genufunc(f, signatures)

        def __call__(self, *args):
            return self.ufunc(*args)

    return Compile


if __name__ == '__main__':
    # test1();     exit()

    mkufunc([int, (float, int, float)])
    
