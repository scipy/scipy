""" mkufunc (make U function)


Author: Ilan Schnell (with help from Travis Oliphant and Eric Jones)
"""
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
    """ C compiled python functions

    >>> def sqr(x):
    ...     return x * x

    >>> signature = [int, int] # only the input arguments are used here
    
    compilation is done upon initialization
    >>> x = Cfunc(sqr, signature)
    >>> x.nin # number of input arguments
    1
    >>> x.nout # number of output arguments (must be 1 for now)
    1
    >>> x.sig
    [<type 'int'>, <type 'int'>]

    Attributes:

        n           -- id number
        sig         -- signature
        nin         -- number of input arguments
        nout        -- number of output arguments
        cname       -- name of the C function

    Methods:

        decl()      -- returns the C declaration for the function
        cfunc()     -- returns the C function (as string)
        ufunc_support_code()
                    -- generate the C support code to make this
                       function part work with PyUFuncGenericFunction
    """
    def __init__(self, f, signature):
        global _cnt
        _cnt += 1
        self.n = _cnt
        self.sig = signature
        self.nin = f.func_code.co_argcount     # input args
        self.nout = len(self.sig) - self.nin
        assert self.nout == 1                  # for now
        
        if not verbose:
            tmp = sys.stderr
            sys.stderr = cStringIO.StringIO()
            
        t = Translation(f, backend='c')
        t.annotate(signature[:self.nin])
        t.source()
    
        if not verbose:
            sys.stderr = tmp
        
        c_source_filename = t.driver.c_source_filename
        assert c_source_filename.endswith('.c')
        src = open(c_source_filename, 'r').read()
        
        self._prefix = 'f%i_' % self.n
        self._allCsrc = src.replace('pypy_', self._prefix + 'pypy_')
        self.cname = self._prefix + 'pypy_g_' + f.__name__
        
    def cfunc(self):
        p = re.compile(r'^\w+[*\s\w]+' + self.cname +
                       r'\s*\([^)]*\)\s*\{.*?[\n\r]\}[\n\r]',
                       re.DOTALL | re.MULTILINE | re.VERBOSE)
        
        found = p.findall(self._allCsrc)
        assert len(found) == 1
        res = found[0]
        res = res.replace(self._prefix + 'pypy_g_ll_math_ll_math_', '')
        return res + '\n'
    
    def decl(self):
        p = re.compile(r'^\w+[*\s\w]+' + self.cname +
                       r'\s*\([^)]*\);',
                       re.DOTALL | re.MULTILINE | re.VERBOSE)
        
        found = p.findall(self._allCsrc)
        assert len(found) == 1
        return found[0]


    def ufunc_support_code(self):
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
    """ Given a list of Cfunc instances, write the C code containing the
    functions into a file.
    """
    fo = open('pypy.c', 'w');
    fo.write('#include "head.c"\n\n')
    for cf in cfuncs:
        fo.write(cf.cfunc())
    fo.close()


def genufunc(f, signatures):
    """ Given a Python function and its signatures, do the following:

    - Compile the function to C for each signature

    - Write the C code for all these functions to a file

    - Generate the support code for weave

    - Generate the code for weave.  This contains the actual call to
      PyUFuncGenericFunction

    - Return the Ufunc Python object
    """
    signatures.sort(key=lambda sig: [numpy.dtype(typ).num for typ in sig])
    
    cfuncs = [Cfunc(f, sig) for sig in signatures]
    
    write_pypyc(cfuncs)
    
    declarations = ''.join('\t%s\n' % cf.decl() for cf in cfuncs)

    func_support = ''.join(cf.ufunc_support_code() for cf in cfuncs)

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
    """ The actual API function, to be used as decorator function.
        
    """
    #print 'signatures', signatures
    
    class Compile(object):
        
        def __init__(self, f):
            self.ufunc = genufunc(f, signatures)

        def __call__(self, *args):
            return self.ufunc(*args)

    return Compile


if __name__ == '__main__':
    import doctest
    doctest.testmod()
