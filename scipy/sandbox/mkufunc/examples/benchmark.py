#!/usr/bin/env python
from math import sin, cos
import time, hashlib

from numpy import arange, vectorize, allclose, empty_like
from scipy import weave

from mkufunc.api import mkufunc


def f(x):
    return 4.2 * x * x + 3.7 * x + 1.5


vfunc = vectorize(f)

mfunc = mkufunc([(float, float)])(f)

#####################################################################
support_code = '''
static void
PyUFunc_0(char **args, npy_intp *dimensions, npy_intp *steps, void *func)
{
        npy_intp i, n;
        npy_intp is0 = steps[0];
        npy_intp os = steps[1];
        char *ip0 = args[0];
        char *op = args[1];
        n = dimensions[0];
        double x;
        
        for(i = 0; i < n; i++) {
                x = *(double *)ip0;
                double *out = (double *)op;
                
                *out = 4.2 * x*x + 3.7 * x + 1.5;

                ip0 += is0;
                op += os;
        }
}

static PyUFuncGenericFunction f_functions[] = {
        PyUFunc_0,
};

static char f_types[] = {
        NPY_DOUBLE, NPY_DOUBLE,
};
'''
ufunc_info = weave.base_info.custom_info()
ufunc_info.add_header('"numpy/ufuncobject.h"')

ufunc = weave.inline('/*' + hashlib.md5(support_code).hexdigest() + '''*/
import_ufunc();

return_val = PyUFunc_FromFuncAndData(
            f_functions,
            NULL,
            f_types,
            1,             /* ntypes */
            1,             /* nin */
            1,             /* nout */
            PyUFunc_None,  /* identity */
            "f",           /* name */
            "doc",         /* doc */
            0);
            ''',
                     support_code=support_code,
                     verbose=0,
                     customize=ufunc_info)
#############################################################


x = arange(0, 1000, 0.0001)    #print "x =", x, x.dtype

start_time = time.time()
b_y = empty_like(x)
weave.blitz("b_y[:] = 4.2 * x[:] * x[:] + 3.7 * x[:] + 1.5")
b_time = time.time() - start_time
print 'blitz: %.6f sec' % b_time

start_time = time.time()
n_y = f(x)
n_time = time.time() - start_time
print 'numpy: %.6f sec' % n_time

start_time = time.time()
v_y = vfunc(x)
v_time = time.time() - start_time
print 'vectorize: %.6f sec' % v_time

start_time = time.time()
m_y = mfunc(x)
m_time = time.time() - start_time
print 'mkufunc: %.6f sec' % m_time

start_time = time.time()
u_y = ufunc(x)
u_time = time.time() - start_time
print 'ufunc: %.6f sec' % u_time

print "speedup over blitz:",     b_time/m_time
print "speedup over numpy:",     n_time/m_time
print "speedup over vectorize:", v_time/m_time
print "speedup over ufunc:", u_time/m_time

assert allclose(b_y, n_y)
assert allclose(v_y, n_y)
assert allclose(m_y, n_y)
assert allclose(u_y, n_y)
