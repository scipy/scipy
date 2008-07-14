#!/usr/bin/env python
"""
>>> mandel(-1, .3)
36
>>> mandel(0, 0)
-1
>>> mandel(10, 10)
1
>>> mandel(array([-1, 0, 10]), array([.3, 0, 10]))
array([36, -1,  1])
"""
import hashlib

from numpy import array
from scipy import weave

support_code = '''
#define D 1000

long iterations(double cr, double ci)
{
        long d = 1;
        double zr=cr, zi=ci, zr2, zi2;
        for(;;) {
                zr2 = zr * zr;
                zi2 = zi * zi;
                if( zr2+zi2 > 16.0 ) return d;
                if( ++d == D ) return -1;
                zi = 2.0 * zr * zi + ci;
                zr = zr2 - zi2 + cr;
        }
}

static void
PyUFunc_0(char **args, npy_intp *dimensions, npy_intp *steps, void *func)
{
        npy_intp i, n;
        npy_intp is0 = steps[0];
        npy_intp is1 = steps[1];
        npy_intp os = steps[2];
        char *ip0 = args[0];
        char *ip1 = args[1];
        char *op = args[2];
        n = dimensions[0];
        
        for(i = 0; i < n; i++) {
                *(long *)op = iterations(*(double *)ip0,
                                         *(double *)ip1);
                ip0 += is0;
                ip1 += is1;
                op += os;
        }
}

static PyUFuncGenericFunction f_functions[] = {
        PyUFunc_0,
};

static char f_types[] = {
        NPY_DOUBLE, NPY_DOUBLE, NPY_LONG,
};
'''
ufunc_info = weave.base_info.custom_info()
ufunc_info.add_header('"numpy/ufuncobject.h"')

mandel = weave.inline('/*' + hashlib.md5(support_code).hexdigest() + '''*/
import_ufunc();

return_val = PyUFunc_FromFuncAndData(
            f_functions,
            NULL,
            f_types,
            1,             /* ntypes */
            2,             /* nin */
            1,             /* nout */
            PyUFunc_None,  /* identity */
            "mandel",      /* name */
            "returns number of iterations from cr, ci",    /* doc */
            0);
            ''',
                     support_code=support_code,
                     verbose=0,
                     customize=ufunc_info)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
