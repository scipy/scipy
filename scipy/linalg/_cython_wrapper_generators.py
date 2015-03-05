"""
Code generator script to make the Cython BLAS and LAPACK wrappers
from the files "cython_blas_signatures.txt" and
"cython_lapack_signatures.txt" which contain the signatures for
all the BLAS/LAPACK routines that should be included in the wrappers.
"""

from operator import itemgetter

fortran_types = {'int':'integer',
                 'c':'complex',
                 'd':'double precision',
                 's':'real',
                 'z':'complex*16',
                 'char':'character',
                 'bint':'logical'}

c_types = {'int':'int',
           'c':'_scipy_linalg_float_complex',
           'd':'double',
           's':'float',
           'z':'_scipy_linalg_double_complex',
           'char':'char',
           'bint':'int',
           'cselect1':'_cselect1',
           'cselect2':'_cselect2',
           'dselect2':'_dselect2',
           'dselect3':'_dselect3',
           'sselect2':'_sselect2',
           'sselect3':'_sselect3',
           'zselect1':'_zselect1',
           'zselect2':'_zselect2'}

def arg_names_and_types(args):
    return zip(*[arg.split(' *') for arg in args.split(', ')])

pyx_func_template = """ctypedef void _w{name}_t({ret_type} *out, {args}) nogil
cdef extern from "{header_name}":
    void _fortran_{name} "F_FUNC({name}wrp, {upname}WRP)"({ret_type} *out, {args}) nogil
cdef {ret_type} _wrap_{name}({args}) nogil:
    cdef {ret_type} out
    _fortran_{name}(&out, {argnames})
    return out
cdef {name}_t *{name}_f = &_wrap_{name}
"""

pyx_func_template = """
cdef extern from "{header_name}":
    void _fortran_{name} "F_FUNC({name}wrp, {upname}WRP)"({ret_type} *out, {args}) nogil
cdef {ret_type} {name}({args}) nogil:
    cdef {ret_type} out
    _fortran_{name}(&out, {argnames})
    return out
"""

def pyx_decl_func(name, ret_type, args, header_name):
    argtypes, argnames = arg_names_and_types(args)
    # Fix the case where one of the arguments has the same name as the
    # abbreviation for the argument type.
    # Otherwise the variable passed as an argument is considered overwrites
    # the previous typedef and Cython compilation fails.
    if ret_type in argnames:
        argnames = [n if n != ret_type else ret_type + '_' for n in argnames]
        args = ', '.join([' *'.join([n, t]) for n, t in zip(argtypes, argnames)])
    argnames = ', '.join(argnames)
    c_ret_type = c_types[ret_type]
    return pyx_func_template.format(name=name, upname=name.upper(), args=args,
                                    ret_type=ret_type, c_ret_type=c_ret_type,
                                    argnames=argnames, header_name=header_name)

pyx_sub_template = """cdef extern from "{header_name}":
    void _fortran_{name} "F_FUNC({name},{upname})"({args}) nogil
cdef {name}_t *{name}_f = &_fortran_{name}
"""

pyx_sub_template = """cdef extern from "{header_name}":
    void _fortran_{name} "F_FUNC({name},{upname})"({args}) nogil
cdef void {name}({args}) nogil:
    _fortran_{name}({argnames})
"""

def pyx_decl_sub(name, args, header_name):
    argtypes, argnames = arg_names_and_types(args)
    argnames = ', '.join(argnames)
    return pyx_sub_template.format(name=name, upname=name.upper(),
                                   args=args, argnames=argnames,
                                   header_name=header_name)

blas_pyx_preamble = '''# cython: boundscheck = False
# cython: wraparound = False
# cython: cdivision = True

"""
BLAS Functions for Cython
=========================

Usable from Cython via::

    cimport scipy.linalg.cython_blas

Raw function pointers (Fortran-style pointer arguments):

- {}


"""

cdef extern from "fortran_defs.h":
    pass

'''

def make_blas_pyx_preamble(all_sigs):
    names = [sig[0] for sig in all_sigs]
    return blas_pyx_preamble.format("\n- ".join(names))

lapack_pyx_preamble = '''"""
LAPACK functions for Cython
===========================

Usable from Cython via::

    cimport scipy.linalg.cython_lapack

This module provides Cython-level wrappers for all primary routines included
in LAPACK 3.1.0 and some of the fixed-api auxiliary routines.

The signature for dcgesv changed from LAPACK 3.1.1 to LAPACK 3.2.0.
The version here is the newer of the two since it matches the signature
from later versions of LAPACK and the version in the CLAPACK included in OSX.

Raw function pointers (Fortran-style pointer arguments):

- {}


"""

cdef extern from "fortran_defs.h":
    pass

'''

def make_lapack_pyx_preamble(all_sigs):
    names = [sig[0] for sig in all_sigs]
    return lapack_pyx_preamble.format("\n- ".join(names))

blas_py_wrappers = """

# Python-accessible wrappers for testing:

cdef inline bint _is_contiguous(double[:,:] a, int axis) nogil:
    return (a.strides[axis] == sizeof(a[0,0]) or a.shape[axis] == 1)

cpdef float complex _test_cdotc(float complex[:] cx, float complex[:] cy) nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
        int incy = cy.strides[0] // sizeof(cy[0])
    return cdotc(&n, &cx[0], &incx, &cy[0], &incy)

cpdef float complex _test_cdotu(float complex[:] cx, float complex[:] cy) nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
        int incy = cy.strides[0] // sizeof(cy[0])
    return cdotu(&n, &cx[0], &incx, &cy[0], &incy)

cpdef double _test_dasum(double[:] dx) nogil:
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
    return dasum(&n, &dx[0], &incx)

cpdef double _test_ddot(double[:] dx, double[:] dy) nogil:
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
        int incy = dy.strides[0] // sizeof(dy[0])
    return ddot(&n, &dx[0], &incx, &dy[0], &incy)

cpdef int _test_dgemm(double alpha, double[:,:] a, double[:,:] b, double beta,
                double[:,:] c) nogil except -1:
    cdef:
        char *transa
        char *transb
        int m, n, k, lda, ldb, ldc
        double *a0=&a[0,0]
        double *b0=&b[0,0]
        double *c0=&c[0,0]
    # In the case that c is C contiguous, swap a and b and
    # swap whether or not each of them is transposed.
    # This can be done because a.dot(b) = b.T.dot(a.T).T.
    if _is_contiguous(c, 1):
        if _is_contiguous(a, 1):
            transb = 'n'
            ldb = (&a[1,0]) - a0 if a.shape[0] > 1 else 1
        elif _is_contiguous(a, 0):
            transb = 't'
            ldb = (&a[0,1]) - a0 if a.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'a' is neither C nor Fortran contiguous.")
        if _is_contiguous(b, 1):
            transa = 'n'
            lda = (&b[1,0]) - b0 if b.shape[0] > 1 else 1
        elif _is_contiguous(b, 0):
            transa = 't'
            lda = (&b[0,1]) - b0 if b.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'b' is neither C nor Fortran contiguous.")
        k = b.shape[0]
        if k != a.shape[1]:
            with gil:
                raise ValueError("Shape mismatch in input arrays.")
        m = b.shape[1]
        n = a.shape[0]
        if n != c.shape[0] or m != c.shape[1]:
            with gil:
                raise ValueError("Output array does not have the correct shape.")
        ldc = (&c[1,0]) - c0 if c.shape[0] > 1 else 1
        dgemm(transa, transb, &m, &n, &k, &alpha, b0, &lda, a0,
                   &ldb, &beta, c0, &ldc)
    elif _is_contiguous(c, 0):
        if _is_contiguous(a, 1):
            transa = 't'
            lda = (&a[1,0]) - a0 if a.shape[0] > 1 else 1
        elif _is_contiguous(a, 0):
            transa = 'n'
            lda = (&a[0,1]) - a0 if a.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'a' is neither C nor Fortran contiguous.")
        if _is_contiguous(b, 1):
            transb = 't'
            ldb = (&b[1,0]) - b0 if b.shape[0] > 1 else 1
        elif _is_contiguous(b, 0):
            transb = 'n'
            ldb = (&b[0,1]) - b0 if b.shape[1] > 1 else 1
        else:
            with gil:
                raise ValueError("Input 'b' is neither C nor Fortran contiguous.")
        m = a.shape[0]
        k = a.shape[1]
        if k != b.shape[0]:
            with gil:
                raise ValueError("Shape mismatch in input arrays.")
        n = b.shape[1]
        if m != c.shape[0] or n != c.shape[1]:
            with gil:
                raise ValueError("Output array does not have the correct shape.")
        ldc = (&c[0,1]) - c0 if c.shape[1] > 1 else 1
        dgemm(transa, transb, &m, &n, &k, &alpha, a0, &lda, b0,
                   &ldb, &beta, c0, &ldc)
    else:
        with gil:
            raise ValueError("Input 'c' is neither C nor Fortran contiguous.")
    return 0

cpdef double _test_dnrm2(double[:] x) nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return dnrm2(&n, &x[0], &incx)

cpdef double _test_dzasum(double complex[:] zx) nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
    return dzasum(&n, &zx[0], &incx)

cpdef double _test_dznrm2(double complex[:] x) nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return dznrm2(&n, &x[0], &incx)

cpdef int _test_icamax(float complex[:] cx) nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
    return icamax(&n, &cx[0], &incx)

cpdef int _test_idamax(double[:] dx) nogil:
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
    return idamax(&n, &dx[0], &incx)

cpdef int _test_isamax(float[:] sx) nogil:
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
    return isamax(&n, &sx[0], &incx)

cpdef int _test_izamax(double complex[:] zx) nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
    return izamax(&n, &zx[0], &incx)

cpdef float _test_sasum(float[:] sx) nogil:
    cdef:
        int n = sx.shape[0]
        int incx = sx.shape[0] // sizeof(sx[0])
    return sasum(&n, &sx[0], &incx)

cpdef float _test_scasum(float complex[:] cx) nogil:
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
    return scasum(&n, &cx[0], &incx)

cpdef float _test_scnrm2(float complex[:] x) nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return scnrm2(&n, &x[0], &incx)

cpdef float _test_sdot(float[:] sx, float[:] sy) nogil:
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
        int incy = sy.strides[0] // sizeof(sy[0])
    return sdot(&n, &sx[0], &incx, &sy[0], &incy)

cpdef float _test_snrm2(float[:] x) nogil:
    cdef:
        int n = x.shape[0]
        int incx = x.shape[0] // sizeof(x[0])
    return snrm2(&n, &x[0], &incx)

cpdef double complex _test_zdotc(double complex[:] zx, double complex[:] zy) nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
        int incy = zy.strides[0] // sizeof(zy[0])
    return zdotc(&n, &zx[0], &incx, &zy[0], &incy)

cpdef double complex _test_zdotu(double complex[:] zx, double complex[:] zy) nogil:
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
        int incy = zy.strides[0] // sizeof(zy[0])
    return zdotu(&n, &zx[0], &incx, &zy[0], &incy)
"""

def generate_blas_pyx(func_sigs, sub_sigs, all_sigs, header_name):
    funcs = "\n".join(pyx_decl_func(*(s+(header_name,))) for s in func_sigs)
    subs = "\n" + "\n".join(pyx_decl_sub(*(s[::2]+(header_name,)))
                            for s in sub_sigs)
    return make_blas_pyx_preamble(all_sigs) + funcs + subs + blas_py_wrappers

lapack_py_wrappers = """

# Python accessible wrappers for testing:

def _test_dlamch(cmach):
    # This conversion is necessary to handle Python 3 strings.
    cmach_bytes = bytes(cmach)
    # Now that it is a bytes representation, a non-temporary variable
    # must be passed as a part of the function call.
    cdef char* cmach_char = cmach_bytes
    return dlamch(cmach_char)

def _test_slamch(cmach):
    # This conversion is necessary to handle Python 3 strings.
    cmach_bytes = bytes(cmach)
    # Now that it is a bytes representation, a non-temporary variable
    # must be passed as a part of the function call.
    cdef char* cmach_char = cmach_bytes
    return slamch(cmach_char)
"""

def generate_lapack_pyx(func_sigs, sub_sigs, all_sigs, header_name):
    funcs = "\n".join(pyx_decl_func(*(s+(header_name,))) for s in func_sigs)
    subs = "\n" + "\n".join(pyx_decl_sub(*(s[::2]+(header_name,)))
                            for s in sub_sigs)
    preamble = make_lapack_pyx_preamble(all_sigs)
    return preamble + funcs + subs + lapack_py_wrappers

pxd_template = """ctypedef {ret_type} {name}_t({args}) nogil
cdef {name}_t *{name}_f
"""
pxd_template = """cdef {ret_type} {name}({args}) nogil
"""

def pxd_decl(name, ret_type, args):
    return pxd_template.format(name=name, ret_type=ret_type, args=args)

blas_pxd_preamble = """ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

"""

def generate_blas_pxd(all_sigs):
    body = '\n'.join(pxd_decl(*sig) for sig in all_sigs)
    return blas_pxd_preamble + body

lapack_pxd_preamble = """ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

# Function pointer type declarations for
# gees and gges families of functions.
ctypedef bint cselect1(c*)
ctypedef bint cselect2(c*, c*)
ctypedef bint dselect2(d*, d*)
ctypedef bint dselect3(d*, d*, d*)
ctypedef bint sselect2(s*, s*)
ctypedef bint sselect3(s*, s*, s*)
ctypedef bint zselect1(z*)
ctypedef bint zselect2(z*, z*)

"""

def generate_lapack_pxd(all_sigs):
    return lapack_pxd_preamble + '\n'.join(pxd_decl(*sig) for sig in all_sigs)

fortran_template = """      subroutine {name}wrp(ret, {argnames})
        external {wrapper}
        {ret_type} {wrapper}
        {ret_type} ret
        {argdecls}
        ret = {wrapper}({argnames})
      end
"""

def process_fortran_name(name):
    if 'inc' in name:
        return name
    if 'x' in name or 'y' in name:
        return name + '(n)'
    return name

def fort_subroutine_wrapper(name, ret_type, args):
    if name[0] in ['c', 's'] or name in ['zladiv', 'zdotu', 'zdotc']:
        wrapper = 'w' + name
    else:
        wrapper = name
    types, names = arg_names_and_types(args)
    argnames = ', '.join(names)
    
    names = [process_fortran_name(n) for n in names]
    argdecls = '\n        '.join('{} {}'.format(fortran_types[t], n)
                                 for n, t in zip(names, types))
    return fortran_template.format(name=name, wrapper=wrapper,
                                   argnames=argnames, argdecls=argdecls,
                                   ret_type=fortran_types[ret_type])

def generate_fortran(func_sigs):
    return "\n".join(fort_subroutine_wrapper(*sig) for sig in func_sigs)

def make_c_args(args):
    types, names = arg_names_and_types(args)
    types = [c_types[arg] for arg in types]
    return ', '.join('{} *{}'.format(t, n) for t, n in zip(types, names))

c_func_template = "void F_FUNC({name}wrp, {upname}WRP)({return_type} *ret, {args});\n"

def c_func_decl(name, return_type, args):
    args = make_c_args(args)
    return_type = c_types[return_type]
    return c_func_template.format(name=name, upname=name.upper(),
                                  return_type=return_type, args=args)

c_sub_template = "void F_FUNC({name},{upname})({args});\n"

def c_sub_decl(name, return_type, args):
    args = make_c_args(args)
    return c_sub_template.format(name=name, upname=name.upper(), args=args)

c_preamble = """#ifndef SCIPY_LINALG_{lib}_FORTRAN_WRAPPERS_H
#define SCIPY_LINALG_{lib}_FORTRAN_WRAPPERS_H
#include "fortran_defs.h"
"""

# Mimic the complex declarations from cython
complex_decls = """
// Mimic the complex number definitions from cython.
#if !defined(_SCIPY_LINALG_CCOMPLEX)
  #if defined(__cplusplus)
    #define _SCIPY_LINALG_CCOMPLEX 1
  #elif defined(_Complex_I)
    #define _SCIPY_LINALG_CCOMPLEX 1
  #else
    #define _SCIPY_LINALG_CCOMPLEX 0
  #endif
#endif

#if _SCIPY_LINALG_CCOMPLEX
  #ifdef __cplusplus
    #include <complex>
  #else
    #include <complex.h>
  #endif
#endif

#if _SCIPY_LINALG_CCOMPLEX && !defined(__cplusplus) && defined(__sun__) && defined(__GNUC__)
  #undef _Complex_I
  #define _Complex_I 1.0fj
#endif

#if _SCIPY_LINALG_CCOMPLEX
  #ifdef __cplusplus
    typedef ::std::complex< float > _scipy_linalg_float_complex;
  #else
    typedef float _Complex _scipy_linalg_float_complex;
  #endif
#else
    typedef struct { float real, imag; } _scipy_linalg_float_complex;
#endif

#if _SCIPY_LINALG_CCOMPLEX
  #ifdef __cplusplus
    typedef ::std::complex< double > _scipy_linalg_double_complex;
  #else
    typedef double _Complex _scipy_linalg_double_complex;
  #endif
#else
    typedef struct { double real, imag; } _scipy_linalg_double_complex;
#endif
"""

lapack_decls = """
typedef int (*_cselect1)(_scipy_linalg_float_complex*);
typedef int (*_cselect2)(_scipy_linalg_float_complex*, _scipy_linalg_float_complex*);
typedef int (*_dselect2)(double*, double*);
typedef int (*_dselect3)(double*, double*, double*);
typedef int (*_sselect2)(float*, float*);
typedef int (*_sselect3)(float*, float*, float*);
typedef int (*_zselect1)(_scipy_linalg_double_complex*);
typedef int (*_zselect2)(_scipy_linalg_double_complex*, _scipy_linalg_double_complex*);
"""

cpp_guard = """
#ifdef __cplusplus
extern "C" {
#endif

"""

c_end = """
#ifdef __cplusplus
}
#endif
#endif
"""

def generate_c_header(func_sigs, sub_sigs, all_sigs, lib_name):
    funcs = "".join(c_func_decl(*sig) for sig in func_sigs)
    subs = "\n" + "".join(c_sub_decl(*sig) for sig in sub_sigs)
    if lib_name == 'LAPACK':
        preamble = (c_preamble.format(lib=lib_name) +
                    complex_decls + lapack_decls)
    else:
        preamble = c_preamble.format(lib=lib_name) + complex_decls
    return "".join([preamble, cpp_guard, funcs, subs, c_end])

def split_signature(sig):
    name_and_type, args = sig[:-1].split('(')
    ret_type, name = name_and_type.split(' ')
    return name, ret_type, args

def filter_lines(ls):
    ls = [l.strip() for l in ls if l != '\n' and l[0] != '#']
    func_sigs = [split_signature(l) for l in ls if l.split(' ')[0] != 'void']
    sub_sigs = [split_signature(l) for l in ls if l.split(' ')[0] == 'void']
    all_sigs = list(sorted(func_sigs + sub_sigs, key=itemgetter(0)))
    return func_sigs, sub_sigs, all_sigs

def make_all(blas_signature_file="cython_blas_signatures.txt",
             lapack_signature_file="cython_lapack_signatures.txt",
             blas_name="cython_blas",
             lapack_name="cython_lapack",
             blas_fortran_name="_blas_subroutine_wrappers.f",
             lapack_fortran_name="_lapack_subroutine_wrappers.f",
             blas_header_name="_blas_subroutines.h",
             lapack_header_name="_lapack_subroutines.h"):
    with open(blas_signature_file, 'r') as f:
        blas_sigs = f.readlines()
    blas_sigs = filter_lines(blas_sigs)
    blas_pyx = generate_blas_pyx(*(blas_sigs + (blas_header_name,)))
    with open(blas_name + '.pyx', 'w') as f:
        f.write(blas_pyx)
    blas_pxd = generate_blas_pxd(blas_sigs[2])
    with open(blas_name + '.pxd', 'w') as f:
        f.write(blas_pxd)
    blas_fortran = generate_fortran(blas_sigs[0])
    with open(blas_fortran_name, 'w') as f:
        f.write(blas_fortran)
    blas_c_header = generate_c_header(*(blas_sigs + ('BLAS',)))
    with open(blas_header_name, 'w') as f:
        f.write(blas_c_header)
    with open(lapack_signature_file, 'r') as f:
        lapack_sigs = f.readlines()
    lapack_sigs = filter_lines(lapack_sigs)
    lapack_pyx = generate_lapack_pyx(*(lapack_sigs + (lapack_header_name,)))
    with open(lapack_name + '.pyx', 'w') as f:
        f.write(lapack_pyx)
    lapack_pxd = generate_lapack_pxd(lapack_sigs[2])
    with open(lapack_name + '.pxd', 'w') as f:
        f.write(lapack_pxd)
    lapack_fortran = generate_fortran(lapack_sigs[0])
    with open(lapack_fortran_name, 'w') as f:
        f.write(lapack_fortran)
    lapack_c_header = generate_c_header(*(lapack_sigs + ('LAPACK',)))
    with open(lapack_header_name, 'w') as f:
        f.write(lapack_c_header)

if __name__ == '__main__':
    make_all()
