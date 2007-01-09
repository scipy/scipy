# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _sparsetools

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


PYA_QS_STACK = _sparsetools.PYA_QS_STACK
SMALL_QUICKSORT = _sparsetools.SMALL_QUICKSORT

def int_aquicksort(*args):
    """int_aquicksort(int v, unsigned int tosort, unsigned int num, void unused)"""
    return _sparsetools.int_aquicksort(*args)

def csrtocsc(*args):
    """
    csrtocsc(int n_row, int n_col, int Ap, int Aj, int Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bi, std::vector<(int)> Bx)
    csrtocsc(int n_row, int n_col, int Ap, int Aj, long Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bi, std::vector<(long)> Bx)
    csrtocsc(int n_row, int n_col, int Ap, int Aj, float Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bi, std::vector<(float)> Bx)
    csrtocsc(int n_row, int n_col, int Ap, int Aj, double Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bi, std::vector<(double)> Bx)
    csrtocsc(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(npy_cfloat)> Bx)
    csrtocsc(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(npy_cdouble)> Bx)
    """
    return _sparsetools.csrtocsc(*args)

def csctocsr(*args):
    """
    csctocsr(int n_row, int n_col, int Ap, int Ai, int Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bj, std::vector<(int)> Bx)
    csctocsr(int n_row, int n_col, int Ap, int Ai, long Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bj, std::vector<(long)> Bx)
    csctocsr(int n_row, int n_col, int Ap, int Ai, float Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bj, std::vector<(float)> Bx)
    csctocsr(int n_row, int n_col, int Ap, int Ai, double Ax, std::vector<(int)> Bp, 
        std::vector<(int)> Bj, std::vector<(double)> Bx)
    csctocsr(int n_row, int n_col, int Ap, int Ai, npy_cfloat Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(npy_cfloat)> Bx)
    csctocsr(int n_row, int n_col, int Ap, int Ai, npy_cdouble Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(npy_cdouble)> Bx)
    """
    return _sparsetools.csctocsr(*args)

def csrtocoo(*args):
    """
    csrtocoo(int n_row, int n_col, int Ap, int Aj, int Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(int)> Bx)
    csrtocoo(int n_row, int n_col, int Ap, int Aj, long Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(long)> Bx)
    csrtocoo(int n_row, int n_col, int Ap, int Aj, float Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(float)> Bx)
    csrtocoo(int n_row, int n_col, int Ap, int Aj, double Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(double)> Bx)
    csrtocoo(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        std::vector<(int)> Bi, std::vector<(int)> Bj, 
        std::vector<(npy_cfloat)> Bx)
    csrtocoo(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        std::vector<(int)> Bi, std::vector<(int)> Bj, 
        std::vector<(npy_cdouble)> Bx)
    """
    return _sparsetools.csrtocoo(*args)

def csctocoo(*args):
    """
    csctocoo(int n_row, int n_col, int Ap, int Ai, int Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(int)> Bx)
    csctocoo(int n_row, int n_col, int Ap, int Ai, long Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(long)> Bx)
    csctocoo(int n_row, int n_col, int Ap, int Ai, float Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(float)> Bx)
    csctocoo(int n_row, int n_col, int Ap, int Ai, double Ax, std::vector<(int)> Bi, 
        std::vector<(int)> Bj, std::vector<(double)> Bx)
    csctocoo(int n_row, int n_col, int Ap, int Ai, npy_cfloat Ax, 
        std::vector<(int)> Bi, std::vector<(int)> Bj, 
        std::vector<(npy_cfloat)> Bx)
    csctocoo(int n_row, int n_col, int Ap, int Ai, npy_cdouble Ax, 
        std::vector<(int)> Bi, std::vector<(int)> Bj, 
        std::vector<(npy_cdouble)> Bx)
    """
    return _sparsetools.csctocoo(*args)

def cootocsr(*args):
    """
    cootocsr(int n_row, int n_col, int NNZ, int Ai, int Aj, int Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(int)> Bx)
    cootocsr(int n_row, int n_col, int NNZ, int Ai, int Aj, long Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(long)> Bx)
    cootocsr(int n_row, int n_col, int NNZ, int Ai, int Aj, float Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(float)> Bx)
    cootocsr(int n_row, int n_col, int NNZ, int Ai, int Aj, double Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(double)> Bx)
    cootocsr(int n_row, int n_col, int NNZ, int Ai, int Aj, npy_cfloat Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(npy_cfloat)> Bx)
    cootocsr(int n_row, int n_col, int NNZ, int Ai, int Aj, npy_cdouble Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bj, 
        std::vector<(npy_cdouble)> Bx)
    """
    return _sparsetools.cootocsr(*args)

def cootocsc(*args):
    """
    cootocsc(int n_row, int n_col, int NNZ, int Ai, int Aj, int Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(int)> Bx)
    cootocsc(int n_row, int n_col, int NNZ, int Ai, int Aj, long Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(long)> Bx)
    cootocsc(int n_row, int n_col, int NNZ, int Ai, int Aj, float Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(float)> Bx)
    cootocsc(int n_row, int n_col, int NNZ, int Ai, int Aj, double Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(double)> Bx)
    cootocsc(int n_row, int n_col, int NNZ, int Ai, int Aj, npy_cfloat Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(npy_cfloat)> Bx)
    cootocsc(int n_row, int n_col, int NNZ, int Ai, int Aj, npy_cdouble Ax, 
        std::vector<(int)> Bp, std::vector<(int)> Bi, 
        std::vector<(npy_cdouble)> Bx)
    """
    return _sparsetools.cootocsc(*args)

def csrplcsr(*args):
    """
    csrplcsr(int n_row, int n_col, int Ap, int Aj, int Ax, int Bp, 
        int Bj, int Bx, std::vector<(int)> Cp, std::vector<(int)> Cj, 
        std::vector<(int)> Cx)
    csrplcsr(int n_row, int n_col, int Ap, int Aj, long Ax, int Bp, 
        int Bj, long Bx, std::vector<(int)> Cp, std::vector<(int)> Cj, 
        std::vector<(long)> Cx)
    csrplcsr(int n_row, int n_col, int Ap, int Aj, float Ax, int Bp, 
        int Bj, float Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(float)> Cx)
    csrplcsr(int n_row, int n_col, int Ap, int Aj, double Ax, int Bp, 
        int Bj, double Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(double)> Cx)
    csrplcsr(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        int Bp, int Bj, npy_cfloat Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(npy_cfloat)> Cx)
    csrplcsr(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        int Bp, int Bj, npy_cdouble Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(npy_cdouble)> Cx)
    """
    return _sparsetools.csrplcsr(*args)

def cscplcsc(*args):
    """
    cscplcsc(int n_row, int n_col, int Ap, int Ai, int Ax, int Bp, 
        int Bi, int Bx, std::vector<(int)> Cp, std::vector<(int)> Ci, 
        std::vector<(int)> Cx)
    cscplcsc(int n_row, int n_col, int Ap, int Ai, long Ax, int Bp, 
        int Bi, long Bx, std::vector<(int)> Cp, std::vector<(int)> Ci, 
        std::vector<(long)> Cx)
    cscplcsc(int n_row, int n_col, int Ap, int Ai, float Ax, int Bp, 
        int Bi, float Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(float)> Cx)
    cscplcsc(int n_row, int n_col, int Ap, int Ai, double Ax, int Bp, 
        int Bi, double Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(double)> Cx)
    cscplcsc(int n_row, int n_col, int Ap, int Ai, npy_cfloat Ax, 
        int Bp, int Bi, npy_cfloat Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(npy_cfloat)> Cx)
    cscplcsc(int n_row, int n_col, int Ap, int Ai, npy_cdouble Ax, 
        int Bp, int Bi, npy_cdouble Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(npy_cdouble)> Cx)
    """
    return _sparsetools.cscplcsc(*args)

def csrmucsr(*args):
    """
    csrmucsr(int n_row, int n_col, int Ap, int Aj, int Ax, int Bp, 
        int Bj, int Bx, std::vector<(int)> Cp, std::vector<(int)> Cj, 
        std::vector<(int)> Cx)
    csrmucsr(int n_row, int n_col, int Ap, int Aj, long Ax, int Bp, 
        int Bj, long Bx, std::vector<(int)> Cp, std::vector<(int)> Cj, 
        std::vector<(long)> Cx)
    csrmucsr(int n_row, int n_col, int Ap, int Aj, float Ax, int Bp, 
        int Bj, float Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(float)> Cx)
    csrmucsr(int n_row, int n_col, int Ap, int Aj, double Ax, int Bp, 
        int Bj, double Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(double)> Cx)
    csrmucsr(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        int Bp, int Bj, npy_cfloat Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(npy_cfloat)> Cx)
    csrmucsr(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        int Bp, int Bj, npy_cdouble Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(npy_cdouble)> Cx)
    """
    return _sparsetools.csrmucsr(*args)

def cscmucsc(*args):
    """
    cscmucsc(int n_row, int n_col, int Ap, int Ai, int Ax, int Bp, 
        int Bi, int Bx, std::vector<(int)> Cp, std::vector<(int)> Ci, 
        std::vector<(int)> Cx)
    cscmucsc(int n_row, int n_col, int Ap, int Ai, long Ax, int Bp, 
        int Bi, long Bx, std::vector<(int)> Cp, std::vector<(int)> Ci, 
        std::vector<(long)> Cx)
    cscmucsc(int n_row, int n_col, int Ap, int Ai, float Ax, int Bp, 
        int Bi, float Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(float)> Cx)
    cscmucsc(int n_row, int n_col, int Ap, int Ai, double Ax, int Bp, 
        int Bi, double Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(double)> Cx)
    cscmucsc(int n_row, int n_col, int Ap, int Ai, npy_cfloat Ax, 
        int Bp, int Bi, npy_cfloat Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(npy_cfloat)> Cx)
    cscmucsc(int n_row, int n_col, int Ap, int Ai, npy_cdouble Ax, 
        int Bp, int Bi, npy_cdouble Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(npy_cdouble)> Cx)
    """
    return _sparsetools.cscmucsc(*args)

def csrmux(*args):
    """
    csrmux(int n_row, int n_col, int Ap, int Aj, int Ax, int Xx, 
        std::vector<(int)> Yx)
    csrmux(int n_row, int n_col, int Ap, int Aj, long Ax, long Xx, 
        std::vector<(long)> Yx)
    csrmux(int n_row, int n_col, int Ap, int Aj, float Ax, float Xx, 
        std::vector<(float)> Yx)
    csrmux(int n_row, int n_col, int Ap, int Aj, double Ax, double Xx, 
        std::vector<(double)> Yx)
    csrmux(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        npy_cfloat Xx, std::vector<(npy_cfloat)> Yx)
    csrmux(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        npy_cdouble Xx, std::vector<(npy_cdouble)> Yx)
    """
    return _sparsetools.csrmux(*args)

def cscmux(*args):
    """
    cscmux(int n_row, int n_col, int Ap, int Ai, int Ax, int Xx, 
        std::vector<(int)> Yx)
    cscmux(int n_row, int n_col, int Ap, int Ai, long Ax, long Xx, 
        std::vector<(long)> Yx)
    cscmux(int n_row, int n_col, int Ap, int Ai, float Ax, float Xx, 
        std::vector<(float)> Yx)
    cscmux(int n_row, int n_col, int Ap, int Ai, double Ax, double Xx, 
        std::vector<(double)> Yx)
    cscmux(int n_row, int n_col, int Ap, int Ai, npy_cfloat Ax, 
        npy_cfloat Xx, std::vector<(npy_cfloat)> Yx)
    cscmux(int n_row, int n_col, int Ap, int Ai, npy_cdouble Ax, 
        npy_cdouble Xx, std::vector<(npy_cdouble)> Yx)
    """
    return _sparsetools.cscmux(*args)

def csrelmulcsr(*args):
    """
    csrelmulcsr(int n_row, int n_col, int Ap, int Aj, int Ax, int Bp, 
        int Bj, int Bx, std::vector<(int)> Cp, std::vector<(int)> Cj, 
        std::vector<(int)> Cx)
    csrelmulcsr(int n_row, int n_col, int Ap, int Aj, long Ax, int Bp, 
        int Bj, long Bx, std::vector<(int)> Cp, std::vector<(int)> Cj, 
        std::vector<(long)> Cx)
    csrelmulcsr(int n_row, int n_col, int Ap, int Aj, float Ax, int Bp, 
        int Bj, float Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(float)> Cx)
    csrelmulcsr(int n_row, int n_col, int Ap, int Aj, double Ax, int Bp, 
        int Bj, double Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(double)> Cx)
    csrelmulcsr(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        int Bp, int Bj, npy_cfloat Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(npy_cfloat)> Cx)
    csrelmulcsr(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        int Bp, int Bj, npy_cdouble Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Cj, std::vector<(npy_cdouble)> Cx)
    """
    return _sparsetools.csrelmulcsr(*args)

def cscelmulcsc(*args):
    """
    cscelmulcsc(int n_row, int n_col, int Ap, int Ai, int Ax, int Bp, 
        int Bi, int Bx, std::vector<(int)> Cp, std::vector<(int)> Ci, 
        std::vector<(int)> Cx)
    cscelmulcsc(int n_row, int n_col, int Ap, int Ai, long Ax, int Bp, 
        int Bi, long Bx, std::vector<(int)> Cp, std::vector<(int)> Ci, 
        std::vector<(long)> Cx)
    cscelmulcsc(int n_row, int n_col, int Ap, int Ai, float Ax, int Bp, 
        int Bi, float Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(float)> Cx)
    cscelmulcsc(int n_row, int n_col, int Ap, int Ai, double Ax, int Bp, 
        int Bi, double Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(double)> Cx)
    cscelmulcsc(int n_row, int n_col, int Ap, int Ai, npy_cfloat Ax, 
        int Bp, int Bi, npy_cfloat Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(npy_cfloat)> Cx)
    cscelmulcsc(int n_row, int n_col, int Ap, int Ai, npy_cdouble Ax, 
        int Bp, int Bi, npy_cdouble Bx, std::vector<(int)> Cp, 
        std::vector<(int)> Ci, std::vector<(npy_cdouble)> Cx)
    """
    return _sparsetools.cscelmulcsc(*args)

def spdiags(*args):
    """
    spdiags(int n_row, int n_col, int n_diag, int offsets, int diags, 
        std::vector<(int)> Ap, std::vector<(int)> Ai, 
        std::vector<(int)> Ax)
    spdiags(int n_row, int n_col, int n_diag, int offsets, long diags, 
        std::vector<(int)> Ap, std::vector<(int)> Ai, 
        std::vector<(long)> Ax)
    spdiags(int n_row, int n_col, int n_diag, int offsets, float diags, 
        std::vector<(int)> Ap, std::vector<(int)> Ai, 
        std::vector<(float)> Ax)
    spdiags(int n_row, int n_col, int n_diag, int offsets, double diags, 
        std::vector<(int)> Ap, std::vector<(int)> Ai, 
        std::vector<(double)> Ax)
    spdiags(int n_row, int n_col, int n_diag, int offsets, npy_cfloat diags, 
        std::vector<(int)> Ap, std::vector<(int)> Ai, 
        std::vector<(npy_cfloat)> Ax)
    spdiags(int n_row, int n_col, int n_diag, int offsets, npy_cdouble diags, 
        std::vector<(int)> Ap, std::vector<(int)> Ai, 
        std::vector<(npy_cdouble)> Ax)
    """
    return _sparsetools.spdiags(*args)

def csrtodense(*args):
    """
    csrtodense(int n_row, int n_col, int Ap, int Aj, int Ax, int Mx)
    csrtodense(int n_row, int n_col, int Ap, int Aj, long Ax, long Mx)
    csrtodense(int n_row, int n_col, int Ap, int Aj, float Ax, float Mx)
    csrtodense(int n_row, int n_col, int Ap, int Aj, double Ax, double Mx)
    csrtodense(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax, 
        npy_cfloat Mx)
    csrtodense(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax, 
        npy_cdouble Mx)
    """
    return _sparsetools.csrtodense(*args)

def densetocsr(*args):
    """
    densetocsr(int n_row, int n_col, int Mx, std::vector<(int)> Ap, 
        std::vector<(int)> Aj, std::vector<(int)> Ax)
    densetocsr(int n_row, int n_col, long Mx, std::vector<(int)> Ap, 
        std::vector<(int)> Aj, std::vector<(long)> Ax)
    densetocsr(int n_row, int n_col, float Mx, std::vector<(int)> Ap, 
        std::vector<(int)> Aj, std::vector<(float)> Ax)
    densetocsr(int n_row, int n_col, double Mx, std::vector<(int)> Ap, 
        std::vector<(int)> Aj, std::vector<(double)> Ax)
    densetocsr(int n_row, int n_col, npy_cfloat Mx, std::vector<(int)> Ap, 
        std::vector<(int)> Aj, std::vector<(npy_cfloat)> Ax)
    densetocsr(int n_row, int n_col, npy_cdouble Mx, std::vector<(int)> Ap, 
        std::vector<(int)> Aj, std::vector<(npy_cdouble)> Ax)
    """
    return _sparsetools.densetocsr(*args)

def ensure_sorted_indices(*args):
    """
    ensure_sorted_indices(int n_row, int n_col, int Ap, int Aj, int Ax)
    ensure_sorted_indices(int n_row, int n_col, int Ap, int Aj, long Ax)
    ensure_sorted_indices(int n_row, int n_col, int Ap, int Aj, float Ax)
    ensure_sorted_indices(int n_row, int n_col, int Ap, int Aj, double Ax)
    ensure_sorted_indices(int n_row, int n_col, int Ap, int Aj, npy_cfloat Ax)
    ensure_sorted_indices(int n_row, int n_col, int Ap, int Aj, npy_cdouble Ax)
    """
    return _sparsetools.ensure_sorted_indices(*args)

