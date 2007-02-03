"""MA: a facility for dealing with missing observations
MA is generally used as a numpy.array look-alike.
by Paul F. Dubois.

Copyright 1999, 2000, 2001 Regents of the University of California.
Released for unlimited redistribution.
Adapted for numpy_core 2005 by Travis Oliphant and
(mainly) Paul Dubois.

Subclassing of the base ndarray 2006 by Pierre Gerard-Marchant.
pgmdevlist_at_gmail_dot_com

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

__all__ = ['MAError', 'MaskType', 'MaskedArray',
           'bool_', 'complex_', 'float_', 'int_', 'object_',
           'abs', 'absolute', 'add', 'all', 'allclose', 'allequal', 'alltrue',
               'amax', 'amin', 'anom', 'anomalies', 'any', 'arange', 
               'arccos', 'arccosh', 'arcsin', 'arcsinh', 'arctan', 'arctan2', 
               'arctanh', 'argmax', 'argmin', 'argsort', 'around', 
               'array', 'asarray',  
           'bitwise_and', 'bitwise_or', 'bitwise_xor',
           'ceil', 'choose', 'compressed', 'concatenate', 'conjugate', 
               'cos', 'cosh', 'count',
           'diagonal', 'divide', 'dump', 'dumps',
           'empty', 'empty_like', 'equal', 'exp', 
           'fabs', 'fmod', 'filled', 'floor', 'floor_divide', 
           'getmask', 'getmaskarray', 'greater', 'greater_equal', 'hypot', 
           'ids', 'inner', 'innerproduct',
               'isMA', 'isMaskedArray', 'is_mask', 'is_masked', 'isarray',
           'left_shift', 'less', 'less_equal', 'load', 'loads', 'log', 'log10', 
               'logical_and', 'logical_not', 'logical_or', 'logical_xor',
           'make_mask', 'make_mask_none', 'mask_or', 'masked', 
               'masked_array', 'masked_equal', 'masked_greater',
               'masked_greater_equal', 'masked_inside', 'masked_less', 
               'masked_less_equal', 'masked_not_equal', 'masked_object', 
               'masked_outside', 'masked_print_option', 'masked_singleton',
               'masked_values', 'masked_where', 'max', 'maximum', 'mean', 'min',
               'minimum', 'multiply', 
           'negative', 'nomask', 'nonzero', 'not_equal',
           'ones', 'outer', 'outerproduct', 
           'power', 'product', 'ptp', 'put', 'putmask', 
           'rank', 'ravel', 'remainder', 'repeat', 'reshape', 'resize', 
               'right_shift', 'round_', 
           'shape', 'sin', 'sinh', 'size', 'sometrue', 'sort', 'sqrt', 'std', 
               'subtract', 'sum', 'swapaxes',
           'take', 'tan', 'tanh', 'transpose', 'true_divide', 
           'var', 'where',
           'zeros']

import sys
import types
import cPickle
#
import numpy
from numpy import bool_, complex_, float_, int_, object_, str_

import numpy.core.umath as umath 
import numpy.core.fromnumeric  as fromnumeric
from numpy.core.numeric import ndarray
from numpy.core.fromnumeric import amax, amin
import numpy.core.numerictypes as ntypes
from numpy.core.numerictypes import bool_, typecodes
from numpy.core.multiarray import dtype
import numpy.core.numeric as numeric
from numpy.lib.shape_base import expand_dims as n_expand_dims
import warnings



MaskType = bool_
nomask = MaskType(0)

divide_tolerance = 1.e-35
numpy.seterr(all='ignore')

#####--------------------------------------------------------------------------
#---- --- Helper functions ---
#####--------------------------------------------------------------------------
def convert_typecode(f,dtchar):
    """Converts the type of `f` to a type compatible with `dtchar`, for inline operations."""
    ftype = f.dtype.char
    if dtchar == ftype:
        return f
    elif dtchar in typecodes['Integer']:
        if ftype in typecodes['Integer']:
            f = f.astype(dtchar)
        else:
            raise TypeError, 'Incorrect type for in-place operation.'
    elif dtchar in typecodes['Float']:
        if ftype in typecodes['Integer']:
            f = f.astype(dtchar)
        elif ftype in typecodes['Float']:
            f = f.astype(dtchar)
        else:
            raise TypeError, 'Incorrect type for in-place operation.'
    elif dtchar in typecodes['Complex']:
        if ftype in typecodes['Integer']:
            f = f.astype(dtchar)
        elif ftype in typecodes['Float']:
            f = f.astype(dtchar)
        elif ftype in typecodes['Complex']:
            f = f.astype(dtchar)
        else:
            raise TypeError, 'Incorrect type for in-place operation.'
    else:
        raise TypeError, 'Incorrect type for in-place operation.'
    return f

#####--------------------------------------------------------------------------
#---- --- Exceptions ---
#####--------------------------------------------------------------------------
class MAError(Exception):
    "Class for MA related errors."
    def __init__ (self, args=None):
        "Creates an exception."
        Exception.__init__(self,args)
        self.args = args
    def __str__(self):
        "Calculates the string representation."
        return str(self.args)
    __repr__ = __str__

#####--------------------------------------------------------------------------
#---- --- Filling options ---
#####--------------------------------------------------------------------------
# b: boolean - c: complex - f: floats - i: integer - O: object - S: string
default_filler = {'b': True,
                  'c' : 1.e20 + 0.0j,
                  'f' : 1.e20,
                  'i' : 999999,
                  'O' : '?',
                  'S' : 'N/A',
                  'u' : 999999,
                  'V' : '???',        
                  }
#{0: <type 'numpy.bool_'>,
# 1: <type 'numpy.int8'>,
# 2: <type 'numpy.uint8'>,
# 3: <type 'numpy.int16'>,
# 4: <type 'numpy.uint16'>,
# 5: <type 'numpy.int32'>,
# 6: <type 'numpy.uint32'>,
# 7: <type 'numpy.int64'>,
# 8: <type 'numpy.uint64'>,
# 9: <type 'numpy.int64'>,
# 10: <type 'numpy.uint64'>,
# 11: <type 'numpy.float32'>,
# 12: <type 'numpy.float64'>,
# 13: <type 'numpy.float128'>,
# 14: <type 'numpy.complex64'>,
# 15: <type 'numpy.complex128'>,
# 16: <type 'numpy.complex256'>,
# 17: <type 'numpy.object_'>,
# 18: <type 'numpy.string_'>,
# 19: <type 'numpy.unicode_'>,
# 20: <type 'numpy.void'>,
max_filler = ntypes._minvals   
max_filler.update([(k,-numeric.inf) for k in [numpy.float32, numpy.float64]])         
min_filler = ntypes._maxvals
min_filler.update([(k,numeric.inf) for k in [numpy.float32, numpy.float64]])    
if 'float128' in ntypes.typeDict:
    max_filler.update([(numpy.float128,-numeric.inf)])
    min_filler.update([(numpy.float128, numeric.inf)])       


def default_fill_value (obj):
    "Calculates the default fill value for an object `obj`."
    if hasattr(obj,'dtype'):
        return default_filler[obj.dtype.kind]
    elif isinstance(obj, float):
        return default_filler['f']
    elif isinstance(obj, int) or isinstance(obj, long):
        return default_filler['i']
    elif isinstance(obj, str):
        return default_filler['S']
    elif isinstance(obj, complex):
        return default_filler['c']
    elif isinstance(obj, numeric.dtype):
        return default_filler[obj.kind]
    else:
        return default_filler['O']

def minimum_fill_value (obj):
    "Calculates the default fill value suitable for taking the minimum of `obj`."
    if hasattr(obj, 'dtype'):
        objtype = obj.dtype
        filler = min_filler[objtype]
        if filler is None:
            raise TypeError, 'Unsuitable type for calculating minimum.'
        return filler
    elif isinstance(obj, float):
        return min_filler[ntypes.typeDict['float_']]
    elif isinstance(obj, int):
        return min_filler[ntypes.typeDict['int_']]
    elif isinstance(obj, long):
        return min_filler[ntypes.typeDict['uint']]
    elif isinstance(obj, numeric.dtype):
        return min_filler[obj]
    else:
        raise TypeError, 'Unsuitable type for calculating minimum.'

def maximum_fill_value (obj):
    "Calculates the default fill value suitable for taking the maximum of `obj`."
    if hasattr(obj, 'dtype'):
        objtype = obj.dtype
        filler = max_filler[objtype]
        if filler is None:
            raise TypeError, 'Unsuitable type for calculating minimum.'
        return filler
    elif isinstance(obj, float):
        return max_filler[ntypes.typeDict['float_']]
    elif isinstance(obj, int):
        return max_filler[ntypes.typeDict['int_']]
    elif isinstance(obj, long):
        return max_filler[ntypes.typeDict['uint']]
    elif isinstance(obj, numeric.dtype):
        return max_filler[obj]
    else:
        raise TypeError, 'Unsuitable type for calculating minimum.'

def set_fill_value (a, fill_value):
    "Sets the fill value of `a` if it is a masked array."
    if isinstance(a, MaskedArray):
        a.set_fill_value(fill_value)

def get_fill_value (a):
    """Returns the fill value of `a`, if any.
    Otherwise, returns the default fill value for that type.
    """
    if isinstance(a, MaskedArray):
        result = a.fill_value
    else:
        result = default_fill_value(a)
    return result

def common_fill_value (a, b):
    "Returns the common fill_value of `a` and `b`, if any, or `None`."
    t1 = get_fill_value(a)
    t2 = get_fill_value(b)
    if t1 == t2: 
        return t1
    return None

#................................................
def filled(a, value = None):
    """Returns `a` as an array with masked data replaced by `value`.
If `value` is `None` or the special element `masked`, `get_fill_value(a)`
is used instead.

If `a` is already a contiguous numeric array, `a` itself is returned.

`filled(a)` can be used to be sure that the result is numeric when passing 
an object a to other software ignorant of MA, in particular to numpy itself.
    """
    if hasattr(a, 'filled'):
        return a.filled(value)
    elif isinstance(a, ndarray): # and a.flags['CONTIGUOUS']:
        return a
    elif isinstance(a, types.DictType):
        return numeric.array(a, 'O')
    else:
        return numeric.array(a)
        
#####--------------------------------------------------------------------------
#---- --- Ufuncs ---
#####--------------------------------------------------------------------------
ufunc_domain = {}
ufunc_fills = {}

class domain_check_interval:
    """Defines a valid interval, 
so that `domain_check_interval(a,b)(x) = true` where `x < a` or `x > b`."""
    def __init__(self, a, b):
        "domain_check_interval(a,b)(x) = true where x < a or y > b"
        if (a > b):
            (a, b) = (b, a)
        self.a = a
        self.b = b

    def __call__ (self, x):
        "Execute the call behavior."
        return umath.logical_or(umath.greater (x, self.b),
                                umath.less(x, self.a))
#............................
class domain_tan:
    """Defines a valid interval for the `tan` function,
so that `domain_tan(eps) = True where `abs(cos(x)) < eps`"""
    def __init__(self, eps):
        "domain_tan(eps) = true where abs(cos(x)) < eps)"
        self.eps = eps
    def __call__ (self, x):
        "Execute the call behavior."
        return umath.less(umath.absolute(umath.cos(x)), self.eps)
#............................
class domain_safe_divide:
    """defines a domain for safe division."""
    def __init__ (self, tolerance=divide_tolerance):
        self.tolerance = tolerance
    def __call__ (self, a, b):
        return umath.absolute(a) * self.tolerance >= umath.absolute(b)
#............................
class domain_greater:
    "domain_greater(v)(x) = true where x <= v"
    def __init__(self, critical_value):
        "domain_greater(v)(x) = true where x <= v"
        self.critical_value = critical_value
        
    def __call__ (self, x):
        "Execute the call behavior."
        return umath.less_equal(x, self.critical_value)
#............................
class domain_greater_equal:
    "domain_greater_equal(v)(x) = true where x < v"
    def __init__(self, critical_value):
        "domain_greater_equal(v)(x) = true where x < v"
        self.critical_value = critical_value

    def __call__ (self, x):
        "Execute the call behavior."
        return umath.less(x, self.critical_value)
#..............................................................................
class masked_unary_operation:
    """Defines masked version of unary operations,
where invalid values are pre-masked. 

:IVariables:
    - `f` : function.
    - `fill` : Default filling value *[0]*.
    - `domain` : Default domain *[None]*. 
    """
    def __init__ (self, mufunc, fill=0, domain=None):
        """ masked_unary_operation(aufunc, fill=0, domain=None)
            aufunc(fill) must be defined
            self(x) returns aufunc(x)
            with masked values where domain(x) is true or getmask(x) is true.
        """
        self.f = mufunc
        self.fill = fill
        self.domain = domain
        self.__doc__ = getattr(mufunc, "__doc__", str(mufunc))
        self.__name__ = getattr(mufunc, "__name__", str(mufunc))
        ufunc_domain[mufunc] = domain
        ufunc_fills[mufunc] = fill
    #
    def __call__ (self, a, *args, **kwargs):
        "Execute the call behavior."
# numeric tries to return scalars rather than arrays when given scalars.
        m = getmask(a)
        d1 = filled(a, self.fill)
        if self.domain is not None:
            m = mask_or(m, numeric.asarray(self.domain(d1)))
        result = self.f(d1, *args, **kwargs)
        #
        if isinstance(result, MaskedArray):
            return result.__class__(result, mask=m)
        return masked_array(result, mask=m)
    #
    def __str__ (self):
        return "Masked version of %s. [Invalid values are masked]" % str(self.f)
#..............................................................................
class masked_binary_operation:
    """Defines masked version of binary operations,
where invalid values are pre-masked. 

:IVariables:
    - `f` : function.
    - `fillx` : Default filling value for first array*[0]*.
    - `filly` : Default filling value for second array*[0]*.
    - `domain` : Default domain *[None]*. 
    """
    def __init__ (self, mbfunc, fillx=0, filly=0):
        """abfunc(fillx, filly) must be defined.
           abfunc(x, filly) = x for all x to enable reduce.
        """
        self.f = mbfunc
        self.fillx = fillx
        self.filly = filly
        self.__doc__ = getattr(mbfunc, "__doc__", str(mbfunc))
        self.__name__ = getattr(mbfunc, "__name__", str(mbfunc))
        ufunc_domain[mbfunc] = None
        ufunc_fills[mbfunc] = (fillx, filly)
    #
    def __call__ (self, a, b, *args, **kwargs):
        "Execute the call behavior."
        m = mask_or(getmask(a), getmask(b))
        d1 = filled(a, self.fillx)
        d2 = filled(b, self.filly)
        result = self.f(d1, d2, *args, **kwargs)
#        if isinstance(result, ndarray) \
#               and m.ndim != 0 \
#               and m.shape != result.shape:
#            m = mask_or(getmaskarray(a), getmaskarray(b))
        if isinstance(result, MaskedArray):
            return result.__class__(result, mask=m)
        return masked_array(result, mask=m)
    #
    def reduce (self, target, axis=0, dtype=None):
        """Reduces `target` along the given `axis`."""
        if isinstance(target, MaskedArray):
            tclass = target.__class__
        else:
            tclass = MaskedArray
        m = getmask(target)
        t = filled(target, self.filly)
        if t.shape == ():
            t = t.reshape(1)
            if m is not nomask:
                m = make_mask(m, copy=1)
                m.shape = (1,)
        if m is nomask:
            return tclass(self.f.reduce (t, axis))
        else:
            t = tclass(t, mask=m)
            # XXX: "or t.dtype" below is a workaround for what appears
            # XXX: to be a bug in reduce.
            t = self.f.reduce(filled(t, self.filly), axis, dtype=dtype or t.dtype)
            m = umath.logical_and.reduce(m, axis)
            if isinstance(t, ndarray):
                return tclass(t, mask=m, fill_value=get_fill_value(target))
            elif m:
                return masked
            else:
                return t

    def outer (self, a, b):
        "Returns the function applied to the outer product of a and b."
        ma = getmask(a)
        mb = getmask(b)
        if ma is nomask and mb is nomask:
            m = nomask
        else:
            ma = getmaskarray(a)
            mb = getmaskarray(b)
            m = umath.logical_or.outer(ma, mb)
        d = self.f.outer(filled(a, self.fillx), filled(b, self.filly))
        if isinstance(d, MaskedArray):
            return d.__class__(d, mask=m)
        return masked_array(d, mask=m)

    def accumulate (self, target, axis=0):
        """Accumulates `target` along `axis` after filling with y fill value."""
        if isinstance(target, MaskedArray):
            tclass = target.__class__
        else:
            tclass = masked_array
        t = filled(target, self.filly)
        return tclass(self.f.accumulate(t, axis))
    
    def __str__ (self):
        return "Masked version of " + str(self.f)
#..............................................................................
class domained_binary_operation:
    """Defines binary operations that have a domain, like divide. 
    
These are complicated so they are a separate class. 
They have no reduce, outer or accumulate.

:IVariables:
    - `f` : function.
    - `fillx` : Default filling value for first array*[0]*.
    - `filly` : Default filling value for second array*[0]*.
    - `domain` : Default domain *[None]*. 
    """
    def __init__ (self, dbfunc, domain, fillx=0, filly=0):
        """abfunc(fillx, filly) must be defined.
           abfunc(x, filly) = x for all x to enable reduce.
        """
        self.f = dbfunc
        self.domain = domain
        self.fillx = fillx
        self.filly = filly
        self.__doc__ = getattr(dbfunc, "__doc__", str(dbfunc))
        self.__name__ = getattr(dbfunc, "__name__", str(dbfunc))
        ufunc_domain[dbfunc] = domain
        ufunc_fills[dbfunc] = (fillx, filly)

    def __call__(self, a, b):
        "Execute the call behavior."
        ma = getmask(a)
        mb = getmask(b)
        d1 = filled(a, self.fillx)
        d2 = filled(b, self.filly)
        t = numeric.asarray(self.domain(d1, d2))

        if fromnumeric.sometrue(t, None):
            d2 = numeric.where(t, self.filly, d2)
            mb = mask_or(mb, t)
        m = mask_or(ma, mb)
        result =  self.f(d1, d2)
        if isinstance(result, MaskedArray):
            return result.__class__(result, mask=m) 
        return masked_array(result, mask=m)

    def __str__ (self):
        return "Masked version of " + str(self.f)

#..............................................................................
# Unary ufuncs
exp = masked_unary_operation(umath.exp)
conjugate = masked_unary_operation(umath.conjugate)
sin = masked_unary_operation(umath.sin)
cos = masked_unary_operation(umath.cos)
tan = masked_unary_operation(umath.tan)
arctan = masked_unary_operation(umath.arctan)
arcsinh = masked_unary_operation(umath.arcsinh)
sinh = masked_unary_operation(umath.sinh)
cosh = masked_unary_operation(umath.cosh)
tanh = masked_unary_operation(umath.tanh)
abs = absolute = masked_unary_operation(umath.absolute)
fabs = masked_unary_operation(umath.fabs)
negative = masked_unary_operation(umath.negative)
floor = masked_unary_operation(umath.floor)
ceil = masked_unary_operation(umath.ceil)
around = masked_unary_operation(fromnumeric.round_)
logical_not = masked_unary_operation(umath.logical_not)
# Domained unary ufuncs
sqrt = masked_unary_operation(umath.sqrt, 0.0, domain_greater_equal(0.0))
log = masked_unary_operation(umath.log, 1.0, domain_greater(0.0))
log10 = masked_unary_operation(umath.log10, 1.0, domain_greater(0.0))
tan = masked_unary_operation(umath.tan, 0.0, domain_tan(1.e-35))
arcsin = masked_unary_operation(umath.arcsin, 0.0, 
                                domain_check_interval(-1.0, 1.0))
arccos = masked_unary_operation(umath.arccos, 0.0, 
                                domain_check_interval(-1.0, 1.0))
arccosh = masked_unary_operation(umath.arccosh, 1.0, domain_greater_equal(1.0))
arctanh = masked_unary_operation(umath.arctanh, 0.0, 
                                 domain_check_interval(-1.0+1e-15, 1.0-1e-15))
# Binary ufuncs
add = masked_binary_operation(umath.add)
subtract = masked_binary_operation(umath.subtract)
multiply = masked_binary_operation(umath.multiply, 1, 1)
arctan2 = masked_binary_operation(umath.arctan2, 0.0, 1.0)
equal = masked_binary_operation(umath.equal)
equal.reduce = None
not_equal = masked_binary_operation(umath.not_equal)
not_equal.reduce = None
less_equal = masked_binary_operation(umath.less_equal)
less_equal.reduce = None
greater_equal = masked_binary_operation(umath.greater_equal)
greater_equal.reduce = None
less = masked_binary_operation(umath.less)
less.reduce = None
greater = masked_binary_operation(umath.greater)
greater.reduce = None
logical_and = masked_binary_operation(umath.logical_and)
alltrue = masked_binary_operation(umath.logical_and, 1, 1).reduce
logical_or = masked_binary_operation(umath.logical_or)
sometrue = logical_or.reduce
logical_xor = masked_binary_operation(umath.logical_xor)
bitwise_and = masked_binary_operation(umath.bitwise_and)
bitwise_or = masked_binary_operation(umath.bitwise_or)
bitwise_xor = masked_binary_operation(umath.bitwise_xor)
hypot = masked_binary_operation(umath.hypot)
# Domained binary ufuncs
divide = domained_binary_operation(umath.divide, domain_safe_divide(), 0, 1)
true_divide = domained_binary_operation(umath.true_divide, 
                                        domain_safe_divide(), 0, 1)
floor_divide = domained_binary_operation(umath.floor_divide, 
                                         domain_safe_divide(), 0, 1)
remainder = domained_binary_operation(umath.remainder, 
                                      domain_safe_divide(), 0, 1)
fmod = domained_binary_operation(umath.fmod, domain_safe_divide(), 0, 1)


#####--------------------------------------------------------------------------
#---- --- Mask creation functions ---
#####--------------------------------------------------------------------------
def getmask(a):
    """Returns the mask of `a`, if any, or `nomask`.
Returns `nomask` if `a` is not a masked array.
To get an array for sure use getmaskarray."""
    if hasattr(a, "_mask"):
        return a._mask
    else:
        return nomask

def getmaskarray(a):
    """Returns the mask of `a`, if any.
Otherwise, returns an array of `False`, with the same shape as `a`.
    """
    m = getmask(a)
    if m is nomask:
        return make_mask_none(fromnumeric.shape(a))
    else:
        return m

def is_mask(m):
    """Returns `True` if `m` is a legal mask.
Does not check contents, only type.
    """
    try:
        return m.dtype.type is MaskType
    except AttributeError:
        return False
#
def make_mask(m, copy=False, small_mask=True, flag=None):
    """make_mask(m, copy=0, small_mask=0)
Returns `m` as a mask, creating a copy if necessary or requested.
The function can accept any sequence of integers or `nomask`. 
Does not check that contents must be 0s and 1s.
If `small_mask=True`, returns `nomask` if `m` contains no true elements.

:Parameters:
    - `m` (ndarray) : Mask.
    - `copy` (boolean, *[False]*) : Returns a copy of `m` if true.
    - `small_mask` (boolean, *[False]*): Flattens mask to `nomask` if `m` is all false.
    """
    if flag is not None:
        warnings.warn("The flag 'flag' is now called 'small_mask'!",
                      DeprecationWarning)
        small_mask = flag
    if m is nomask:
        return nomask
    elif isinstance(m, ndarray):
        m = filled(m, True)
        if m.dtype.type is MaskType:
            if copy:
                result = numeric.array(m, dtype=MaskType, copy=copy)
            else:
                result = m
        else:
            result = numeric.array(m, dtype=MaskType)
    else:
        result = numeric.array(filled(m, True), dtype=MaskType)
    # Bas les masques !
    if small_mask and not result.any():
        return nomask
    else:
        return result

def make_mask_none(s):
    "Returns a mask of shape `s`, filled with `False`."
    result = numeric.zeros(s, dtype=MaskType)
    return result

def mask_or (m1, m2, copy=False, small_mask=True):
    """Returns the combination of two masks `m1` and `m2`.
The masks are combined with the `logical_or` operator, treating `nomask` as false.
The result may equal m1 or m2 if the other is nomask.

:Parameters:
    - `m` (ndarray) : Mask.
    - `copy` (boolean, *[False]*) : Returns a copy of `m` if true.
    - `small_mask` (boolean, *[False]*): Flattens mask to `nomask` if `m` is all false.
     """
    if m1 is nomask: 
        return make_mask(m2, copy=copy, small_mask=small_mask)
    if m2 is nomask: 
        return make_mask(m1, copy=copy, small_mask=small_mask)
    if m1 is m2 and is_mask(m1): 
        return m1
    return make_mask(umath.logical_or(m1, m2), copy=copy, small_mask=small_mask)

#####--------------------------------------------------------------------------
#--- --- Masking functions ---
#####--------------------------------------------------------------------------
def masked_where(condition, x, copy=True):
    """Returns `x` as an array masked where `condition` is true.
Masked values of `x` or `condition` are kept.

:Parameters:
    - `condition` (ndarray) : Masking condition.
    - `x` (ndarray) : Array to mask.
    - `copy` (boolean, *[False]*) : Returns a copy of `m` if true.
    """
    cm = filled(condition,1)
    if isinstance(x,MaskedArray):
        m = mask_or(x._mask, cm)
        return x.__class__(x._data, mask=m, copy=copy)
    else:
        return MaskedArray(fromnumeric.asarray(x), copy=copy, mask=cm)

def masked_greater(x, value, copy=1):
    "Shortcut to `masked_where`, with ``condition = (x > value)``."
    return masked_where(greater(x, value), x, copy=copy)

def masked_greater_equal(x, value, copy=1):
    "Shortcut to `masked_where`, with ``condition = (x >= value)``."
    return masked_where(greater_equal(x, value), x, copy=copy)

def masked_less(x, value, copy=True):
    "Shortcut to `masked_where`, with ``condition = (x < value)``."
    return masked_where(less(x, value), x, copy=copy)

def masked_less_equal(x, value, copy=True):
    "Shortcut to `masked_where`, with ``condition = (x <= value)``."
    return masked_where(less_equal(x, value), x, copy=copy)

def masked_not_equal(x, value, copy=True):
    "Shortcut to `masked_where`, with ``condition = (x != value)``."
    return masked_where((x != value), x, copy=copy)

#
def masked_equal(x, value, copy=True):
    """Shortcut to `masked_where`, with ``condition = (x == value)``.
For floating point, consider `masked_values(x, value)` instead.
    """
    return masked_where((x == value), x, copy=copy)
#    d = filled(x, 0)
#    c = umath.equal(d, value)
#    m = mask_or(c, getmask(x))
#    return array(d, mask=m, copy=copy)

def masked_inside(x, v1, v2, copy=True):
    """Shortcut to `masked_where`, where `condition` is True for x inside 
the interval `[v1,v2]` ``(v1 <= x <= v2)``.
The boundaries `v1` and `v2` can be given in either order.
    """
    if v2 < v1:
        (v1, v2) = (v2, v1)
    xf = filled(x)
    condition = (xf >= v1) & (xf <= v2)
    return masked_where(condition, x, copy=copy)

def masked_outside(x, v1, v2, copy=True):
    """Shortcut to `masked_where`, where `condition` is True for x outside 
the interval `[v1,v2]` ``(x < v1)|(x > v2)``.
The boundaries `v1` and `v2` can be given in either order.
    """
    if v2 < v1:
        (v1, v2) = (v2, v1)
    xf = filled(x)
    condition = (xf < v1) | (xf > v2)
    return masked_where(condition, x, copy=copy)

#
def masked_object(x, value, copy=True):
    """Masks the array `x` where the data are exactly equal to `value`.
This function is suitable only for `object` arrays: for floating point, 
please use `masked_values` instead.
The mask is set to `nomask` if posible.

:parameter copy (Boolean, *[True]*):  Returns a copy of `x` if true. """
    if isMaskedArray(x):
        condition = umath.equal(x._data, value)
        mask = x._mask
    else:
        condition = umath.equal(fromnumeric.asarray(x), value)
        mask = nomask
    mask = mask_or(mask, make_mask(condition, small_mask=True))
    return masked_array(x, mask=mask, copy=copy, fill_value=value)

def masked_values(x, value, rtol=1.e-5, atol=1.e-8, copy=True):
    """Masks the array `x` where the data are approximately equal to `value`
(that is, ``abs(x - value) <= atol+rtol*abs(value)``).
Suitable only for floating points. For integers, please use `masked_equal`.
The mask is set to `nomask` if posible.

:Parameters:
    - `rtol` (Float, *[1e-5]*): Tolerance parameter.
    - `atol` (Float, *[1e-8]*): Tolerance parameter.
    - `copy` (boolean, *[False]*) : Returns a copy of `x` if True.
    """
    abs = umath.absolute
    xnew = filled(x, value)
    if issubclass(xnew.dtype.type, numeric.floating):
        condition = umath.less_equal(abs(xnew-value), atol+rtol*abs(value))
        try:
            mask = x._mask
        except AttributeError:
            mask = nomask
    else:
        condition = umath.equal(xnew, value)
        mask = nomask
    mask = mask_or(mask, make_mask(condition, small_mask=True))
    return masked_array(xnew, mask=mask, copy=copy, fill_value=value)

#####--------------------------------------------------------------------------
#---- --- Printing options ---
#####--------------------------------------------------------------------------
class _MaskedPrintOption:
    """Handles the string used to represent missing data in a masked array."""
    def __init__ (self, display):
        "Creates the masked_print_option object."
        self._display = display
        self._enabled = True

    def display(self):
        "Displays the string to print for masked values."
        return self._display

    def set_display (self, s):
        "Sets the string to print for masked values."
        self._display = s

    def enabled(self):
        "Is the use of the display value enabled?"
        return self._enabled

    def enable(self, small_mask=1):
        "Set the enabling small_mask to `small_mask`."
        self._enabled = small_mask

    def __str__ (self):
        return str(self._display)

    __repr__ = __str__

#if you single index into a masked location you get this object.
masked_print_option = _MaskedPrintOption('--')

#####--------------------------------------------------------------------------
#---- --- MaskedArray class ---
#####--------------------------------------------------------------------------
class MaskedArray(numeric.ndarray):
    """Arrays with possibly masked values.
Masked values of True exclude the corresponding element from any computation.

Construction:
    x = array(data, dtype=None, copy=True, order=False,
              mask = nomask, fill_value=None, small_mask=True)

If copy=False, every effort is made not to copy the data:
If `data` is a MaskedArray, and argument mask=nomask, then the candidate data 
is `data._data` and the mask used is `data._mask`. 
If `data` is a numeric array, it is used as the candidate raw data.
If `dtype` is not None and is different from data.dtype.char then a data copy is required.
Otherwise, the candidate is used.

If a data copy is required, the raw (unmasked) data stored is the result of:
numeric.array(data, dtype=dtype.char, copy=copy)

If `mask` is `nomask` there are no masked values. 
Otherwise mask must be convertible to an array of booleans with the same shape as x.
If `small_mask` is True, a mask consisting of zeros (False) only is compressed to `nomask`.
Otherwise, the mask is not compressed.

fill_value is used to fill in masked values when necessary, such as when 
printing and in method/function filled().
The fill_value is not used for computation within this module.
    """
    __array_priority__ = 10.1
    _defaultmask = nomask
    _defaulthardmask = False
    #TODO: There some reorganization to do round here
    def __new__(cls, data=None, mask=nomask, dtype=None, copy=False, fill_value=None,
                keep_mask=True, small_mask=True, hard_mask=False, flag=None,
                **options):
        """array(data, dtype=None, copy=True, mask=nomask, fill_value=None)
        
If `data` is already a ndarray, its dtype becomes the default value of dtype.
        """  
        if flag is not None:
            warnings.warn("The flag 'flag' is now called 'small_mask'!",
                          DeprecationWarning)
            small_mask = flag
        # 1. Argument is MA ...........
        if isinstance(data, MaskedArray) or\
           (hasattr(data,"_mask") and hasattr(data,"_data")) :
            if data is masked:
                return masked.view(cls)
            if keep_mask:
                if mask is nomask:
                    if copy:
                        cls._defaultmask = data._mask.copy()
                    else:
                        cls._defaultmask = data._mask
                else:
                    cls._defaultmask = mask_or(data._mask, mask, 
                                                copy=copy, small_mask=small_mask)
            else:
                cls._defaultmask = make_mask(mask, copy=copy, small_mask=small_mask)
            # Update fille_value
            if fill_value is None:
                cls._fill_value = data._fill_value
            else:
                cls._fill_value = fill_value
            cls._defaulthardmask = hard_mask
            _data = data._data
            if dtype is not None and _data.dtype != numeric.dtype(dtype):
                return _data.astype(dtype).view(cls)
            elif copy:
                return _data.copy().view(cls)
            else:
                return _data.view(cls)
        # 2. Argument is not MA .......
        if isinstance(data, ndarray):
            if dtype is not None and data.dtype != numeric.dtype(dtype):
                _data = data.astype(dtype)
            elif copy:
                _data = data.copy()
            else:
                _data = data
        else:
            _data = numeric.array(data, dtype=dtype, copy=copy)
#            try:
#                _data = numeric.array(data, dtype=dtype, copy=copy)
#            except TypeError:
#                _data = empty(len(data), dtype=dtype)
#                for (k,v) in enumerate(data):
#                    _data[k] = v
#                if mask is nomask:
#                    cls.__defaultmask = getmask(_data)
#                    return _data.view(cls)
        # Define mask .................
        mask = make_mask(mask, copy=copy, small_mask=small_mask)
        #....Check shapes compatibility
        if mask is not nomask:
            (nd, nm) = (_data.size, mask.size)
            if (nm != nd):
                # We need to resize w/ a function, in case _data is only a reference
                if nm == 1:
                    mask = fromnumeric.resize(mask, _data.shape)
                elif nd == 1:
                    _data = fromnumeric.resize(_data, mask.shape)
                else:
                    msg = "Mask and data not compatible: data size is %i, "+\
                          "mask size is %i."
                    raise MAError, msg % (nd, nm)
            elif (mask.shape != _data.shape):
                mask = mask.reshape(_data.shape)
#                mask = _data.shape
        #....
        cls._fill_value = fill_value
        cls._defaulthardmask = hard_mask
        cls._defaultmask = mask
        cls._defaultoptions = options
        return numeric.asanyarray(_data).view(cls)
    #..................................
    def __array_wrap__(self, obj, context=None):
        """Special hook for ufuncs.
Wraps the numpy array and sets the mask according to context.
        """
#        mclass = self.__class__
        #..........
        if context is None:
#            return mclass(obj, mask=self._mask, copy=False)
            return MaskedArray(obj, mask=self._mask, copy=False,
                               dtype=obj.dtype,
                               fill_value=self.fill_value, )
        #..........
        (func, args) = context[:2]
        m = reduce(mask_or, [getmask(arg) for arg in args])
        # Get domain mask
        domain = ufunc_domain.get(func, None)
        if domain is not None:
            m = mask_or(m, domain(*[getattr(arg, '_data', arg) for arg in args]))
        # Update mask
        if m is not nomask:
            try:
                dshape = obj.shape
            except AttributeError:
                pass
            else:
                if m.shape != dshape:
                    m = reduce(mask_or, [getmaskarray(arg) for arg in args])
#        return mclass(obj, copy=False, mask=m)
        return MaskedArray(obj, copy=False, mask=m,)
#                           dtype=obj.dtype, fill_value=self._fill_value)
    #........................
    #TODO: there should be some reorganization to do round here.
    def __array_finalize__(self,obj):
        """Finalizes the masked array.
        """
        #
        if isinstance(obj, MaskedArray):
            # We came here from a MaskedArray
            self._data = obj._data
            self._mask = obj._mask
            self._hardmask = obj._hardmask
            self._fill_value = obj._fill_value
            self.options = obj.options
        else:
            # We came here from a .view()
            if hasattr(obj,'_data') and hasattr(obj, '_mask'):
                # obj is an old masked array or a smart record
                self._data = obj._data
                self._mask = obj._mask
            else:
                # obj is anything but...
                self._data = obj
                self._mask = self._defaultmask
            # Set the instance default
            self._hardmask = self._defaulthardmask
            self.fill_value = self._fill_value
            self.options = self._defaultoptions
            # Reset the class default
            MaskedArray._defaultmask = nomask
            MaskedArray._defaulthardmask = False
            MaskedArray._fill_value = None
#        #
        return    
    #............................................
    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]
Returns the item described by i. Not a copy as in previous versions.
        """
        if getmask(indx) is not nomask:
            msg = "Masked arrays must be filled before they can be used as indices!"
            raise IndexError, msg
        dout = self._data[indx]
        m = self._mask
        scalardout = (len(numeric.shape(dout))==0)
        # 
        if m is nomask:
            if scalardout:
                return dout
            else:
                return self.__class__(dout, mask=nomask, keep_mask=True,
                                      fill_value=self._fill_value,
                                      **self.options)
        #....
        mi = m[indx]
        if mi.size == 1:
            if mi:
                return masked
            return dout
        else:
            return self.__class__(dout, mask=mi, fill_value=self._fill_value,
                                  **self.options)
    #........................
    def __setitem__(self, index, value):
        """x.__setitem__(i, y) <==> x[i]=y
Sets item described by index. If value is masked, masks those locations.
        """
        if self is masked:
            raise MAError, 'Cannot alter the masked element.'
        if getmask(index) is not nomask:
            msg = "Masked arrays must be filled before they can be used as indices!"
            raise IndexError, msg
        #....
        (d, m) = (self._data, self._mask)
        #....
        if value is masked:
            if m is nomask:
                m = make_mask_none(d.shape)
            else:
                m = m.copy()
            m[index] = True
            self._mask = m
            return
        #....   
        if m is nomask:
            d[index] = filled(value)
            valmask = getmask(value)
            if valmask is not nomask:
                m = make_mask_none(d.shape)
                m[index] = valmask
        elif not self._hardmask:
            d[index] = filled(value)
            valmask = getmask(value)
            m = m.copy()
            if valmask is nomask:
                m[index] = False
            else:
                m[index] = valmask
        elif hasattr(index, 'dtype') and (index.dtype==bool_):
            index *= ~m
            d[index] = filled(value)
#        elif isinstance(index, int):
        else:            
            mindx = m[index]
            value = masked_array(value, mask=mindx, keep_mask=True)
            valdata = filled(value)
            valmask = getmask(value)
            if valmask is nomask:
                d[index] = valdata
            elif valmask.size > 1:
                dindx = d[index]
                numeric.putmask(dindx, ~valmask, valdata)
                d[index] = dindx
                numeric.putmask(mindx, valmask, True)
                m[index] = mindx
        #.....    
        if not m.any():
            self._mask = nomask
        else:
            self._mask = m
    #............................................
    def __getslice__(self, i, j):
        """x.__getslice__(i, j) <==> x[i:j]
Returns the slice described by i, j.
The use of negative indices is not supported."""
        m = self._mask
        dout = self._data[i:j]
        if m is nomask:
            return self.__class__(dout, fill_value=self._fill_value,
                                  **self.options)
        else:
            return self.__class__(dout, mask=m[i:j], fill_value=self._fill_value,
                                  **self.options)
    #........................
    def __setslice__(self, i, j, value):
        """x.__setslice__(i, j, value) <==> x[i:j]=value
Sets a slice i:j to `value`. 
If `value` is masked, masks those locations."""
        if self is masked:
            #TODO: Well, maybe we could/should
            raise MAError, "Cannot alter the 'masked' object."
        #....
        (d, m) = (self._data, self._mask)
        #....
        if value is masked:
            if m is nomask:
                m = make_mask_none(d.shape)
            m[i:j] = True
            self._mask = m
            return
        #....   
        if m is nomask:
            valmask = getmask(value)
            valdata = filled(value)
            d[i:j] = valdata
            if valmask is not nomask:
                m = make_mask_none(d.shape)
                m[i:j] = valmask
        elif not self._hardmask:
            valmask = getmask(value)
            valdata = filled(value)
            d[i:j] = valdata
            if valmask is nomask:
                m[i:j] = False
            else:
                m[i:j] = valmask
        else:
            mindx = m[i:j]
            value = masked_array(value, mask=mindx, keep_mask=True)
            valmask = value._mask
            if valmask is nomask:
                d[i:j][~mindx] = filled(value)
            elif valmask.size > 1:
                d[i:j][~mindx] = value[~valmask]
                m[i:j][valmask] = True
        #.....    
        if not m.any():
            self._mask = nomask
        else:
            self._mask = m
    #............................................
    # If we don't want to crash the performance, we better leave __getattribute__ alone...
#    def __getattribute__(self, name):
#        """x.__getattribute__('name') = x.name
#Returns the chosen attribute. 
#If the attribute cannot be directly accessed, checks the _data section.
#        """
#        try:
#            return ndarray.__getattribute__(self, name)
#        except AttributeError:
#            pass
#        try:
#            return self._data.__getattribute__(name)
#        except AttributeError:
#            raise AttributeError
    #............................................
    def __str__(self):
        """x.__str__() <==> str(x)
Calculates the string representation, using masked for fill if it is enabled. 
Otherwise, fills with fill value.
        """
        if masked_print_option.enabled():
            f = masked_print_option
            if self is masked:
                return str(f)
            m = self._mask
            if m is nomask:
                res = self._data
            else:
                if m.shape == () and m:
                    return str(f)
                # convert to object array to make filled work
#CHECK: the two lines below seem more robust than the self._data.astype
#                res = numeric.empty(self._data.shape, object_)
#                numeric.putmask(res,~m,self._data)
                res = self._data.astype("|O8")
                res[self._mask] = f
        else:
            res = self.filled(self.fill_value)
        return str(res)

    def __repr__(self):
        """x.__repr__() <==> repr(x)
Calculates the repr representation, using masked for fill if it is enabled. 
Otherwise fill with fill value.
        """
        with_mask = """\
masked_%(name)s(data =
 %(data)s,
      mask =
 %(mask)s,
      fill_value=%(fill)s)
"""
        with_mask1 = """\
masked_%(name)s(data = %(data)s,
      mask = %(mask)s,
      fill_value=%(fill)s)
"""
        n = len(self.shape)
        name = repr(self._data).split('(')[0]
        if n <= 1:
            return with_mask1 % {
                'name': name,
                'data': str(self),
                'mask': str(self._mask),
                'fill': str(self.fill_value),
                }
        return with_mask % {
            'name': name,
            'data': str(self),
            'mask': str(self._mask),
            'fill': str(self.fill_value),
            }
    #............................................
    def __abs__(self):
        """x.__abs__() <==> abs(x)
Returns a masked array of the current subclass, with the new `_data` 
the absolute of the inital `_data`.
        """
        return self.__class__(self._data.__abs__(), mask=self._mask,
                              fill_value = self._fill_value, **self.options)
    #
    def __neg__(self):
        """x.__abs__() <==> neg(x)
Returns a masked array of the current subclass, with the new `_data` 
the negative of the inital `_data`."""
        try:
            return self.__class__(self._data.__neg__(), mask=self._mask,
                                  fill_value = self._fill_value, **self.options)
        except MAError:
            return negative(self)
    #
    def __iadd__(self, other):
        "Adds other to self in place."
        f = convert_typecode(filled(other, 0), self._data.dtype.char) 
        m = getmask(other)
        self._data += f
        if self._mask is nomask:
            self._mask = m
        elif m is not nomask:
            self._mask += m
        return self
    #
    def __isub__(self, other):
        "Subtracts other from self in place."
        f = convert_typecode(filled(other, 0), self._data.dtype.char) 
        m = getmask(other)
        self._data -= f
        if self._mask is nomask:
            self._mask = m
        elif m is not nomask:
            self._mask += m
        return self
    #
    def __imul__(self, other):
        "Multiplies self by other in place."
        f = convert_typecode(filled(other, 1), self._data.dtype.char) 
        m = getmask(other)
        self._data *= f
        if self._mask is nomask:
            self._mask = m
        elif m is not nomask:
            self._mask += m
        return self
    #
    def __idiv__(self, other):
        "Divides self by other in place."
        f = convert_typecode(filled(other, 0), self._data.dtype.char) 
        mo = getmask(other)
        result = divide(self, masked_array(f, mask=mo))
        self._data = result._data
        dm = result._mask
        if dm is not self._mask:
            self._mask = dm
        return self    

#    #
#    def __eq__(self, other):
#        return equal(self,other)
#
#    def __ne__(self, other):
#        return not_equal(self,other)
#
#    def __lt__(self, other):
#        return less(self,other)
#
#    def __le__(self, other):
#        return less_equal(self,other)
#
#    def __gt__(self, other):
#        return greater(self,other)
#
#    def __ge__(self, other):
#        return greater_equal(self,other)    

    #............................................
    def __float__(self):
        "Converts self to float."
        if self._mask is not nomask:
            warnings.warn("Warning: converting a masked element to nan.")
            return numpy.nan
            #raise MAError, 'Cannot convert masked element to a Python float.'
        return float(self._data.item())

    def __int__(self):
        "Converts self to int."
        if self._mask is not nomask:
            raise MAError, 'Cannot convert masked element to a Python int.'
        return int(self._data.item())    
    
    @property
    def dtype(self):
        """returns the data type of `_data`."""
        return self._data.dtype
    
    def astype (self, tc):
        """Returns self as an array of given type.
Subclassing is preserved."""
        if tc == self._data.dtype:
            return self
        try:
            return self.__class__(self, mask=self._mask, dtype=tc, copy=True,
                                  **self.options)
        except:
#            d = self._data.astype(tc) 
            return self.__class__(self._data.astype(tc), mask=self._mask, 
                                  dtype=tc, **self.options)
#        
#        
    #............................................
    def harden_mask(self):
        "Forces the mask to hard"
        self._hardmask = True
    def soften_mask(self):
        "Forces the mask to soft"
        self._hardmask = False
    #............................................        
    #TODO: FIX THAT: THAT"S NOT A REAL FLATITER
    def _get_flat(self):
        """Calculates the flat value.
        """
        if self._mask is nomask:
            return masked_array(self._data.ravel(), mask=nomask, copy=False,
                                fill_value = self._fill_value,
                                **self.options)
        else:
            return masked_array(self._data.ravel(), mask=self._mask.ravel(),
                                copy=False, fill_value = self._fill_value,
                                **self.options)
    #
    def _set_flat (self, value):
        "x.flat = value"
        y = self.ravel()
        y[:] = value
    #
    flat = property(fget=_get_flat, fset=_set_flat, doc="Flat version")
    #   
    #............................................
    def _get_real(self):
        "Returns the real part of a complex array."
        return self.__class__(self._data.real, mask=self.mask,
                              fill_value = self._fill_value, **self.options)
    def _set_real (self, value):
        "Sets the real part of a complex array to `value`."
        y = self.real
        y[...] = value
        
    real = property(fget=_get_real, fset=_set_real, doc="Get it real!")

    def _get_imaginary(self):
        "Returns the imaginary part of a complex array."
        return self.__class__(self._data.imag, mask=self.mask,
                              fill_value = self._fill_value, **self.options)

    def _set_imaginary (self, value):
        "Sets the imaginary part of a complex array to `value`."
        y = self.imaginary
        y[...] = value
        
    imag = property(fget=_get_imaginary, fset=_set_imaginary, 
                    doc="Imaginary part.")
    imaginary = imag
    #............................................
    def _get_mask(self):
        """Returns the current mask."""
        return self._mask
    def _set_mask(self, mask):
        """Sets the mask to `mask`."""
        mask = make_mask(mask, copy=False, small_mask=True)
        if mask is not nomask:
            if mask.size != self._data.size:
                raise ValueError, "Inconsistent shape between data and mask!"
            if mask.shape != self._data.shape:
                mask.shape = self._data.shape
            self._mask = mask  
        else:
            self._mask = nomask
    mask = property(fget=_get_mask, fset=_set_mask, doc="Mask")
    #............................................
    def get_fill_value(self):
        "Returns the filling value."
        return self._fill_value
    
    def set_fill_value(self, value=None):
        """Sets the filling value to `value`. 
If None, uses the default, based on the data type."""
        if value is None:
            value = default_fill_value(self._data)
        self._fill_value = value
        
    fill_value = property(fget=get_fill_value, fset=set_fill_value,
                          doc="Filling value")
    
    def filled(self, fill_value=None):
        """Returns an array of the same class as `_data`,
 with masked values filled with `fill_value`.
Subclassing is preserved.
        
If `fill_value` is None, uses self.fill_value.
        """
        d = self._data
        m = self._mask
        if m is nomask:
            return d
        #
        if fill_value is None:
            value = self.fill_value
        else:
            value = fill_value
        #
        if self is masked_singleton:
            result = numeric.asanyarray(value)
        else:
            result = d.copy()
            try:
                result[m] = value
            except (TypeError, AttributeError):
                value = numeric.array(value, dtype=object)
                d = d.astype(object)
                result = fromnumeric.choose(m, (d, value))
            except IndexError:
                #ok, if scalar
                if d.shape:
                    raise
                elif m:
                    result = numeric.array(value, dtype=d.dtype)
                else:
                    result = d
        return result
    
    def compressed(self):
        "A 1-D array of all the non-masked data."
        d = self._data.ravel()
        if self._mask is nomask:
            return d
        else:
            return d[~self._mask.ravel()]
    #............................................
    def count(self, axis=None):
        """Counts the non-masked elements of the array along a given axis,
and returns a masked array where the mask is True where all data are masked.
If `axis` is None, counts all the non-masked elements, and returns either a
scalar or the masked singleton."""
        m = self._mask
        s = self._data.shape
        ls = len(s)
        if m is nomask:
            if ls == 0:
                return 1
            if ls == 1:
                return s[0]
            if axis is None:
                return self._data.size
            else:
                n = s[axis]
                t = list(s)
                del t[axis]
                return numeric.ones(t) * n      
        n1 = fromnumeric.size(m, axis)
        n2 = m.astype(int_).sum(axis)
        if axis is None:
            return (n1-n2)
        else:
            return masked_array(n1 - n2)
    #............................................
    def _get_shape(self):
        "Returns the current shape."
        return self._data.shape
    #
    def _set_shape (self, newshape):
        "Sets the array's shape."
        self._data.shape = newshape
        if self._mask is not nomask:
            #self._mask = self._mask.copy()
            self._mask.shape = newshape
    #
    shape = property(fget=_get_shape, fset=_set_shape, 
                     doc="Shape of the array, as a tuple.")
    #
    def _get_size(self):
        "Returns the current size."
        return self._data.size
    size = property(fget=_get_size,
                    doc="Size (number of elements) of the array.")
    #
    def _get_ndim(self):
        "Returns the number of dimensions."
        return self._data.ndim
    ndim = property(fget=_get_ndim,
                    doc="Number of dimensions of the array.")
    #
    def reshape (self, *s):
        """Reshapes the array to shape s.
Returns a new masked array. 
If you want to modify the shape in place, please use `a.shape = s`"""
        if self._mask is not nomask:
            return self.__class__(self._data.reshape(*s), 
                                  mask=self._mask.reshape(*s),
                                  fill_value=self.fill_value, **self.options)
        else:
            return self.__class__(self._data.reshape(*s),
                                  fill_value=self.fill_value, **self.options)
    #
    def repeat(self, repeats, axis=None):
        """Repeat elements of `a` `repeats` times along `axis`.
`repeats` is a sequence of length `a.shape[axis]` telling how many times 
each element should be repeated.
The mask is repeated accordingly.
        """
        f = self.filled()
        if isinstance(repeats, types.IntType):
            if axis is None:
                num = f.size
            else:
                num = f.shape[axis]
            repeats = tuple([repeats]*num)
    
        m = self._mask
        if m is not nomask:
            m = fromnumeric.repeat(m, repeats, axis)
        d = fromnumeric.repeat(f, repeats, axis)
        return self.__class__(d, mask=m, fill_value=self.fill_value, 
                              **self.options)
    #
    def resize(self, newshape, refcheck=True, order=False):
        """Attempts to modify size and shape of self inplace.  
        The array must own its own memory and not be referenced by other arrays.    
        Returns None.
        """
        try:
            self._data.resize(newshape,)
            if self.mask is not nomask:
                self._mask.resize(newshape,)
        except ValueError:
            msg = "Cannot resize an array that has been referenced or "+\
                  "is referencing another array in this way.\n"+\
                  "Use the resize function."
            raise ValueError, msg
        return None
    #
    def flatten(self):
        """Flattens the array in place.
        """
        flatsize = self.size
        self._data.resize((flatsize,))
        if self.mask is not nomask:
            self._mask.resize((flatsize,))
        return self

    #
    def put(self, indices, values, mode='raise'):
        """Sets storage-indexed locations to corresponding values.
a.put(values, indices, mode) sets a.flat[n] = values[n] for each n in indices.
`values` can be scalar or an array shorter than indices, and it will be repeat,
if necessary.
If `values` has some masked values, the initial mask is updated in consequence,
else the corresponding values are unmasked.
        """
        #TODO: Check that
        (d, m) = (self._data, self._mask)
        ind = filled(indices)
        v = filled(values)
        d.put(ind, v, mode=mode)
        if m is not nomask:
            if getmask(values) is not nomask:
                m.put(ind, values._mask, mode=mode)
            else:
                m.put(ind, False, mode=mode)
            self._mask = make_mask(m, copy=False, small_mask=True)
    #............................................
    def ids (self):
        """Return the ids of the data and mask areas."""
        return (id(self._data), id(self._mask))
    #............................................
    def all(self, axis=None):
        """a.all(axis) returns True if all entries along the axis are True.
    Returns False otherwise. If axis is None, uses the flatten array.
    Masked data are considered as True during computation.
    Outputs a masked array, where the mask is True if all data are masked along the axis.
        """
        d = filled(self, True).all(axis)
        m = self._mask.all(axis)
        return self.__class__(d, mask=m, dtype=bool_,
                              fill_value=self._fill_value, **self.options)
    def any(self, axis=None):
        """a.any(axis) returns True if some or all entries along the axis are True.
    Returns False otherwise. If axis is None, uses the flatten array.
    Masked data are considered as False during computation.
    Outputs a masked array, where the mask is True if all data are masked along the axis.
        """
        d = filled(self, False).any(axis)
        m = self._mask.all(axis)
        return self.__class__(d, mask=m, dtype=bool_,
                              fill_value=self._fill_value, **self.options) 
    def nonzero(self):
        """a.nonzero() returns a tuple of arrays

    Returns a tuple of arrays, one for each dimension of a,
    containing the indices of the non-zero elements in that
    dimension.  The corresponding non-zero values can be obtained
    with
        a[a.nonzero()].

    To group the indices by element, rather than dimension, use
        transpose(a.nonzero())
    instead. The result of this is always a 2d array, with a row for
    each non-zero element."""
        return numeric.asarray(self.filled(0)).nonzero()
    #............................................
    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        """a.trace(offset=0, axis1=0, axis2=1, dtype=None, out=None)
Returns the sum along the offset diagonal of the array's indicated `axis1` and `axis2`.
        """
        #TODO: What are we doing with `out`?
        (d,m) = (self._data, self._mask)
        if m is nomask:
            return d.trace(offset=offset, axis1=axis1, axis2=axis2, 
                           out=out).astype(dtype)
        else:
            D = self.diagonal(offset=offset, axis1=axis1, axis2=axis2, 
                              ).astype(dtype)
            return D.sum(axis=None)
    #............................................
    def sum(self, axis=None, dtype=None):
        """a.sum(axis=None, dtype=None) 
Sums the array `a` over the given axis `axis`.
Masked values are set to 0.
If `axis` is None, applies to a flattened version of the array.
    """
        if self._mask is nomask:
#            if axis is None:
#                return self._data.sum(None, dtype=dtype)
            return self.__class__(self._data.sum(axis, dtype=dtype),
                                  mask=nomask, fill_value=self.fill_value,
                                  **self.options)
        else:
#            if axis is None:
#                return self.filled(0).sum(None, dtype=dtype)
            return self.__class__(self.filled(0).sum(axis, dtype=dtype),
                                  mask=self._mask.all(axis),
                                  fill_value=self.fill_value, **self.options)
            
    def cumsum(self, axis=None, dtype=None):
        """a.cumprod(axis=None, dtype=None)
Returns the cumulative sum of the elements of array `a` along the given axis `axis`. 
Masked values are set to 0.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
#            if axis is None:
#                return self._data.cumsum(None, dtype=dtype)
            return self.__class__(self._data.cumsum(axis=axis, dtype=dtype),
                                  fill_value=self.fill_value, **self.options)
        else:
#            if axis is None:
#                return self.filled(0).cumsum(None, dtype=dtype)
            return self.__class__(self.filled(0).cumsum(axis=axis, dtype=dtype),
                                  mask=self._mask, fill_value=self.fill_value,
                                  **self.options)
        
    def prod(self, axis=None, dtype=None):
        """a.prod(axis=None, dtype=None)
Returns the product of the elements of array `a` along the given axis `axis`. 
Masked elements are set to 1.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
#            if axis is None:
#                return self._data.prod(None, dtype=dtype)
            return self.__class__(self._data.prod(axis, dtype=dtype),
                                  mask=nomask, fill_value=self.fill_value,
                                  **self.options)
#            return self.__class__(self._data.prod(axis=axis, dtype=dtype))
        else:
#            if axis is None:
#                return self.filled(1).prod(None, dtype=dtype)
            return self.__class__(self.filled(1).prod(axis=axis, dtype=dtype),
                                  mask=self._mask.all(axis),
                                  fill_value=self.fill_value,
                                  **self.options)
    product = prod
            
    def cumprod(self, axis=None, dtype=None):
        """a.cumprod(axis=None, dtype=None)
Returns the cumulative product of ethe lements of array `a` along the given axis `axis`. 
Masked values are set to 1.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
#            if axis is None:
#                return self._data.cumprod(None, dtype=dtype)
            return self.__class__(self._data.cumprod(axis=axis, dtype=dtype),
                                  mask=nomask, fill_value=self.fill_value,
                                  **self.options)
        else:
#            if axis is None:
#                return self.filled(1).cumprod(None, dtype=dtype)
            return self.__class__(self.filled(1).cumprod(axis=axis, dtype=dtype),
                                  mask=self._mask,
                                  fill_value=self.fill_value, **self.options)        
            
    def mean(self, axis=None, dtype=None):
        """a.mean(axis=None, dtype=None)
        
    Averages the array over the given axis.  If the axis is None,
    averages over all dimensions of the array.  Equivalent to

      a.sum(axis, dtype) / size(a, axis).

    The optional dtype argument is the data type for intermediate
    calculations in the sum.

    Returns a masked array, of the same class as a.
        """
        if self._mask is nomask:
#            if axis is None:
#                return self._data.mean(axis=None, dtype=dtype)
            return self.__class__(self._data.mean(axis=axis, dtype=dtype),
                                  mask=nomask, fill_value=self.fill_value,
                                  **self.options)
        else:
            dsum = fromnumeric.sum(self.filled(0), axis=axis, dtype=dtype)
            cnt = self.count(axis=axis)
            mask = self._mask.all(axis)
            if axis is None and mask:
                return masked
            return self.__class__(dsum*1./cnt, mask=mask,
                                  fill_value=self.fill_value, **self.options)
    
    def anom(self, axis=None, dtype=None):
        """a.anom(axis=None, dtype=None)
    Returns the anomalies, or deviation from the average.
            """    
        m = self.mean(axis, dtype)
        if not axis:
            return (self - m)
        else:
            return (self - expand_dims(m,axis))
 
    def var(self, axis=None, dtype=None):
        """a.var(axis=None, dtype=None)
Returns the variance, a measure of the spread of a distribution.

The variance is the average of the squared deviations from the mean,
i.e. var = mean((x - x.mean())**2).
        """
        if self._mask is nomask:
#            if axis is None:
#                return self._data.var(axis=None, dtype=dtype)
            return self.__class__(self._data.var(axis=axis, dtype=dtype),
                                  mask=nomask,
                                  fill_value=self.fill_value, **self.options)
        else:
            cnt = self.count(axis=axis)
            danom = self.anom(axis=axis, dtype=dtype)
            danom *= danom
            dvar = danom.sum(axis) / cnt
#            dvar /= cnt
            if axis is None:
                return dvar
            return self.__class__(dvar, 
                                  mask=mask_or(self._mask.all(axis), (cnt==1)),
                                  fill_value=self.fill_value, **self.options)
            
    def std(self, axis=None, dtype=None):
        """a.std(axis=None, dtype=None)
Returns the standard deviation, a measure of the spread of a distribution.

The standard deviation is the square root of the average of the squared
deviations from the mean, i.e. std = sqrt(mean((x - x.mean())**2)).
        """
        dvar = self.var(axis,dtype)
        if axis is None:
            if dvar is masked:
                return masked
            else:
                # Should we use umath.sqrt instead ?
                return sqrt(dvar)
        return self.__class__(sqrt(dvar._data), mask=dvar._mask,
                              dtype = self.dtype,
                              fill_value=self.fill_value, **self.options)
    #............................................
    def argsort(self, axis=None, fill_value=None, kind='quicksort'):
        """Returns an array of indices that sort 'a' along the specified axis.
    Masked values are filled beforehand to `fill_value`.        
    If `fill_value` is None, uses the default for the data type.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `kind` : String *['quicksort']*
        Sorting algorithm (default 'quicksort')
        Possible values: 'quicksort', 'mergesort', or 'heapsort'

    Returns: array of indices that sort 'a' along the specified axis.

    This method executes an indirect sort along the given axis using the
    algorithm specified by the kind keyword. It returns an array of indices of
    the same shape as 'a' that index data along the given axis in sorted order.

    The various sorts are characterized by average speed, worst case
    performance, need for work space, and whether they are stable. A stable
    sort keeps items with the same key in the same relative order. The three
    available algorithms have the following properties:

    |------------------------------------------------------|
    |    kind   | speed |  worst case | work space | stable|
    |------------------------------------------------------|
    |'quicksort'|   1   | O(n^2)      |     0      |   no  |
    |'mergesort'|   2   | O(n*log(n)) |    ~n/2    |   yes |
    |'heapsort' |   3   | O(n*log(n)) |     0      |   no  |
    |------------------------------------------------------|

    All the sort algorithms make temporary copies of the data when the sort is not
    along the last axis. Consequently, sorts along the last axis are faster and use
    less space than sorts along other axis.
        """
        if fill_value is None:
            fill_value = default_fill_value(self._data)
        d = self.filled(fill_value)
        if axis is None:
            return d.argsort(kind=kind)
        return d.argsort(axis, kind)

    def argmin(self, axis=None, fill_value=None):
        """Returns the array of indices for the minimum values of `a` along the 
    specified axis.
    Masked values are treated as if they had the value `fill_value`.        
    If `fill_value` is None, the default for the data type is used.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `fill_value` : var *[None]*
        Default filling value. If None, uses the data type default.
        """
        if fill_value is None:
            fill_value = default_fill_value(self._data)
        d = self.filled(fill_value)
        if axis is None:
            return d.argmin()
        return d.argmin(axis)

    def argmax(self, axis=None, fill_value=None):
        """Returns the array of indices for the maximum values of `a` along the 
    specified axis.
    Masked values are treated as if they had the value `fill_value`.        
    If `fill_value` is None, the default for the data type is used.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `fill_value` : var *[None]*
        Default filling value. If None, uses the data type default.
        """
        if fill_value is None:
            fill_value = default_fill_value(self._data)
            try:
                fill_value = - fill_value
            except:
                pass
        d = self.filled(fill_value)
        if axis is None:
            return d.argmax()
        return d.argmax(axis)
    
    def sort(self, axis=-1, kind='quicksort', order=None, endwith=True):
        """
        Sort a along the given axis.

    Keyword arguments:

    axis  -- axis to be sorted (default -1)
    kind  -- sorting algorithm (default 'quicksort')
             Possible values: 'quicksort', 'mergesort', or 'heapsort'.
    order -- If a has fields defined, then the order keyword can be the
             field name to sort on or a list (or tuple) of field names
             to indicate the order that fields should be used to define
             the sort.
    endwith--Boolean flag indicating whether missing values (if any) should
             be forced in the upper indices (at the end of the array) or
             lower indices (at the beginning).

    Returns: None.

    This method sorts 'a' in place along the given axis using the algorithm
    specified by the kind keyword.

    The various sorts may characterized by average speed, worst case
    performance, need for work space, and whether they are stable. A stable
    sort keeps items with the same key in the same relative order and is most
    useful when used with argsort where the key might differ from the items
    being sorted. The three available algorithms have the following properties:

    |------------------------------------------------------|
    |    kind   | speed |  worst case | work space | stable|
    |------------------------------------------------------|
    |'quicksort'|   1   | O(n^2)      |     0      |   no  |
    |'mergesort'|   2   | O(n*log(n)) |    ~n/2    |   yes |
    |'heapsort' |   3   | O(n*log(n)) |     0      |   no  |
    |------------------------------------------------------|

    All the sort algorithms make temporary copies of the data when the sort is
    not along the last axis. Consequently, sorts along the last axis are faster
    and use less space than sorts along other axis.

    """
        
        if endwith:
            filler = minimum_fill_value(self.dtype)
        else:
            filler = maximum_fill_value(self.dtype)
        indx = self.filled(filler).argsort(axis=axis,kind=kind,order=order)
        self._data = self._data[indx]
        m = self._mask
        if m is not nomask:
            self._mask = m[indx]
        return 
    #............................................
    # Backwards Compatibility. Heck...
    @property
    def data(self):
        """Returns the `_data` part of the MaskedArray. 
You should really use `_data` instead..."""
        return self._data
    def raw_data(self):
        """Returns the `_data` part of the MaskedArray. 
You should really use `_data` instead..."""
        return self._data

##..............................................................................



#class _arithmethods:
#    """Defines a wrapper for arithmetic methods.
#Instead of directly calling a ufunc, the corresponding method of  the `array._data` 
#object is called instead.
#    """
#    def __init__ (self, methodname, fill_self=0, fill_other=0, domain=None):
#        """
#:Parameters:
#    - `methodname` (String) : Method name.
#    - `fill_self` (Float *[0]*) : Fill value for the instance.
#    - `fill_other` (Float *[0]*) : Fill value for the target.
#    - `domain` (Domain object *[None]*) : Domain of non-validity. 
#        """
#        self.methodname = methodname
#        self.fill_self = fill_self
#        self.fill_other = fill_other
#        self.domain = domain
#    #
#    def __call__ (self, instance, other, *args):
#        "Execute the call behavior."
#        m_self = instance._mask
#        m_other = getmask(other)
#        base = filled(instance,self.fill_self)
#        target = filled(other, self.fill_other)
#        if self.domain is not None:
#            # We need to force the domain to a ndarray only.
#            if self.fill_other > self.fill_self:
#                domain = self.domain(base, target)
#            else:
#                domain = self.domain(target, base)
#            if domain.any():
#                #If `other` is a subclass of ndarray, `filled` must have the 
#                # same subclass, else we'll lose some info.
#                #The easiest then is to fill `target` instead of creating
#                # a pure ndarray.
#                #Oh, and we better make a copy!
#                if isinstance(other, ndarray):
#                    if target is other:
#                        # We don't want to modify other: let's copy target, then
#                        target = target.copy()
#                    target[:] = numeric.where(fromnumeric.asarray(domain), 
#                                              self.fill_other, target)
#                else:
#                    target = numeric.where(fromnumeric.asarray(domain), 
#                                           self.fill_other, target)
#                m_other = mask_or(m_other, domain)
#        m = mask_or(m_self, m_other)      
#        method = getattr(base, self.methodname)       
#        return instance.__class__(method(target, *args), mask=m)
#    #
#    def patch(self):
#        """Applies the method `func` from class `method` to MaskedArray"""
#        return types.MethodType(self,None,MaskedArray)
#..............................................................................
class _arithmethods(object):
    """Defines a wrapper for arithmetic methods.
Instead of directly calling a ufunc, the corresponding method of  the `array._data` 
object is called instead.
    """
    def __init__ (self, methodname, fill_self=0, fill_other=0, domain=None):
        """
:Parameters:
    - `methodname` (String) : Method name.
    - `fill_self` (Float *[0]*) : Fill value for the instance.
    - `fill_other` (Float *[0]*) : Fill value for the target.
    - `domain` (Domain object *[None]*) : Domain of non-validity. 
        """
        self.methodname = methodname
        self.fill_self = fill_self
        self.fill_other = fill_other
        self.domain = domain
        self.obj = None
        self.__doc__ = self.getdoc()
    #
    def getdoc(self):
        "Returns the doc of the function (from the doc of the method)."
        try:
            return getattr(MaskedArray, self.methodname).__doc__
        except:
            return getattr(numpy, self.methodname).__doc__
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, other, *args):
        "Execute the call behavior."
        instance = self.obj
        m_self = instance._mask
        m_other = getmask(other)
        base = filled(instance,self.fill_self)
        target = filled(other, self.fill_other)
        if self.domain is not None:
            # We need to force the domain to a ndarray only.
            if self.fill_other > self.fill_self:
                domain = self.domain(base, target)
            else:
                domain = self.domain(target, base)
            if domain.any():
                #If `other` is a subclass of ndarray, `filled` must have the 
                # same subclass, else we'll lose some info.
                #The easiest then is to fill `target` instead of creating
                # a pure ndarray.
                #Oh, and we better make a copy!
                if isinstance(other, ndarray):
                    if target is other or target is base:
                        # We don't want to modify other: let's copy target, then
                        # Same if target us base, instead...
                        target = target.copy()
                    target[:] = numeric.where(fromnumeric.asarray(domain), 
                                              self.fill_other, target)
                else:
                    target = numeric.where(fromnumeric.asarray(domain), 
                                           self.fill_other, target)
                m_other = mask_or(m_other, domain)
        m = mask_or(m_self, m_other)      
        method = getattr(base, self.methodname)       
        return instance.__class__(method(target, *args), mask=m,
                                  fill_value=instance.fill_value,
                                  **instance.options)
#......................................
class _compamethods(object):
    """Defines comparison methods (eq, ge, gt...).
Instead of calling a ufunc, the method of the masked object is called.
    """
    def __init__ (self, methodname, fill_self=0, fill_other=0):
        """
:Parameters:
    - `methodname` (String) : Method name.
    - `fill_self` (Float *[0]*) : Fill value for the instance.
    - `fill_other` (Float *[0]*) : Fill value for the target.
    - `domain` (Domain object *[None]*) : Domain of non-validity. 
        """
        self.methodname = methodname
        self.fill_self = fill_self
        self.fill_other = fill_other
        self.obj = None
        self.__doc__ = self.getdoc()
    #
    def getdoc(self):
        "Returns the doc of the function (from the doc of the method)."
        try:
            return getattr(MaskedArray, self.methodname).__doc__
        except:
            return getattr(numpy, self.methodname).__doc__
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__ (self, other, *args):
        "Execute the call behavior."
        instance = self.obj
        m = mask_or(instance._mask, getmask(other), small_mask=False)
        base = instance.filled(self.fill_self)
        target = filled(other, self.fill_other)
        method = getattr(base, self.methodname)      
        return instance.__class__(method(target, *args), mask=m,
                                  **instance.options)  
#..........................................................
MaskedArray.__add__ = _arithmethods('__add__')
MaskedArray.__radd__ = _arithmethods('__add__')
MaskedArray.__sub__ = _arithmethods('__sub__')
MaskedArray.__rsub__ = _arithmethods('__rsub__')
MaskedArray.__pow__ = _arithmethods('__pow__')
MaskedArray.__mul__ = _arithmethods('__mul__', 1, 1)
MaskedArray.__rmul__ = _arithmethods('__mul__', 1, 1)
MaskedArray.__div__ = _arithmethods('__div__', 0, 1, 
                                    domain_safe_divide())
MaskedArray.__rdiv__ = _arithmethods('__rdiv__', 1, 0, 
                                     domain_safe_divide())
MaskedArray.__truediv__ = _arithmethods('__truediv__', 0, 1, 
                                        domain_safe_divide())
MaskedArray.__rtruediv__ = _arithmethods('__rtruediv__', 1, 0, 
                                         domain_safe_divide())
MaskedArray.__floordiv__ = _arithmethods('__floordiv__', 0, 1, 
                                         domain_safe_divide())
MaskedArray.__rfloordiv__ = _arithmethods('__rfloordiv__', 1, 0, 
                                          domain_safe_divide())
MaskedArray.__eq__ = _compamethods('__eq__')
MaskedArray.__ne__ = _compamethods('__ne__')
MaskedArray.__le__ = _compamethods('__le__')
MaskedArray.__lt__ = _compamethods('__lt__')
MaskedArray.__ge__ = _compamethods('__ge__')
MaskedArray.__gt__ = _compamethods('__gt__')
#####--------------------------------------------------------------------------
#---- --- Shortcuts ---
#####---------------------------------------------------------------------------
def isMaskedArray (x):
    "Is x a masked array, that is, an instance of MaskedArray?"
    return isinstance(x, MaskedArray)
isarray = isMaskedArray
isMA = isMaskedArray  #backward compatibility
#masked = MaskedArray(0, int, mask=1)
masked_singleton = MaskedArray(0, dtype=int_, mask=True)
masked = masked_singleton

masked_array = MaskedArray
def array(data, dtype=None, copy=False, order=False, mask=nomask, 
          keep_mask=True, small_mask=True, hard_mask=None, fill_value=None):
    """array(data, dtype=None, copy=True, order=False, mask=nomask, 
             keep_mask=True, small_mask=True, fill_value=None)
Acts as shortcut to MaskedArray, with options in a different order for convenience.
And backwards compatibility...
    """
    return MaskedArray(data, mask=mask, dtype=dtype, copy=copy, 
                       keep_mask=keep_mask, small_mask=small_mask, 
                       hard_mask=hard_mask, fill_value=fill_value)

def is_masked(x):
    """Returns whether x has some masked values."""
    m = getmask(x)
    if m is nomask:
        return False
    elif m.any():
        return True
    return False


#####--------------------------------------------------------------------------
#---- --- Patch methods ---
#####--------------------------------------------------------------------------
#class _arraymethod:
#    """Defines a wrapper for basic array methods.
#Upon call, returns a masked array, where the new `_data` array is the output
#of the corresponding method called on the original `_data`.
#
#If `onmask` is True, the new mask is the output of the method calld on the initial mask.
#If `onmask` is False, the new mask is just a reference to the initial mask.
#    
#:Parameters:
#    `funcname` : String
#        Name of the function to apply on data. 
#    `onmask` : Boolean *[True]*
#        Whether the mask must be processed also (True) or left alone (False).
#    """
#    def __init__(self, funcname, onmask=True):
#        self._name = funcname
#        self._onmask = onmask
#        self.__doc__ = getattr(ndarray, self._name).__doc__
#    def __call__(self, instance, *args, **params):
#        methodname = self._name
#        (d,m) = (instance._data, instance._mask)
#        C = instance.__class__
#        if m is nomask:
#            return C(getattr(d,methodname).__call__(*args, **params))
#        elif self._onmask:
#            return C(getattr(d,methodname).__call__(*args, **params),
#                     mask=getattr(m,methodname)(*args, **params) )
#        else:
#            return C(getattr(d,methodname).__call__(*args, **params), mask=m) 
#        
#    def patch(self):
#        "Adds the new method to MaskedArray."
#        return types.MethodType(self, None, MaskedArray)
##......................................
#MaskedArray.conj = MaskedArray.conjugate = _arraymethod('conjugate').patch()  
#MaskedArray.diagonal = _arraymethod('diagonal').patch()
#MaskedArray.take = _arraymethod('take').patch()
#MaskedArray.ravel = _arraymethod('ravel').patch()
#MaskedArray.transpose = _arraymethod('transpose').patch()
#MaskedArray.T = _arraymethod('transpose').patch()
#MaskedArray.swapaxes = _arraymethod('swapaxes').patch()
#MaskedArray.clip = _arraymethod('clip', onmask=False).patch()
#MaskedArray.compress = _arraymethod('compress').patch()
#MaskedArray.resize = _arraymethod('resize').patch()
#MaskedArray.copy = _arraymethod('copy').patch()

class _arraymethod(object):
    """Defines a wrapper for basic array methods.
Upon call, returns a masked array, where the new `_data` array is the output
of the corresponding method called on the original `_data`.

If `onmask` is True, the new mask is the output of the method calld on the initial mask.
If `onmask` is False, the new mask is just a reference to the initial mask.
    
:Parameters:
    `funcname` : String
        Name of the function to apply on data. 
    `onmask` : Boolean *[True]*
        Whether the mask must be processed also (True) or left alone (False).
    """
    def __init__(self, funcname, onmask=True):
        self._name = funcname
        self._onmask = onmask
        self.obj = None
        self.__doc__ = self.getdoc()
    #
    def getdoc(self):
        "Returns the doc of the function (from the doc of the method)."
        try:
            return getattr(MaskedArray, self._name).__doc__
        except:
            return getattr(numpy, self._name).__doc__
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__(self, *args, **params):
        methodname = self._name
        obj = self.obj
        (d, m) = (obj._data, obj._mask)
        (t, f) = (obj.dtype, obj._fill_value)
        C = self.obj.__class__
        if m is nomask:
            return C(getattr(d,methodname).__call__(*args, **params),
                     dtype=t, fill_value=f)
        elif self._onmask:
            return C(getattr(d,methodname).__call__(*args, **params),
                     mask=getattr(m,methodname)(*args, **params),
                     dtype=t, fill_value=f, **obj.options)
        else:
            return C(getattr(d,methodname).__call__(*args, **params), mask=m,
                     dtype=t, fill_value=f, **obj.options) 
#......................................
MaskedArray.conj = MaskedArray.conjugate = _arraymethod('conjugate') 
MaskedArray.copy = _arraymethod('copy') 
MaskedArray.diagonal = _arraymethod('diagonal')
MaskedArray.take = _arraymethod('take')
MaskedArray.ravel = _arraymethod('ravel')
MaskedArray.transpose = _arraymethod('transpose')
MaskedArray.T = property(fget=lambda self:self.transpose())
MaskedArray.swapaxes = _arraymethod('swapaxes')
MaskedArray.clip = _arraymethod('clip', onmask=False)
MaskedArray.compress = _arraymethod('compress')
MaskedArray.copy = _arraymethod('copy')
MaskedArray.squeeze = _arraymethod('squeeze')

#####--------------------------------------------------------------------------
#---- --- Extrema functions ---
#####--------------------------------------------------------------------------
class _minimum_operation:
    "Object to calculate minima"
    def __init__ (self):
        """minimum(a, b) or minimum(a)
In one argument case, returns the scalar minimum.
        """
        pass
    #.........
    def __call__ (self, a, b=None):
        "Execute the call behavior."
        if b is None:
            m = getmask(a)
            if m is nomask:
                d = amin(filled(a).ravel())
                return d
            ac = a.compressed()
            if len(ac) == 0:
                return masked
            else:
                return amin(ac)
        else:
            return where(less(a, b), a, b)
    #.........
    def reduce(self, target, axis=0):
        """Reduces `target` along the given `axis`."""
        m = getmask(target)
        if m is nomask:
            t = filled(target)
            return masked_array (umath.minimum.reduce (t, axis))
        else:
            t = umath.minimum.reduce(filled(target, minimum_fill_value(target)), 
                                     axis)
            m = umath.logical_and.reduce(m, axis)
#            return masked_array(t, mask=m, fill_value=get_fill_value(target))
            try:
                return target.__class__(t, mask=m, dtype=t.dtype,
                                        fill_value=get_fill_value(target))
            except AttributeError:
                return masked_array(t, mask=m, dtype=t.dtype,
                                    fill_value=get_fill_value(target))
    #.........
    def outer(self, a, b):
        "Returns the function applied to the outer product of a and b."
        ma = getmask(a)
        mb = getmask(b)
        if ma is nomask and mb is nomask:
            m = nomask
        else:
            ma = getmaskarray(a)
            mb = getmaskarray(b)
            m = logical_or.outer(ma, mb)
        d = umath.minimum.outer(filled(a), filled(b))
        return masked_array(d, mask=m)

def min(array, axis=None, out=None):
    """Returns the minima along the given axis.
If `axis` is None, applies to the flattened array."""
    if out is not None:
        raise TypeError("Output arrays Unsupported for masked arrays")
    if axis is None:
        return minimum(array)
    else:
        return minimum.reduce(array, axis)
#................................................
class _maximum_operation:
    "Object to calculate maxima"
    def __init__ (self):
        """maximum(a, b) or maximum(a)
           In one argument case returns the scalar maximum.
        """
        pass
    #.........
    def __call__ (self, a, b=None):
        "Executes the call behavior."
        if b is None:
            m = getmask(a)
            if m is nomask:
                d = amax(filled(a).ravel())
                return d
            ac = a.compressed()
            if len(ac) == 0:
                return masked
            else:
                return amax(ac)
        else:
            return where(greater(a, b), a, b)
    #.........
    def reduce (self, target, axis=0):
        """Reduces target along the given axis."""
        m = getmask(target)
        if m is nomask:
            t = filled(target)
            return masked_array(umath.maximum.reduce (t, axis))
        else:
            t = umath.maximum.reduce(filled(target, maximum_fill_value(target)), 
                                     axis)
            m = umath.logical_and.reduce(m, axis)
            try:
                return target.__class__(t, mask=m, dtype=t.dtype,
                                        fill_value=get_fill_value(target))
            except AttributeError:
                return masked_array(t, mask=m, dtype=t.dtype,
                                    fill_value=get_fill_value(target))
    #.........
    def outer (self, a, b):
        "Returns the function applied to the outer product of a and b."
        ma = getmask(a)
        mb = getmask(b)
        if ma is nomask and mb is nomask:
            m = nomask
        else:
            ma = getmaskarray(a)
            mb = getmaskarray(b)
            m = logical_or.outer(ma, mb)
        d = umath.maximum.outer(filled(a), filled(b))
        return masked_array(d, mask=m)

def max(obj, axis=None, out=None):
    """Returns the maxima along the given axis.
If `axis` is None, applies to the flattened array."""
    if out is not None:
        raise TypeError("Output arrays Unsupported for masked arrays")
    if axis is None:
        return maximum(obj)
    else:
        return maximum.reduce(obj, axis)
#................................................
def ptp(obj, axis=None):
    """a.ptp(axis=None) =  a.max(axis)-a.min(axis)"""
    try:
        return obj.max(axis)-obj.min(axis)
    except AttributeError:
        return max(obj, axis=axis) - min(obj, axis=axis)
#................................................
MaskedArray.min = min
MaskedArray.max = max
MaskedArray.ptp = ptp

#####---------------------------------------------------------------------------
#---- --- Definition of functions from the corresponding methods ---
#####---------------------------------------------------------------------------
class _frommethod:
    """Defines functions from existing MaskedArray methods.
:ivar _methodname (String): Name of the method to transform.
    """
    def __init__(self, methodname):
        self._methodname = methodname
        self.__doc__ = self.getdoc()
    def getdoc(self):
        "Returns the doc of the function (from the doc of the method)."
        try:
            return getattr(MaskedArray, self._methodname).__doc__
        except:
            return getattr(numpy, self._methodname).__doc__
    def __call__(self, x, *args, **params):
        if isinstance(x, MaskedArray):
            return getattr(x, self._methodname).__call__(*args, **params)
        #FIXME: As x is not a MaskedArray, we transform it to a ndarray with asarray
        #FIXME: ... and call the corresponding method.
        #FIXME: Except that sometimes it doesn't work (try reshape([1,2,3,4],(2,2)))
        #FIXME: we end up with a "SystemError: NULL result without error in PyObject_Call"
        #FIXME: A dirty trick is then to call the initial numpy function...
        method = getattr(fromnumeric.asarray(x), self._methodname)
        try:
            return method(*args, **params)
        except SystemError:
            return getattr(numpy,self._methodname).__call__(x, *args, **params)

all = _frommethod('all')
anomalies = anom = _frommethod('anom')
any = _frommethod('any')
conjugate = _frommethod('conjugate')
ids = _frommethod('ids')
nonzero = _frommethod('nonzero')
diagonal = _frommethod('diagonal')
maximum = _maximum_operation()
mean = _frommethod('mean')
minimum = _minimum_operation ()
product = _frommethod('prod')
ptp = _frommethod('ptp')
ravel = _frommethod('ravel')
repeat = _frommethod('repeat')
std = _frommethod('std')
sum = _frommethod('sum')
swapaxes = _frommethod('swapaxes')
take = _frommethod('take')
var = _frommethod('var')

#..............................................................................
def power(a, b, third=None):
    """Computes a**b elementwise.
    Masked values are set to 1."""
    if third is not None:
        raise MAError, "3-argument power not supported."
    ma = getmask(a)
    mb = getmask(b)
    m = mask_or(ma, mb)
    fa = filled(a, 1)
    fb = filled(b, 1)
    if fb.dtype.char in typecodes["Integer"]:
        return masked_array(umath.power(fa, fb), m)
    md = make_mask((fa < 0), small_mask=1)
    m = mask_or(m, md)
    if m is nomask:
        return masked_array(umath.power(fa, fb))
    else:
        fa[m] = 1
        return masked_array(umath.power(fa, fb), m)

#..............................................................................
def argsort(a, axis=None, kind='quicksort', fill_value=None):
    """Returns an array of indices that sort 'a' along the specified axis.
    Masked values are filled beforehand to `fill_value`.        
    If `fill_value` is None, uses the default for the data type.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `kind` : String *['quicksort']*
        Sorting algorithm (default 'quicksort')
        Possible values: 'quicksort', 'mergesort', or 'heapsort'

    Returns: array of indices that sort 'a' along the specified axis.

    This method executes an indirect sort along the given axis using the
    algorithm specified by the kind keyword. It returns an array of indices of
    the same shape as 'a' that index data along the given axis in sorted order.

    The various sorts are characterized by average speed, worst case
    performance, need for work space, and whether they are stable. A stable
    sort keeps items with the same key in the same relative order. The three
    available algorithms have the following properties:

    |------------------------------------------------------|
    |    kind   | speed |  worst case | work space | stable|
    |------------------------------------------------------|
    |'quicksort'|   1   | O(n^2)      |     0      |   no  |
    |'mergesort'|   2   | O(n*log(n)) |    ~n/2    |   yes |
    |'heapsort' |   3   | O(n*log(n)) |     0      |   no  |
    |------------------------------------------------------|

    All the sort algorithms make temporary copies of the data when the sort is not
    along the last axis. Consequently, sorts along the last axis are faster and use
    less space than sorts along other axis.
    """
    if fill_value is None:
        fill_value = default_fill_value(a)
    d = filled(a, fill_value)
    if axis is None:
        return d.argsort(kind=kind)
    return d.argsort(axis, kind)

def argmin(a, axis=None, fill_value=None):
    """Returns the array of indices for the minimum values of `a` along the 
    specified axis.
    Masked values are treated as if they had the value `fill_value`.        
    If `fill_value` is None, the default for the data type is used.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `fill_value` : var *[None]*
        Default filling value. If None, uses the data type default.
    """
    if fill_value is None:
        fill_value = default_fill_value(a)
    d = filled(a, fill_value)
    if axis is None:
        return d.argmin(axis=None)
    return d.argmin(axis=axis)

def argmax(a, axis=None, fill_value=None):
    """Returns the array of indices for the maximum values of `a` along the 
    specified axis.
    Masked values are treated as if they had the value `fill_value`.        
    If `fill_value` is None, the default for the data type is used.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `fill_value` : var *[None]*
        Default filling value. If None, uses the data type default.
    """
    if fill_value is None:
        fill_value = default_fill_value(a)
        try:
            fill_value = - fill_value
        except:
            pass
    d = filled(a, fill_value)
    if axis is None:
        return d.argmax(axis=None)
    return d.argmax(axis=axis)

def sort(a, axis=-1, kind='quicksort', order=None, endwith=True):
    """
    Sort a along the given axis.

Keyword arguments:

axis  -- axis to be sorted (default -1)
kind  -- sorting algorithm (default 'quicksort')
         Possible values: 'quicksort', 'mergesort', or 'heapsort'.
order -- If a has fields defined, then the order keyword can be the
         field name to sort on or a list (or tuple) of field names
         to indicate the order that fields should be used to define
         the sort.
endwith--Boolean flag indicating whether missing values (if any) should
         be forced in the upper indices (at the end of the array) or
         lower indices (at the beginning).

Returns: None.

This method sorts 'a' in place along the given axis using the algorithm
specified by the kind keyword.

The various sorts may characterized by average speed, worst case
performance, need for work space, and whether they are stable. A stable
sort keeps items with the same key in the same relative order and is most
useful when used with argsort where the key might differ from the items
being sorted. The three available algorithms have the following properties:

|------------------------------------------------------|
|    kind   | speed |  worst case | work space | stable|
|------------------------------------------------------|
|'quicksort'|   1   | O(n^2)      |     0      |   no  |
|'mergesort'|   2   | O(n*log(n)) |    ~n/2    |   yes |
|'heapsort' |   3   | O(n*log(n)) |     0      |   no  |
|------------------------------------------------------|

All the sort algorithms make temporary copies of the data when the sort is
not along the last axis. Consequently, sorts along the last axis are faster
and use less space than sorts along other axis.

"""
    a = numeric.asanyarray(a)
    if endwith:
        filler = minimum_fill_value(a)
    else:
        filler = maximum_fill_value(a)
    indx = filled(a,filler).argsort(axis=axis,kind=kind,order=order)
    return a[indx]   

def compressed(x):
    """Returns a compressed version of a masked array (or just the array if it
    wasn't masked first)."""
    if getmask(x) is None:
        return x
    else:
        return x.compressed()

def count(a, axis = None):
    "Count of the non-masked elements in a, or along a certain axis."
    a = masked_array(a)
    return a.count(axis)

def concatenate(arrays, axis=0):
    "Concatenates the arrays along the given axis"
    #TODO: We lose the subclass, here! We should keep track of the classes...
    #TODO: ...and find the max ? the lowest according to MRO?
    d = []
    for x in arrays:
        d.append(filled(x))
    d = numeric.concatenate(d, axis)
    for x in arrays:
        if getmask(x) is not nomask: 
            break
    else:
        return masked_array(d)
    dm = []
    for x in arrays:
        dm.append(getmaskarray(x))
    dm = make_mask(numeric.concatenate(dm, axis), copy=False, small_mask=True)
    return masked_array(d, mask=dm)

def expand_dims(x,axis):
    """Expand the shape of a by including newaxis before given axis."""
    if isinstance(x, MaskedArray):
        (d,m) = (x._data, x._mask)
        if m is nomask:
            return masked_array(n_expand_dims(d,axis),                                
                                dtype=d.dtype, fill_value=x._fill_value)
        else:
            return masked_array(n_expand_dims(d,axis), 
                                mask=n_expand_dims(m,axis),
                                dtype=d.dtype, fill_value=x._fill_value)
    else:
        return n_expand_dims(x,axis)

#......................................
def left_shift (a, n):
    "Left shift n bits"
    m = getmask(a)
    if m is nomask:
        d = umath.left_shift(filled(a), n)
        return masked_array(d)
    else:
        d = umath.left_shift(filled(a, 0), n)
        return masked_array(d, mask=m)

def right_shift (a, n):
    "Right shift n bits"
    m = getmask(a)
    if m is nomask:
        d = umath.right_shift(filled(a), n)
        return masked_array(d)
    else:
        d = umath.right_shift(filled(a, 0), n)
        return masked_array(d, mask=m)
#......................................
def put(a, indices, values, mode='raise'):
    """Sets storage-indexed locations to corresponding values.
    Values and indices are filled if necessary."""
    # We can't use 'frommethod', the order of arguments is different
    try:
        return a.put(indices, values, mode=mode)
    except AttributeError:
        return fromnumeric.asarray(a).put(indices, values, mode=mode)

def putmask(a, mask, values): #, mode='raise'):
    """`putmask(a, mask, v)` results in `a = v` for all places where `mask` is true.
If `v` is shorter than `mask`, it will be repeated as necessary.
In particular `v` can be a scalar or length 1 array."""
    # We can't use 'frommethod', the order of arguments is different
    try:
        return a.putmask(values, mask)
    except AttributeError:
        return fromnumeric.asarray(a).putmask(values, mask)

def transpose(a,axes=None):
    """Returns a view of the array with dimensions permuted according to axes.  
If `axes` is None (default), returns array with dimensions reversed.
    """
    #We can't use 'frommethod', as 'transpose' doesn't take keywords
    try:
        return a.transpose(axes)
    except AttributeError:
        return fromnumeric.asarray(a).transpose(axes)
    
def reshape(a, new_shape):
    """Changes the shape of the array `a` to `new_shape`."""
    #We can't use 'frommethod', it whine about some parameters. Dmmit.
    try:
        return a.reshape(new_shape)
    except AttributeError:
        return fromnumeric.asarray(a).reshape(new_shape)

def resize(x, new_shape):
    """resize(a,new_shape) returns a new array with the specified shape.
    The total size of the original array can be any size. 
    The new array is filled with repeated copies of a. If a was masked, the new
    array will be masked, and the new mask will be a repetition of the old one.
    """
    # We can't use _frommethods here, as N.resize is notoriously whiny.
    m = getmask(x)
    if m is not nomask:
        m = fromnumeric.resize(m, new_shape)
    if isinstance(x, MaskedArray):
        result = x.__class__(fromnumeric.resize(filled(x), new_shape), mask=m)
    else:
        result = masked_array(fromnumeric.resize(filled(x), new_shape), mask=m)
    result.set_fill_value(get_fill_value(x))
    return result


#................................................
def rank(obj):
    """Gets the rank of sequence a (the number of dimensions, not a matrix rank)
The rank of a scalar is zero."""
    return fromnumeric.rank(filled(obj))
#
def shape(obj):
    """Returns the shape of `a` (as a function call which also works on nested sequences).
    """
    return fromnumeric.shape(filled(obj))
#
def size(obj, axis=None):
    """Returns the number of elements in the array along the given axis, 
or in the sequence if `axis` is None.
    """
    return fromnumeric.size(filled(obj), axis)    
#................................................

#####--------------------------------------------------------------------------
#---- --- Extra functions ---
#####--------------------------------------------------------------------------
def where (condition, x, y):
    """where(condition, x, y) is x where condition is nonzero, y otherwise.
       condition must be convertible to an integer array.
       Answer is always the shape of condition.
       The type depends on x and y. It is integer if both x and y are
       the value masked.
    """
    fc = filled(not_equal(condition, 0), 0)
    xv = filled(x)
    xm = getmask(x)
    yv = filled(y)
    ym = getmask(y)
    d = numeric.choose(fc, (yv, xv))
    md = numeric.choose(fc, (ym, xm))
    m = getmask(condition)
    m = make_mask(mask_or(m, md), copy=False, small_mask=True)
    return masked_array(d, mask=m)

def choose (indices, t, out=None, mode='raise'):
    "Returns array shaped like indices with elements chosen from t"
    #TODO: implement options `out` and `mode`, if possible.
    def fmask (x):
        "Returns the filled array, or True if ``masked``."
        if x is masked: 
            return 1
        return filled(x)
    def nmask (x):
        "Returns the mask, True if ``masked``, False if ``nomask``."
        if x is masked: 
            return 1
        m = getmask(x)
        if m is nomask: 
            return 0
        return m
    c = filled(indices, 0)
    masks = [nmask(x) for x in t]
    a = [fmask(x) for x in t]
    d = numeric.choose(c, a)
    m = numeric.choose(c, masks)
    m = make_mask(mask_or(m, getmask(indices)), copy=0, small_mask=1)
    return masked_array(d, mask=m)

def round_(a, decimals=0, out=None):
    """Returns reference to result. Copies a and rounds to 'decimals' places.

    Keyword arguments:
        decimals -- number of decimals to round to (default 0). May be negative.
        out -- existing array to use for output (default copy of a).

    Return:
        Reference to out, where None specifies a copy of the original array a.

    Round to the specified number of decimals. When 'decimals' is negative it
    specifies the number of positions to the left of the decimal point. The
    real and imaginary parts of complex numbers are rounded separately.
    Nothing is done if the array is not of float type and 'decimals' is greater
    than or equal to 0."""
    if not hasattr(a, "_mask"):
        mask = nomask
    else:
        mask = a._mask
    if out is None:
        return a.__class__(fromnumeric.round_(a, decimals, None), mask=mask)
    else:
        out = a.__class__(fromnumeric.round_(a, decimals, out), mask=mask)
        return out

def arange(start, stop=None, step=1, dtype=None):
    """Just like range() except it returns a array whose type can be specified
    by the keyword argument dtype.
    """
    return array(numeric.arange(start, stop, step, dtype), mask=nomask)
    
def inner(a, b):
    """inner(a,b) returns the dot product of two arrays, which has
    shape a.shape[:-1] + b.shape[:-1] with elements computed by summing the
    product of the elements from the last dimensions of a and b.
    Masked elements are replace by zeros.
    """
    fa = filled(a, 0)
    fb = filled(b, 0)
    if len(fa.shape) == 0: 
        fa.shape = (1,)
    if len(fb.shape) == 0: 
        fb.shape = (1,)
    return masked_array(numeric.inner(fa, fb))
innerproduct = inner

def outer(a, b):
    """outer(a,b) = {a[i]*b[j]}, has shape (len(a),len(b))"""
    fa = filled(a, 0).ravel()
    fb = filled(b, 0).ravel()
    d = numeric.outer(fa, fb)
    ma = getmask(a)
    mb = getmask(b)
    if ma is nomask and mb is nomask:
        return masked_array(d)
    ma = getmaskarray(a)
    mb = getmaskarray(b)
    m = make_mask(1-numeric.outer(1-ma, 1-mb), copy=0)
    return masked_array(d, mask=m)
outerproduct = outer    

def allequal (a, b, fill_value=True):
    """
Returns `True` if all entries of  a and b are equal, using
fill_value as a truth value where either or both are masked.
    """
    m = mask_or(getmask(a), getmask(b))
    if m is nomask:
        x = filled(a)
        y = filled(b)
        d = umath.equal(x, y)
        return d.all()
    elif fill_value:
        x = filled(a)
        y = filled(b)
        d = umath.equal(x, y)
        dm = array(d, mask=m, copy=False)
        return dm.filled(True).all(None)
    else:
        return False
    
def allclose (a, b, fill_value=True, rtol=1.e-5, atol=1.e-8):
    """ Returns `True` if all elements of `a` and `b` are equal subject to given tolerances.
If `fill_value` is True, masked values are considered equal.
If `fill_value` is False, masked values considered unequal.
The relative error rtol should be positive and << 1.0
The absolute error `atol` comes into play for those elements of `b` 
 that are very small or zero; it says how small `a` must be also.
    """
    m = mask_or(getmask(a), getmask(b))
    d1 = filled(a)
    d2 = filled(b)
    x = filled(array(d1, copy=0, mask=m), fill_value).astype(float)
    y = filled(array(d2, copy=0, mask=m), 1).astype(float)
    d = umath.less_equal(umath.absolute(x-y), atol + rtol * umath.absolute(y))
    return fromnumeric.alltrue(fromnumeric.ravel(d))

#..............................................................................
def asarray(a, dtype=None):
    """asarray(data, dtype) = array(data, dtype, copy=0)
Returns `a` as an masked array.
No copy is performed if `a` is already an array.
Subclasses are converted to base class MaskedArray.
    """
    return masked_array(a, dtype=dtype, copy=False, keep_mask=True)

def empty(new_shape, dtype=float):
    """empty((d1,...,dn),dtype=float,order='C')
Returns a new array of shape (d1,...,dn) and given type with all its
entries uninitialized. This can be faster than zeros."""
    return masked_array(numeric.empty(new_shape, dtype), mask=nomask)

def empty_like(a):
    """empty_like(a)
Returns an empty (uninitialized) array of the shape and typecode of a.
Note that this does NOT initialize the returned array.  
If you require your array to be initialized, you should use zeros_like()."""
    return masked_array(numeric.empty_like(a), mask=nomask)

def ones(new_shape, dtype=float):
    """ones(shape, dtype=None) 
Returns an array of the given dimensions, initialized to all ones."""
    return masked_array(numeric.ones(new_shape, dtype), mask=nomask)

def zeros(new_shape, dtype=float):
    """zeros(new_shape, dtype=None) 
Returns an array of the given dimensions, initialized to all zeros."""
    return masked_array(numeric.zeros(new_shape, dtype), mask=nomask)

#####--------------------------------------------------------------------------
#---- --- Pickling ---
#####--------------------------------------------------------------------------
#FIXME: We're kinda stuck with forcing the mask to have the same shape as the data
def _mareconstruct(subtype, baseshape, basetype,):
    """Internal function that builds a new MaskedArray from the information stored
in a pickle."""
    _data = ndarray.__new__(ndarray, baseshape, basetype)
    _mask = ndarray.__new__(ndarray, baseshape, basetype)
    return MaskedArray.__new__(subtype, _data, mask=_mask, dtype=basetype, small_mask=False)

def _getstate(a):
    "Returns the internal state of the masked array, for pickling purposes."
    state = (1,
             a.shape, 
             a.dtype,
             a.flags.fnc,
             (a._data).__reduce__()[-1][-1],
             getmaskarray(a).__reduce__()[-1][-1])
    return state
    
def _setstate(a, state):
    """Restores the internal state of the masked array, for pickling purposes.
`state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:

    - class name
    - a tuple giving the shape of the data
    - a typecode for the data
    - a binary string for the data
    - a binary string for the mask.
        """
    (ver, shp, typ, isf, raw, msk) = state
    (a._data).__setstate__((shp, typ, isf, raw))
    (a._mask).__setstate__((shp, dtype('|b1'), isf, msk))
        
def _reduce(a):
    """Returns a 3-tuple for pickling a MaskedArray."""
    return (_mareconstruct,
            (a.__class__, (0,), 'b', ),
            a.__getstate__())
    
def dump(a,F):
    """Pickles the MaskedArray `a` to the file `F`.
`F` can either be the handle of an exiting file, or a string representing a file name.
    """
    if not hasattr(F,'readline'):
        F = open(F,'w')
    return cPickle.dump(a,F)
    
def dumps(a):
    """Returns a string corresponding to the pickling of the MaskedArray."""
    return cPickle.dumps(a)  

def load(F):
    """Wrapper around ``cPickle.load`` which accepts either a file-like object or
 a filename."""
    if not hasattr(F, 'readline'): 
        F = open(F,'r')
    return cPickle.load(F)

def loads(strg):
    "Loads a pickle from the current string."""
    return cPickle.loads(strg)

MaskedArray.__getstate__ = _getstate
MaskedArray.__setstate__ = _setstate
MaskedArray.__reduce__ = _reduce
MaskedArray.__dump__ = dump
MaskedArray.__dumps__ = dumps

################################################################################

if __name__ == '__main__':
    if 1:
        a = masked_array([0,0])
        b = a/a
        assert (b._mask == [1,1]).all()
        assert (a._data == [0,0]).all()
    if 1:
        a = arange(4)
        a[1:-1] = masked
        b = a[:-5]
    if 1:
        assert(masked_array(masked) is masked)
