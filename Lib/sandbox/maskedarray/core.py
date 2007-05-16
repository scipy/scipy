# pylint: disable-msg=E1002
"""MA: a facility for dealing with missing observations
MA is generally used as a numpy.array look-alike.
by Paul F. Dubois.

Copyright 1999, 2000, 2001 Regents of the University of California.
Released for unlimited redistribution.
Adapted for numpy_core 2005 by Travis Oliphant and
(mainly) Paul Dubois.

Subclassing of the base ndarray 2006 by Pierre Gerard-Marchant.
pgmdevlist_AT_gmail_DOT_com
Improvements suggested by Reggie Dugard (reggie_AT_merfinllc_DOT_com)

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
import numpy.core.numeric as numeric
import numpy.core.numerictypes as ntypes
from numpy import bool_, dtype, typecodes, amax, amin, ndarray
from numpy import expand_dims as n_expand_dims
import warnings


MaskType = bool_
nomask = MaskType(0)

divide_tolerance = 1.e-35
numpy.seterr(all='ignore')

# TODO: There's still a problem with N.add.reduce not working...
# TODO: ...neither does N.add.accumulate

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
max_filler = ntypes._minvals
max_filler.update([(k,-numeric.inf) for k in [numpy.float32, numpy.float64]])
min_filler = ntypes._maxvals
min_filler.update([(k,numeric.inf) for k in [numpy.float32, numpy.float64]])
if 'float128' in ntypes.typeDict:
    max_filler.update([(numpy.float128,-numeric.inf)])
    min_filler.update([(numpy.float128, numeric.inf)])


def default_fill_value(obj):
    "Calculates the default fill value for an object `obj`."
    if hasattr(obj,'dtype'):
        defval = default_filler[obj.dtype.kind]
    elif isinstance(obj, numeric.dtype):
        defval = default_filler[obj.kind]
    elif isinstance(obj, float):
        defval = default_filler['f']
    elif isinstance(obj, int) or isinstance(obj, long):
        defval = default_filler['i']
    elif isinstance(obj, str):
        defval = default_filler['S']
    elif isinstance(obj, complex):
        defval = default_filler['c']
    else:
        defval = default_filler['O']
    return defval

def minimum_fill_value(obj):
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

def maximum_fill_value(obj):
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

def set_fill_value(a, fill_value):
    "Sets the fill value of `a` if it is a masked array."
    if isinstance(a, MaskedArray):
        a.set_fill_value(fill_value)

def get_fill_value(a):
    """Returns the fill value of `a`, if any.
    Otherwise, returns the default fill value for that type.
    """
    if isinstance(a, MaskedArray):
        result = a.fill_value
    else:
        result = default_fill_value(a)
    return result

def common_fill_value(a, b):
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
    elif isinstance(a, dict):
        return numeric.array(a, 'O')
    else:
        return numeric.array(a)

def get_masked_subclass(*arrays):
    """Returns the youngest subclass of MaskedArray from a list of arrays,
 or MaskedArray. In case of siblings, the first takes over."""
    if len(arrays) == 1:
        arr = arrays[0]
        if isinstance(arr, MaskedArray):
            rcls = type(arr)
        else:
            rcls = MaskedArray
    else:
        arrcls = [type(a) for a in arrays]
        rcls = arrcls[0]
        if not issubclass(rcls, MaskedArray):
            rcls = MaskedArray
        for cls in arrcls[1:]:
            if issubclass(cls, rcls):
                rcls = cls
    return rcls

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
        # Take care of the masked singletong first ...
        if m.ndim == 0 and m:
            return masked
        # Get the result....
        if isinstance(a, MaskedArray):
            result = self.f(d1, *args, **kwargs).view(type(a))
        else:
            result = self.f(d1, *args, **kwargs).view(MaskedArray)
        # Fix the mask if we don't have a scalar
        if result.ndim > 0:
            result._mask = m
        return result
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
        if (not m.ndim) and m:
            return masked
        d1 = filled(a, self.fillx)
        d2 = filled(b, self.filly)
# CHECK : Do we really need to fill the arguments ? Pro'ly not        
#        result = self.f(a, b, *args, **kwargs).view(get_masked_subclass(a,b))
        result = self.f(d1, d2, *args, **kwargs).view(get_masked_subclass(a,b))
        if result.ndim > 0:
            result._mask = m
        return result
    #
    def reduce (self, target, axis=0, dtype=None):
        """Reduces `target` along the given `axis`."""
        if isinstance(target, MaskedArray):
            tclass = type(target)
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
            return self.f.reduce(t, axis).view(tclass)
        t = t.view(tclass)
        t._mask = m
        # XXX: "or t.dtype" below is a workaround for what appears
        # XXX: to be a bug in reduce.
        tr = self.f.reduce(filled(t, self.filly), axis, dtype=dtype or t.dtype)
        mr = umath.logical_and.reduce(m, axis)
        tr = tr.view(tclass)
        if mr.ndim > 0:
            tr._mask = mr
            return tr
        elif mr:
            return masked
        return tr

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
        if (not m.ndim) and m:
            return masked
        rcls = get_masked_subclass(a,b)
        d = self.f.outer(filled(a, self.fillx), filled(b, self.filly)).view(rcls)
        if d.ndim > 0:
            d._mask = m
        return d

    def accumulate (self, target, axis=0):
        """Accumulates `target` along `axis` after filling with y fill value."""
        if isinstance(target, MaskedArray):
            tclass = type(target)
        else:
            tclass = masked_array
        t = filled(target, self.filly)
        return self.f.accumulate(t, axis).view(tclass)

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
        if (not m.ndim) and m:
            return masked       
        result =  self.f(d1, d2).view(get_masked_subclass(a,b))
        if result.ndim > 0:
            result._mask = m
        return result

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
def masked_where(condition, a, copy=True):
    """Returns `x` as an array masked where `condition` is true.
Masked values of `x` or `condition` are kept.

:Parameters:
    - `condition` (ndarray) : Masking condition.
    - `x` (ndarray) : Array to mask.
    - `copy` (boolean, *[False]*) : Returns a copy of `m` if true.
    """
    cond = filled(condition,1)
    a = numeric.array(a, copy=copy, subok=True)
    if hasattr(a, '_mask'):
        cond = mask_or(cond, a._mask)
        cls = type(a)
    else:
        cls = MaskedArray
    result = a.view(cls)
    result._mask = cond
    return result

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
##def _getoptions(a_out, a_in):
##    "Copies standards options of a_in to a_out."
##    for att in [']
#class _mathmethod(object):
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
#        self.obj = None
#        self.__doc__ = self.getdoc()
#    #
#    def getdoc(self):
#        "Returns the doc of the function (from the doc of the method)."
#        try:
#            return getattr(MaskedArray, self.methodname).__doc__
#        except:
#            return getattr(ndarray, self.methodname).__doc__
#    #
#    def __get__(self, obj, objtype=None):
#        self.obj = obj
#        return self
#    #
#    def __call__ (self, other, *args):
#        "Execute the call behavior."
#        instance = self.obj
#        m_self = instance._mask
#        m_other = getmask(other)
#        base = instance.filled(self.fill_self)
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
#                    # We don't want to modify other: let's copy target, then
#                    target = target.copy()
#                    target[fromnumeric.asarray(domain)] = self.fill_other
#                else:
#                    target = numeric.where(fromnumeric.asarray(domain),
#                                           self.fill_other, target)
#                m_other = mask_or(m_other, domain)
#        m = mask_or(m_self, m_other)
#        method = getattr(base, self.methodname)
#        result = method(target, *args).view(type(instance))
#        try:
#            result._mask = m
#        except AttributeError:
#            if m:
#                result = masked
#        return result
#...............................................................................
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
        methdoc = getattr(ndarray, self._name, None)
        methdoc = getattr(numpy, self._name, methdoc)
#        methdoc = getattr(MaskedArray, self._name, methdoc)
        if methdoc is not None:
            return methdoc.__doc__
#        try:
#            return getattr(MaskedArray, self._name).__doc__
#        except:
#            try:
#                return getattr(numpy, self._name).__doc__
#            except:
#                return getattr(ndarray, self._name).__doc
    #
    def __get__(self, obj, objtype=None):
        self.obj = obj
        return self
    #
    def __call__(self, *args, **params):
        methodname = self._name
        data = self.obj._data
        mask = self.obj._mask
        cls = type(self.obj)
        result = getattr(data, methodname)(*args, **params).view(cls)
        result._smallmask = self.obj._smallmask
        if result.ndim:
            if not self._onmask:
                result._mask = mask
            elif mask is not nomask:
                result.__setmask__(getattr(mask, methodname)(*args, **params))
        return result
#..........................................................

class flatiter(object):
    "Defines an interator."
    def __init__(self, ma):
        self.ma = ma
        self.ma_iter = numpy.asarray(ma).flat

        if ma._mask is nomask:
            self.maskiter = None
        else:
            self.maskiter = ma._mask.flat

    def __iter__(self):
        return self

    ### This won't work is ravel makes a copy
    def __setitem__(self, index, value):
        a = self.ma.ravel()
        a[index] = value

    def next(self):
        d = self.ma_iter.next()
        if self.maskiter is not None and self.maskiter.next():
            d = masked
        return d


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
    _baseclass =  numeric.ndarray
    def __new__(cls, data=None, mask=nomask, dtype=None, copy=False, fill_value=None,
                keep_mask=True, small_mask=True, hard_mask=False, flag=None,
                subok=True, **options):
        """array(data, dtype=None, copy=True, mask=nomask, fill_value=None)

If `data` is already a ndarray, its dtype becomes the default value of dtype.
        """
        if flag is not None:
            warnings.warn("The flag 'flag' is now called 'small_mask'!",
                          DeprecationWarning)
            small_mask = flag
        # Process data............
        _data = numeric.array(data, dtype=dtype, copy=copy, subok=subok)
        _baseclass = getattr(data, '_baseclass', type(_data))
        _basedict = getattr(data, '_basedict', getattr(data, '__dict__', None))
        if not isinstance(data, MaskedArray): 
            _data = _data.view(cls)
        elif not subok:
            _data = data.view(cls)
        else:
            _data = _data.view(type(data))
        # Backwards compat .......
        if hasattr(data,'_mask') and not isinstance(data, ndarray):
            _data._mask = data._mask
            _sharedmask = True
        # Process mask ...........
        if mask is nomask:
            if not keep_mask:
                _data._mask = nomask
            if copy:
                _data._mask = _data._mask.copy()
        else:
            mask = numeric.array(mask, dtype=MaskType, copy=copy)
            if mask.shape != _data.shape:
                (nd, nm) = (_data.size, mask.size) 
                if nm == 1:
                    mask = numeric.resize(mask, _data.shape)
                elif nm == nd:
                    mask = fromnumeric.reshape(mask, _data.shape)
                else:
                    msg = "Mask and data not compatible: data size is %i, "+\
                          "mask size is %i."
                    raise MAError, msg % (nd, nm)
            if _data._mask is nomask:
                _data._mask = mask
                _data._sharedmask = True
            else:
                # Make a copy of the mask to avoid propagation
                _data._sharedmask = False
                if not keep_mask:
                    _data._mask = mask
                else:
                    _data._mask = umath.logical_or(mask, _data._mask) 
                    
                    
        # Update fill_value.......
        _data._fill_value = getattr(data, '_fill_value', fill_value)
        if _data._fill_value is None:
            _data._fill_value = default_fill_value(_data)
        # Process extra options ..
        _data._hardmask = hard_mask
        _data._smallmask = small_mask
        _data._baseclass = _baseclass
        _data._basedict = _basedict
        return _data
    #........................
    def __array_finalize__(self,obj):
        """Finalizes the masked array.
        """
        # Finalize mask ...............
        self._mask = getattr(obj, '_mask', nomask)
        if self._mask is not nomask:
            self._mask.shape = self.shape
        # Get the remaining options ...
        self._hardmask = getattr(obj, '_hardmask', self._defaulthardmask)
        self._smallmask = getattr(obj, '_smallmask', True)
        self._sharedmask = True
        self._baseclass = getattr(obj, '_baseclass', type(obj))
        self._fill_value = getattr(obj, '_fill_value', None)
        # Update special attributes ...
        self._basedict = getattr(obj, '_basedict', getattr(obj, '__dict__', None))
        if self._basedict is not None:
            self.__dict__.update(self._basedict)
        return
    #..................................
    def __array_wrap__(self, obj, context=None):
        """Special hook for ufuncs.
Wraps the numpy array and sets the mask according to context.
        """
        #TODO : Should we check for type result 
        result = obj.view(type(self))
        #..........
        if context is not None:
            result._mask = result._mask.copy()
            (func, args, _) = context
            m = reduce(mask_or, [getmask(arg) for arg in args])
            # Get domain mask
            domain = ufunc_domain.get(func, None)
            if domain is not None:
                if len(args) > 2:
                    d = reduce(domain, args)
                else:
                    d = domain(*args)
                if m is nomask:
                    if d is not nomask:
                        m = d
                else:
                    m |= d
            if not m.ndim and m:
                if m:
                    if result.shape == ():
                        return masked
                    result._mask = numeric.ones(result.shape, bool_)
            else:
                result._mask = m
        #....
#        result._mask = m
        result._fill_value = self._fill_value
        result._hardmask = self._hardmask
        result._smallmask = self._smallmask
        result._baseclass = self._baseclass
        return result
    #.............................................
    def __getitem__(self, indx):
        """x.__getitem__(y) <==> x[y]
Returns the item described by i. Not a copy as in previous versions.
        """
        # This test is useful, but we should keep things light...
#        if getmask(indx) is not nomask:
#            msg = "Masked arrays must be filled before they can be used as indices!"
#            raise IndexError, msg
        # super() can't work here if the underlying data is a matrix...
        dout = (self._data).__getitem__(indx)
        m = self._mask
        if hasattr(dout, 'shape') and len(dout.shape) > 0:
            # Not a scalar: make sure that dout is a MA
            dout = dout.view(type(self))
            dout._smallmask = self._smallmask
            if m is not nomask:
                # use _set_mask to take care of the shape
                dout.__setmask__(m[indx])
        elif m is not nomask and m[indx]:
            return masked
        return dout
    #........................
    def __setitem__(self, indx, value):
        """x.__setitem__(i, y) <==> x[i]=y
Sets item described by index. If value is masked, masks those locations.
        """
        if self is masked:
            raise MAError, 'Cannot alter the masked element.'
#        if getmask(indx) is not nomask:
#            msg = "Masked arrays must be filled before they can be used as indices!"
#            raise IndexError, msg
        #....
        if value is masked:
            if self._mask is nomask:
                self._mask = make_mask_none(self.shape)
            else:
                self._mask = self._mask.copy()
            self._mask[indx] = True
            return
        #....
        dval = numeric.asarray(value).astype(self.dtype)
        valmask = getmask(value)
        if self._mask is nomask:
            if valmask is not nomask:
                self._mask = make_mask_none(self.shape)
                self._mask[indx] = valmask
        elif not self._hardmask:
            _mask = self._mask.copy()
            if valmask is nomask:
                _mask[indx] = False
            else:
                _mask[indx] = valmask
            self._set_mask(_mask)
        elif hasattr(indx, 'dtype') and (indx.dtype==bool_):
            indx = indx * umath.logical_not(self._mask)
        else:
            mindx = mask_or(self._mask[indx], valmask, copy=True)
            dindx = self._data[indx]
            if dindx.size > 1:
                dindx[~mindx] = dval
            elif mindx is nomask:
                dindx = dval
            dval = dindx
            self._mask[indx] = mindx
        # Set data ..........
        #dval = filled(value).astype(self.dtype)
        ndarray.__setitem__(self._data,indx,dval)
    #............................................
    def __getslice__(self, i, j):
        """x.__getslice__(i, j) <==> x[i:j]
Returns the slice described by i, j.
The use of negative indices is not supported."""
        return self.__getitem__(slice(i,j))
    #........................
    def __setslice__(self, i, j, value):
        """x.__setslice__(i, j, value) <==> x[i:j]=value
Sets a slice i:j to `value`.
If `value` is masked, masks those locations."""
        self.__setitem__(slice(i,j), value)
    #............................................
    def __setmask__(self, mask, copy=False):
        newmask = make_mask(mask, copy=copy, small_mask=self._smallmask)
#        self.unshare_mask()
        if self._mask is nomask:
            self._mask = newmask
        elif self._hardmask:
            if newmask is not nomask:
                self._mask.__ior__(newmask)
        else:
            # This one is tricky: if we set the mask that way, we may break the
            # propagation. But if we don't, we end up with a mask full of False
            # and a test on nomask fails...
            if newmask is nomask:
                self._mask = nomask
            else:
                self._mask.flat = newmask
        if self._mask.shape:
            self._mask = numeric.reshape(self._mask, self.shape)
    _set_mask = __setmask__
    
    def _get_mask(self):
        """Returns the current mask."""
        return self._mask

    mask = property(fget=_get_mask, fset=__setmask__, doc="Mask")
    #............................................
    def harden_mask(self):
        "Forces the mask to hard."
        self._hardmask = True
        
    def soften_mask(self):
        "Forces the mask to soft."
        self._hardmask = False     
        
    def unshare_mask(self):
        "Copies the mask and set the sharedmask flag to False."
        if self._sharedmask:
            self._mask = self._mask.copy()
            self._sharedmask = False
        
    #............................................
    def _get_data(self):
        "Returns the current data (as a view of the original underlying data)>"
        return self.view(self._baseclass)
    _data = property(fget=_get_data)        
    #............................................
    def _get_flat(self):
        """Calculates the flat value.
        """
        return flatiter(self)
    #
    def _set_flat (self, value):
        "x.flat = value"
        y = self.ravel()
        y[:] = value
    #
    flat = property(fget=_get_flat, fset=_set_flat, doc="Flat version")
    #............................................
    def get_fill_value(self):
        "Returns the filling value."
        if self._fill_value is None:
            self._fill_value = default_fill_value(self)
        return self._fill_value

    def set_fill_value(self, value=None):
        """Sets the filling value to `value`.
If None, uses the default, based on the data type."""
        if value is None:
            value = default_fill_value(self)
        self._fill_value = value

    fill_value = property(fget=get_fill_value, fset=set_fill_value,
                          doc="Filling value")

    def filled(self, fill_value=None):
        """Returns an array of the same class as `_data`,
 with masked values filled with `fill_value`.
Subclassing is preserved.

If `fill_value` is None, uses self.fill_value.
        """
        m = self._mask
        if m is nomask:
            return self._data
        #
        if fill_value is None:
            fill_value = self.fill_value
        #
        if self is masked_singleton:
            result = numeric.asanyarray(fill_value)
        else:
            result = self._data.copy()
            try:
                result[m] = fill_value
            except (TypeError, AttributeError):
                fill_value = numeric.array(fill_value, dtype=object)
                d = result.astype(object)
                result = fromnumeric.choose(m, (d, fill_value))
            except IndexError:
                #ok, if scalar
                if self._data.shape:
                    raise
                elif m:
                    result = numeric.array(fill_value, dtype=self.dtype)
                else:
                    result = self._data
        return result

    def compressed(self):
        "A 1-D array of all the non-masked data."
        d = self.ravel()
        if self._mask is nomask:
            return d
        elif not self._smallmask and not self._mask.any():
            return d
        else:
            return d[numeric.logical_not(d._mask)]
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
                if m.shape == ():
                    if m:
                        return str(f)
                    else:
                        return str(self._data)
                # convert to object array to make filled work
#CHECK: the two lines below seem more robust than the self._data.astype
#                res = numeric.empty(self._data.shape, object_)
#                numeric.putmask(res,~m,self._data)
                res = self._data.astype("|O8")
                res[m] = f
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
    def __iadd__(self, other):
        "Adds other to self in place."
        ndarray.__iadd__(self._data,other)
        m = getmask(other)
        if self._mask is nomask:
            self._mask = m
        elif m is not nomask:
            self._mask += m
        return self
    #....
    def __isub__(self, other):
        "Subtracts other from self in place."
        ndarray.__isub__(self._data,other)
        m = getmask(other)
        if self._mask is nomask:
            self._mask = m
        elif m is not nomask:
            self._mask += m
        return self
    #....
    def __imul__(self, other):
        "Multiplies self by other in place."
        ndarray.__imul__(self._data,other)
        m = getmask(other)
        if self._mask is nomask:
            self._mask = m
        elif m is not nomask:
            self._mask += m
        return self
    #....
    def __idiv__(self, other):
        "Divides self by other in place."
        dom_mask = domain_safe_divide().__call__(self, filled(other,1))
        other_mask = getmask(other)
        new_mask = mask_or(other_mask, dom_mask)
        ndarray.__idiv__(self._data, other)
        self._mask = mask_or(self._mask, new_mask)
        return self
    #............................................
    def __float__(self):
        "Converts self to float."
        if self._mask is not nomask:
            warnings.warn("Warning: converting a masked element to nan.")
            return numpy.nan
            #raise MAError, 'Cannot convert masked element to a Python float.'
        return float(self.item())

    def __int__(self):
        "Converts self to int."
        if self._mask is not nomask:
            raise MAError, 'Cannot convert masked element to a Python int.'
        return int(self.item())
    #............................................
    def count(self, axis=None):
        """Counts the non-masked elements of the array along a given axis,
and returns a masked array where the mask is True where all data are masked.
If `axis` is None, counts all the non-masked elements, and returns either a
scalar or the masked singleton."""
        m = self._mask
        s = self.shape
        ls = len(s)
        if m is nomask:
            if ls == 0:
                return 1
            if ls == 1:
                return s[0]
            if axis is None:
                return self.size
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
    def reshape (self, *s):
        """Reshapes the array to shape s.
Returns a new masked array.
If you want to modify the shape in place, please use `a.shape = s`"""
        result = self._data.reshape(*s).view(type(self))
        result.__dict__.update(self.__dict__)
        if result._mask is not nomask:
            result._mask = self._mask.copy()
            result._mask.shape = result.shape
        return result
    #
    repeat = _arraymethod('repeat')
    #
    def resize(self, newshape, refcheck=True, order=False):
        """Attempts to modify size and shape of self inplace.
        The array must own its own memory and not be referenced by other arrays.
        Returns None.
        """
        try:
            self._data.resize(newshape, refcheck, order)
            if self.mask is not nomask:
                self._mask.resize(newshape, refcheck, order)
        except ValueError:
            raise ValueError("Cannot resize an array that has been referenced "
                             "or is referencing another array in this way.\n"
                             "Use the resize function.")
        return None
    #
    flatten = _arraymethod('flatten')
    #
    def put(self, indices, values, mode='raise'):
        """Sets storage-indexed locations to corresponding values.
a.put(values, indices, mode) sets a.flat[n] = values[n] for each n in indices.
`values` can be scalar or an array shorter than indices, and it will be repeated,
if necessary.
If `values` has some masked values, the initial mask is updated in consequence,
else the corresponding values are unmasked.
        """
        m = self._mask
        # Hard mask: Get rid of the values/indices that fall on masked data
        if self._hardmask and self._mask is not nomask:
            mask = self._mask[indices]
            indices = numeric.asarray(indices)
            values = numeric.asanyarray(values)
            values.resize(indices.shape)
            indices = indices[~mask]
            values = values[~mask]
        #....
        self._data.put(indices, values, mode=mode)
        #....
        if m is nomask:
            m = getmask(values)
        else:
            m = m.copy()
            if getmask(values) is nomask:
                m.put(indices, False, mode=mode)
            else:
                m.put(indices, values._mask, mode=mode)
            m = make_mask(m, copy=False, small_mask=True)
        self._mask = m
    #............................................
    def ids (self):
        """Return the address of the data and mask areas."""
        return (self.ctypes.data, self._mask.ctypes.data)    
    #............................................
    def all(self, axis=None, out=None):
        """a.all(axis) returns True if all entries along the axis are True.
    Returns False otherwise. If axis is None, uses the flatten array.
    Masked data are considered as True during computation.
    Outputs a masked array, where the mask is True if all data are masked along the axis.
    Note: the out argument is not really operational...
        """
        d = self.filled(True).all(axis=axis, out=out).view(type(self))
        if d.ndim > 0:
            d.__setmask__(self._mask.all(axis))
        return d

    def any(self, axis=None, out=None):
        """a.any(axis) returns True if some or all entries along the axis are True.
    Returns False otherwise. If axis is None, uses the flatten array.
    Masked data are considered as False during computation.
    Outputs a masked array, where the mask is True if all data are masked along the axis.
    Note: the out argument is not really operational...
        """
        d = self.filled(False).any(axis=axis, out=out).view(type(self))
        if d.ndim > 0:
            d.__setmask__(self._mask.all(axis))
        return d
    
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
        # TODO: What are we doing with `out`?
        m = self._mask
        if m is nomask:
            result = super(MaskedArray, self).trace(offset=offset, axis1=axis1,
                                                    axis2=axis2, out=out)
            return result.astype(dtype)
        else:
            D = self.diagonal(offset=offset, axis1=axis1, axis2=axis2)
            return D.astype(dtype).sum(axis=None)
    #............................................
    def sum(self, axis=None, dtype=None):
        """a.sum(axis=None, dtype=None)
Sums the array `a` over the given axis `axis`.
Masked values are set to 0.
If `axis` is None, applies to a flattened version of the array.
    """
        if self._mask is nomask:
            mask = nomask
        else:
            mask = self._mask.all(axis)
            if (not mask.ndim) and mask:
                return masked
        result = self.filled(0).sum(axis, dtype=dtype).view(type(self))
        if result.ndim > 0:
            result.__setmask__(mask)
        return result

    def cumsum(self, axis=None, dtype=None):
        """a.cumprod(axis=None, dtype=None)
Returns the cumulative sum of the elements of array `a` along the given axis `axis`.
Masked values are set to 0.
If `axis` is None, applies to a flattened version of the array.
        """
        result = self.filled(0).cumsum(axis=axis, dtype=dtype).view(type(self))
        result.__setmask__(self.mask)
        return result

    def prod(self, axis=None, dtype=None):
        """a.prod(axis=None, dtype=None)
Returns the product of the elements of array `a` along the given axis `axis`.
Masked elements are set to 1.
If `axis` is None, applies to a flattened version of the array.
        """
        if self._mask is nomask:
            mask = nomask
        else:
            mask = self._mask.all(axis)
            if (not mask.ndim) and mask:
                return masked
        result = self.filled(1).prod(axis=axis, dtype=dtype).view(type(self))
        if result.ndim:
            result.__setmask__(mask)
        return result
    product = prod

    def cumprod(self, axis=None, dtype=None):
        """a.cumprod(axis=None, dtype=None)
Returns the cumulative product of ethe lements of array `a` along the given axis `axis`.
Masked values are set to 1.
If `axis` is None, applies to a flattened version of the array.
        """
        result = self.filled(1).cumprod(axis=axis, dtype=dtype).view(type(self))
        result.__setmask__(self.mask)
        return result

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
            return super(MaskedArray, self).mean(axis=axis, dtype=dtype)
        else:
            dsum = self.sum(axis=axis, dtype=dtype)
            cnt = self.count(axis=axis)
            return dsum*1./cnt

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
            # TODO: Do we keep super, or var _data and take a view ?
            return super(MaskedArray, self).var(axis=axis, dtype=dtype)
        else:
            cnt = self.count(axis=axis)
            danom = self.anom(axis=axis, dtype=dtype)
            danom *= danom
            dvar = numeric.array(danom.sum(axis) / cnt).view(type(self))
            if axis is not None:
                dvar._mask = mask_or(self._mask.all(axis), (cnt==1))
            return dvar

    def std(self, axis=None, dtype=None):
        """a.std(axis=None, dtype=None)
Returns the standard deviation, a measure of the spread of a distribution.

The standard deviation is the square root of the average of the squared
deviations from the mean, i.e. std = sqrt(mean((x - x.mean())**2)).
        """
        dvar = self.var(axis,dtype)
        if axis is not None or dvar is not masked:
            dvar = sqrt(dvar)
        return dvar
    #............................................
    def argsort(self, axis=None, fill_value=None, kind='quicksort',
                order=None):
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
            fill_value = default_fill_value(self)
        d = self.filled(fill_value).view(ndarray)
        return d.argsort(axis=axis, kind=kind, order=order)
    #........................
    def argmin(self, axis=None, fill_value=None):
        """Returns a ndarray of indices for the minimum values of `a` along the
    specified axis.
    Masked values are treated as if they had the value `fill_value`.
    If `fill_value` is None, the default for the data type is used.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `fill_value` : var *[None]*
        Default filling value. If None, uses the minimum default for the data type.
        """
        if fill_value is None:
            fill_value = minimum_fill_value(self)
        d = self.filled(fill_value)
        return d.argmin(axis)
    #........................
    def argmax(self, axis=None, fill_value=None):
        """Returns the array of indices for the maximum values of `a` along the
    specified axis.
    Masked values are treated as if they had the value `fill_value`.
    If `fill_value` is None, the maximum default for the data type is used.
    Returns a numpy array.

:Keywords:
    `axis` : Integer *[None]*
        Axis to be indirectly sorted (default -1)
    `fill_value` : var *[None]*
        Default filling value. If None, uses the data type default.
        """
        if fill_value is None:
            fill_value = maximum_fill_value(self._data)
        d = self.filled(fill_value)
        return d.argmax(axis)

    def sort(self, axis=-1, kind='quicksort', order=None, 
             endwith=True, fill_value=None):
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

    """
        if self._mask is nomask:
            ndarray.sort(self,axis=axis, kind=kind, order=order)
        else:
            if fill_value is None:
                if endwith:
                    filler = minimum_fill_value(self)
                else:
                    filler = maximum_fill_value(self)
            else:
                filler = fill_value
            idx = numpy.indices(self.shape)
            idx[axis] = self.filled(filler).argsort(axis=axis,kind=kind,order=order)
            idx_l = idx.tolist()
            tmp_mask = self._mask[idx_l].flat
            tmp_data = self._data[idx_l].flat
            self.flat = tmp_data
            self._mask.flat = tmp_mask
        return
    #............................................
    def min(self, axis=None, fill_value=None):
        """Returns the minimum/a along the given axis.
If `axis` is None, applies to the flattened array. Masked values are filled 
with `fill_value` during processing. If `fill_value is None, it is set to the
maximum_fill_value corresponding to the data type."""
        mask = self._mask
        # Check all/nothing case ......
        if mask is nomask:
            return super(MaskedArray, self).min(axis=axis)
        elif (not mask.ndim) and mask:
            return masked
        # Get the mask ................
        if axis is None:
            mask = umath.logical_and.reduce(mask.flat)
        else:
            mask = umath.logical_and.reduce(mask, axis=axis)
        # Get the fil value ...........
        if fill_value is None:
            fill_value = minimum_fill_value(self)
        # Get the data ................
        result = self.filled(fill_value).min(axis=axis).view(type(self))
        if result.ndim > 0:
            result._mask = mask
        return result
    #........................
    def max(self, axis=None, fill_value=None):
        """Returns the maximum/a along the given axis.
If `axis` is None, applies to the flattened array. Masked values are filled 
with `fill_value` during processing. If `fill_value is None, it is set to the
maximum_fill_value corresponding to the data type."""
        mask = self._mask
        # Check all/nothing case ......
        if mask is nomask:
            return super(MaskedArray, self).max(axis=axis)
        elif (not mask.ndim) and mask:
            return masked
        # Check the mask ..............
        if axis is None:
            mask = umath.logical_and.reduce(mask.flat)
        else:
            mask = umath.logical_and.reduce(mask, axis=axis)
        # Get the fill value ..........
        if fill_value is None:
            fill_value = maximum_fill_value(self)
        # Get the data ................
        result = self.filled(fill_value).max(axis=axis).view(type(self))
        if result.ndim > 0:
            result._mask = mask
        return result
    #........................
    def ptp(self, axis=None, fill_value=None):
        """Returns the visible data range (max-min) along the given axis.
If the axis is `None`, applies on a flattened array. Masked values are filled
with `fill_value` for processing. If `fill_value` is None, the maximum is uses
the maximum default, the minimum uses the minimum default."""
        return self.max(axis, fill_value) - self.min(axis, fill_value)

    # Array methods ---------------------------------------
    conj = conjugate = _arraymethod('conjugate')
    copy = _arraymethod('copy')
    diagonal = _arraymethod('diagonal')
    take = _arraymethod('take')
    ravel = _arraymethod('ravel')
    transpose = _arraymethod('transpose')
    T = property(fget=lambda self:self.transpose())
    swapaxes = _arraymethod('swapaxes')
    clip = _arraymethod('clip', onmask=False)
    compress = _arraymethod('compress')
    copy = _arraymethod('copy')
    squeeze = _arraymethod('squeeze')
    #--------------------------------------------
    def tolist(self, fill_value=None):
        """Copies the data portion of the array to a hierarchical python list and
    returns that list. Data items are converted to the nearest compatible Python 
    type. Masked values are filled with `fill_value`"""
        return self.filled(fill_value).tolist()
    #........................
    def tostring(self, fill_value=None):
        """a.tostring(order='C', fill_value=None) -> raw copy of array data as a Python string.

    Keyword arguments:
        order      : order of the data item in the copy {"C","F","A"} (default "C")
        fill_value : value used in lieu of missing data 

    Construct a Python string containing the raw bytes in the array. The order
    of the data in arrays with ndim > 1 is specified by the 'order' keyword and
    this keyword overrides the order of the array. The
    choices are:

        "C"       -- C order (row major)
        "Fortran" -- Fortran order (column major)
        "Any"     -- Current order of array.
        None      -- Same as "Any"
    
    Masked data are filled with fill_value. If fill_value is None, the data-type-
    dependent default is used."""
        return self.filled(fill_value).tostring()   
    #--------------------------------------------
    # Backwards Compatibility. Heck...
    @property
    def data(self):
        """Returns the `_data` part of the MaskedArray."""
        return self._data
    def raw_data(self):
        """Returns the `_data` part of the MaskedArray.
You should really use `data` instead..."""
        return self._data
    #--------------------------------------------
    # Pickling
    def __getstate__(self):
        "Returns the internal state of the masked array, for pickling purposes."
        state = (1,
                 self.shape,
                 self.dtype,
                 self.flags.fnc,
                 self._data.tostring(),
                 getmaskarray(self).tostring(),
                 self._fill_value,
                 )
        return state    
    #
    def __setstate__(self, state):
        """Restores the internal state of the masked array, for pickling purposes.
    `state` is typically the output of the ``__getstate__`` output, and is a 5-tuple:
    
        - class name
        - a tuple giving the shape of the data
        - a typecode for the data
        - a binary string for the data
        - a binary string for the mask.
            """
        (ver, shp, typ, isf, raw, msk, flv) = state
        ndarray.__setstate__(self, (shp, typ, isf, raw))
        self._mask.__setstate__((shp, dtype(bool), isf, msk))
        self.fill_value = flv
    #
    def __reduce__(self):
        """Returns a 3-tuple for pickling a MaskedArray."""
        return (_mareconstruct,
                (self.__class__, self._baseclass, (0,), 'b', ),
                self.__getstate__())
    
    
def _mareconstruct(subtype, baseclass, baseshape, basetype,):
    """Internal function that builds a new MaskedArray from the information stored
in a pickle."""
    _data = ndarray.__new__(baseclass, baseshape, basetype)
    _mask = ndarray.__new__(ndarray, baseshape, 'b1')
    return subtype.__new__(subtype, _data, mask=_mask, dtype=basetype, small_mask=False)
#MaskedArray.__dump__ = dump
#MaskedArray.__dumps__ = dumps
    
    

#####--------------------------------------------------------------------------
#---- --- Shortcuts ---
#####---------------------------------------------------------------------------
def isMaskedArray(x):
    "Is x a masked array, that is, an instance of MaskedArray?"
    return isinstance(x, MaskedArray)
isarray = isMaskedArray
isMA = isMaskedArray  #backward compatibility
#masked = MaskedArray(0, int, mask=1)
masked_singleton = MaskedArray(0, dtype=int_, mask=True)
masked = masked_singleton

masked_array = MaskedArray
def array(data, dtype=None, copy=False, order=False, mask=nomask, subok=True,
          keep_mask=True, small_mask=True, hard_mask=None, fill_value=None):
    """array(data, dtype=None, copy=True, order=False, mask=nomask,
             keep_mask=True, small_mask=True, fill_value=None)
Acts as shortcut to MaskedArray, with options in a different order for convenience.
And backwards compatibility...
    """
    #TODO: we should try to put 'order' somwehere
    return MaskedArray(data, mask=mask, dtype=dtype, copy=copy, subok=subok,
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


#####---------------------------------------------------------------------------
#---- --- Extrema functions ---
#####---------------------------------------------------------------------------
class _extrema_operation(object):
    "Generic class for maximum/minimum functions."
    def __call__(self, a, b=None):
        "Executes the call behavior."
        if b is None:
            return self.reduce(a)
        return where(self.compare(a, b), a, b)
    #.........
    def reduce(self, target, axis=None):
        """Reduces target along the given axis."""
        m = getmask(target)
        if axis is not None:
            kargs = { 'axis' : axis }
        else:
            kargs = {}
            target = target.ravel()

        if m is nomask:
            t = self.ufunc.reduce(target, **kargs)
        else:
            target = target.filled(self.fill_value_func(target)).view(type(target))
            t = self.ufunc.reduce(target, **kargs)
            m = umath.logical_and.reduce(m, **kargs)
            if hasattr(t, '_mask'):
                t._mask = m
            elif m:
                t = masked
        return t
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
        result = self.ufunc.outer(filled(a), filled(b))
        result._mask = m
        return result
#............................
class _minimum_operation(_extrema_operation):
    "Object to calculate minima"
    def __init__ (self):
        """minimum(a, b) or minimum(a)
In one argument case, returns the scalar minimum.
        """
        self.ufunc = umath.minimum
        self.afunc = amin
        self.compare = less
        self.fill_value_func = minimum_fill_value
#............................
class _maximum_operation(_extrema_operation):
    "Object to calculate maxima"
    def __init__ (self):
        """maximum(a, b) or maximum(a)
           In one argument case returns the scalar maximum.
        """
        self.ufunc = umath.maximum
        self.afunc = amax
        self.compare = greater
        self.fill_value_func = maximum_fill_value
#..........................................................
def min(array, axis=None, out=None):
    """Returns the minima along the given axis.
If `axis` is None, applies to the flattened array."""
    if out is not None:
        raise TypeError("Output arrays Unsupported for masked arrays")
    if axis is None:
        return minimum(array)
    else:
        return minimum.reduce(array, axis)
#............................
def max(obj, axis=None, out=None):
    """Returns the maxima along the given axis.
If `axis` is None, applies to the flattened array."""
    if out is not None:
        raise TypeError("Output arrays Unsupported for masked arrays")
    if axis is None:
        return maximum(obj)
    else:
        return maximum.reduce(obj, axis)
#.............................
def ptp(obj, axis=None):
    """a.ptp(axis=None) =  a.max(axis)-a.min(axis)"""
    try:
        return obj.max(axis)-obj.min(axis)
    except AttributeError:
        return max(obj, axis=axis) - min(obj, axis=axis)


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
    def __call__(self, a, *args, **params):
        if isinstance(a, MaskedArray):
            return getattr(a, self._methodname).__call__(*args, **params)
        #FIXME ----
        #As x is not a MaskedArray, we transform it to a ndarray with asarray
        #... and call the corresponding method.
        #Except that sometimes it doesn't work (try reshape([1,2,3,4],(2,2)))
        #we end up with a "SystemError: NULL result without error in PyObject_Call"
        #A dirty trick is then to call the initial numpy function...
        method = getattr(fromnumeric.asarray(a), self._methodname)
        try:
            return method(*args, **params)
        except SystemError:
            return getattr(numpy,self._methodname).__call__(a, *args, **params)

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
def argsort(a, axis=None, kind='quicksort', order=None, fill_value=None):
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
        return d.argsort(kind=kind, order=order)
    return d.argsort(axis, kind=kind, order=order)

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
    return d.argmax(axis=axis)

def sort(a, axis=-1, kind='quicksort', order=None, endwith=True, fill_value=None):
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
    if fill_value is None:
        if endwith:
            filler = minimum_fill_value(a)
        else:
            filler = maximum_fill_value(a)
    else:
        filler = fill_value
#    return
    indx = numpy.indices(a.shape).tolist()
    indx[axis] = filled(a,filler).argsort(axis=axis,kind=kind,order=order)
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
    d = numeric.concatenate([filled(a) for a in arrays], axis)
    rcls = get_masked_subclass(*arrays)
    data = d.view(rcls)
    for x in arrays:
        if getmask(x) is not nomask:
            break
    else:
        return data
    dm = numeric.concatenate([getmaskarray(a) for a in arrays], axis)
    dm = make_mask(dm, copy=False, small_mask=True)
    data._mask = dm
    return data

def expand_dims(x,axis):
    """Expand the shape of a by including newaxis before given axis."""
    result = n_expand_dims(x,axis)
    if isinstance(x, MaskedArray):
        new_shape = result.shape
        result = x.view()
        result.shape = new_shape
        if result._mask is not nomask:
            result._mask.shape = new_shape
    return result

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
    result = fromnumeric.resize(x, new_shape).view(get_masked_subclass(x))
    if result.ndim:
        result._mask = m
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
    result = fromnumeric.round_(filled(a), decimals, out)
    if isinstance(a,MaskedArray):
        result = result.view(type(a))
        result._mask = a._mask
    else:
        result = result.view(MaskedArray)
    return result

def arange(start, stop=None, step=1, dtype=None):
    """Just like range() except it returns a array whose type can be specified
    by the keyword argument dtype.
    """
    return array(numeric.arange(start, stop, step, dtype),mask=nomask)

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
    return numeric.empty(new_shape, dtype).view(MaskedArray)

def empty_like(a):
    """empty_like(a)
Returns an empty (uninitialized) array of the shape and typecode of a.
Note that this does NOT initialize the returned array.
If you require your array to be initialized, you should use zeros_like()."""
    return numeric.empty_like(a).view(MaskedArray)

def ones(new_shape, dtype=float):
    """ones(shape, dtype=None)
Returns an array of the given dimensions, initialized to all ones."""
    return numeric.ones(new_shape, dtype).view(MaskedArray)

def zeros(new_shape, dtype=float):
    """zeros(new_shape, dtype=None)
Returns an array of the given dimensions, initialized to all zeros."""
    return numeric.zeros(new_shape, dtype).view(MaskedArray)

#####--------------------------------------------------------------------------
#---- --- Pickling ---
#####--------------------------------------------------------------------------
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


################################################################################

if __name__ == '__main__':
    if 1:
        x = arange(10)
        assert(x.ctypes.data == x.filled().ctypes.data)
    if 1:
        a = array([1,2,3,4],mask=[0,0,0,0],small_mask=False)
        assert(a.ravel()._mask, [0,0,0,0])
        assert(a.compressed(), a)
        a[0] = masked
        assert(a.compressed()._mask, [0,0,0])