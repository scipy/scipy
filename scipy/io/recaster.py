# Author: Matthew Brett

"""
Recaster class for recasting numeric arrays
"""

from numpy import *
from numpy.lib.utils import deprecate

# deprecated in 0.8, will be removed in 0.9.
@deprecate
def sctype_attributes():
    """Return dictionary describing numpy scalar types

    .. deprecated:: sctype_attributes is deprecated in scipy 0.8 and
                    will be removed in scipy 0.9.
    """
    return _sctype_attributes()


def _sctype_attributes():
    d_dict = {}
    for sc_type in ('complex','float'):
        t_list = sctypes[sc_type]
        for T in t_list:
            F = finfo(T)
            dt = dtype(T)
            d_dict[T] = {
                'kind': dt.kind,
                'size': dt.itemsize,
                'max': F.max,
                'min': F.min}
    for T in sctypes['int']:
        dt = dtype(T)
        sz = dt.itemsize
        bits = sz*8-1
        end = 2**bits
        d_dict[T] = {
            'kind': dt.kind,
            'size': sz,
            'min': -end,
            'max': end-1
            }
    for T in sctypes['uint']:
        dt = dtype(T)
        sz = dt.itemsize
        bits = sz*8
        end = 2**bits
        d_dict[T] = {
            'kind': dt.kind,
            'size': sz,
            'min': 0,
            'max': end
        }
    return d_dict

class RecastError(ValueError):
    pass

# deprecated in 0.8, will be removed in 0.9.
class Recaster(object):
    ''' Class to recast arrays to one of acceptable scalar types

    .. deprecated:: Recaster is deprecated in scipy 0.8 and will be
                    removed in scipy 0.9.

    Initialization specifies acceptable types (ATs)

    Implements recast method - returns array that may be of different
    storage type to the input array, where the new type is one of the
    ATs. Recast method will return a larger type if no smaller type
    will contain the data without loss of precision greater than
    specified in options at object creation.
    '''

    _sctype_attributes = _sctype_attributes()
    _k = 2**10
    _option_defaults = {
        'only_if_none': {
        'fp_to_int': 'if_none',
        'fp_to_fp': 'if_none',
        'int_to_int': 'if_none',
        'int_to_fp': 'if_none',
        'downcast_only': False,
        'downcast_within_fp': False,
        'guarantee_fp_to_fp_precision': False,
        'prefer_input_at_threshold': 0,
        'prefer_int_type': 'i',
        },
        'smallest': {
        'fp_to_int': 'always',
        'fp_to_fp': 'always',
        'int_to_int': 'always',
        'int_to_fp': 'always',
        'downcast_only': False,
        'downcast_within_fp': True,
        'guarantee_fp_to_fp_precision': False,
        'prefer_input_at_threshold': 0,
        'prefer_int_type': 'i',
        },
        'fairly_small': {
        'fp_to_int': 'always',
        'fp_to_fp': 'if_none',
        'int_to_int': 'always',
        'int_to_fp': 'if_none',
        'downcast_only': False,
        'downcast_within_fp': False,
        'guarantee_fp_to_fp_precision': False,
        'prefer_input_at_threshold': 2 * _k,
        'prefer_int_type': 'i',
        },
        'preserve_precision': {
        'fp_to_int': 'never',
        'fp_to_fp': 'if_none',
        'int_to_int': 'if_none',
        'int_to_fp': 'never',
        'downcast_only': False,
        'downcast_within_fp': False,
        'guarantee_fp_to_fp_precision': True,
        'prefer_input_at_threshold': 0,
        'prefer_int_type': 'i',
        }
        }

    @deprecate
    def __init__(self, sctype_list=None,
                 sctype_tols=None,
                 recast_options='only_if_none'):
        ''' Set types for which we are attempting to downcast

        Input
        sctype_list  - list of acceptable scalar types
                     If None defaults to all system types
        sctype_tols  - dictionary key datatype, values rtol, tol
                     to specify tolerances for checking near equality in
                     downcasting. Note that tolerance values for integers
                     are used for upcasting integers to floats
        recast_options - dictionary of options for recasting or string
                     specifying one of default options dictionaries.

        recast_option strings can be:
        only_if_none - only attempts recast if the type is not in
                       acceptable types
        smallest     - return array of smallest possible type within tolerance
        fairly_small - compromise set of options between speed of downcast and
                       size of output
        preserve_precision - recasts arrays only to types that preserve precision

        Elements in recast_options dictionary:
        fp_to_int     - "always" or "if_none" or "never"
                         When to attempt cast of floating point to int
        fp_to_fp      - "always" or "if_none" or "never"
                         When to attempt cast of floating point to floating point
        int_to_int    - "always" or "if_none" or "never"
                         When to attempt cast of int to int
        int_to_fp     - "always" or "if_none" or "never"
                         When to attempt cast of int to floating point
        downcast_only - if True, only return datatype of same size or less
        downcast_within_fp - if True, tries downcasting within fp types, even
                             if there is an fp type that already matches
        guarantee_fp_to_fp_precision - if True, will only do fp to fp array
                        casting to type of same or higher precision. Note that
                        if fp_to_int recasting is allowed this will allow
                        precision loss of fp values
        prefer_input_at_threshold - number of bytes. If input array size
                        is less than or equal to this number, and in valid
                        types list, return the array without attempting
                        recasting
        prefer_int_type - if 'i', when recasting to integer type, prefer int
                        when equal sized uint is also available. Prefer
                        uint otherwise.
        '''
        if sctype_list is None:
            sctype_list = self._sctype_attributes.keys()
        self.sctype_list = sctype_list
        # Tolerances
        self.sctype_tols = self.default_sctype_tols()
        if sctype_tols is not None:
            self.sctype_tols.update(sctype_tols)
        # Casting options
        if recast_options is None:
            recast_options = 'only_if_none'
        if isinstance(recast_options, basestring):
            try:
                self.recast_options = self._option_defaults[recast_options]
            except KeyError:
                raise ValueError, \
                      'Did not recognize option string %s' % recast_options
        else:
            self.recast_options = self._option_defaults['only_if_none']
            self.recast_options.update(recast_options)
        # Cache sctype sizes,
        self.sized_sctypes = {}
        for k in ('c', 'f', 'i', 'u'):
            self.sized_sctypes[k] = self.sctypes_by_size(k)
        # Cache all integer sizes
        self.ints_sized_sctypes = []
        for k, v in self.sized_sctypes.items():
            if k in ('u', 'i'):
                for e in v:
                    self.ints_sized_sctypes.append(e)
        if self.ints_sized_sctypes:
            self.ints_sized_sctypes.sort(lambda x, y: cmp(y[1], x[1]))
        # Cache capable types list and sizes
        self._capable_sctypes = {}
        self._capable_sctype_sizes = {}
        self._c2f_capable_sctype_sizes = {}
        flts = self.sized_sctypes['f']
        for k in self._sctype_attributes:
            sct = self.get_capable_sctype(k)
            self._capable_sctypes[k] = sct
            if sct is None:
                self._capable_sctype_sizes[k] = inf
                if dtype(k).type == 'c':
                    self._c2f_capable_sctype_sizes[k] = inf
                continue
            dtp = dtype(sct)
            self._capable_sctype_sizes[k] = dtp.itemsize
            fsz = inf
            min_sz = ceil(dtp.itemsize / 2.0)
            if dtp.kind == 'c':
                for T, sz in flts:
                    if sz < min_sz:
                        break
                    fsz = sz
                self._c2f_capable_sctype_sizes[k] = fsz

    def default_sctype_tols(self):
        ''' Default allclose tolerance values for all dtypes '''
        t_dict = {}
        for sc_type in ('complex','float'):
            t_list = sctypes[sc_type]
            for T in t_list:
                dt = dtype(T)
                F = finfo(dt)
                t_dict[T] = {
                    'rtol': F.eps,
                    'atol': F.tiny}
        F = finfo(float64)
        for sc_type in ('int', 'uint'):
            t_list = sctypes[sc_type]
            for T in t_list:
                dt = dtype(T)
                t_dict[T] = {
                    'rtol': F.eps,
                    'atol': F.tiny}
        return t_dict

    def sctypes_by_size(self, kind):
        ''' Returns storage size ordered list of entries of scalar type sctype

        Input
        kind   - one of  "c",  "f", "i" or "u"
                 (for complex, float, integer, unsigned integer)
        '''
        D = []
        for t in self.sctype_list:
            dt = dtype(t)
            if dt.kind == kind:
                D.append([t, dt.itemsize])
        D.sort(lambda x, y: cmp(y[1], x[1]))
        return D

    def get_capable_sctype(self, sct):
        ''' Return smallest scalar type containing sct type without precision loss

        Input
        sct     - scalar type

        ID = input type. AT = acceptable type.  Return ID if ID is
        in ATs. Otherwise return smallest AT that is larger than or
        same size as ID.

        If the desired sctype is an integer, returns the smallest
        integer (int or uint) that can contain the range of the input
        integer type

        If there is no type that can contain sct without loss of
        precision, return None
        '''
        if sct in self.sctype_list:
            return sct
        out_t = None
        # Unsigned and signed integers
        # Precision loss defined by max min outside datatype range
        D = self._sctype_attributes[sct]
        if D['kind'] in ('u', 'i'):
            out_t = self.smallest_int_sctype(D['max'], D['min'])
        else:
            # Complex and float types
            # Precision loss defined by data size < sct
            sctypes = self.sized_sctypes[D['kind']]
            if not sctypes:
                return None
            dti = D['size']
            out_t = None
            for i, t in enumerate(sctypes):
                if t[1] >= dti:
                    out_t = t[0]
                else:
                    break
        return out_t

    def cast_to_fp(self, arr, kind,
                   max_size=inf,
                   continue_down=False):
        ''' Return fp arr maybe recast to specified kind, different sctype

        Inputs
        arr         - array to possibly recast
        kind        - kind of array to recast within
                      (one of "c", "f", "u", "i")
        max_size    - maximum size of sctype to return (in bytes)
        continue_down - if False, return array of largest sctype
                        within tolerance and >= max_size
                        if True, continue downcasting within kind
                        to find smallest possible within tolerance

        If arr cannot be recast within given tolerances, and size,
        return None
        '''
        tols = self.sctype_tols[arr.dtype.type]
        rtol, atol = tols['rtol'], tols['atol']
        ret_arr = None
        for T, sz in self.sized_sctypes[kind]:
            if sz > max_size:
                continue
            test_arr = arr.astype(T)
            if allclose(test_arr, arr, rtol, atol):
                ret_arr = test_arr
                if not continue_down:
                    break
            else:
                break
        return ret_arr

    def smallest_int_sctype(self, mx, mn, prefer='i'):
        ''' Return integer type with smallest storage containing mx and mn

        Inputs
        mx      - maximum value
        mn      - minumum value
        prefer  - if == 'i' prefer int for range also compatible
                  uint, else prefer uint in same situation

        Returns None if no integer can contain this range
        '''
        sct = None
        sz = inf
        for T, tsz in self.ints_sized_sctypes:
            t_dict = self._sctype_attributes[T]
            if t_dict['max'] >= mx and t_dict['min'] <= mn:
                if tsz < sz:
                    sct = T
                    sz = tsz
                elif tsz == sz:
                    if t_dict['kind'] == prefer:
                        sct = T
        return sct

    def cast_to_integer(self, arr, prefer='i'):
        ''' Casts arr to smallest integer containing range

        Returns None if range of arr cannot be contained in acceptable
        integer types

        prefer  - if == 'i' prefer int for range also compatible
                  uint, else prefer uint in same situation

        '''
        mx = amax(arr)
        mn = amin(arr)
        idt = self.smallest_int_sctype(mx, mn, prefer)
        if idt is not None:
            return arr.astype(idt)
        return None

    def recast(self, arr):
        ''' Recast array to type in type list

        If cannot recast to  an array within tolerance,
        raise error
        '''
        dtp = arr.dtype
        dtk = dtp.kind
        dti = dtp.itemsize
        dtt = dtp.type
        opts = self.recast_options
        curr_size = inf
        ret_arr = None
        valid_input_arr = dtt in self.sctype_list
        if valid_input_arr:
            if opts['prefer_input_at_threshold'] > arr.nbytes:
                return arr
            ret_arr = arr
        if opts['downcast_only'] or valid_input_arr:
            curr_size = dti
        tols = self.sctype_tols[dtt]
        rtol, atol = tols['rtol'], tols['atol']
        if dtk in ('c', 'f'):
            if opts['fp_to_int'] == 'always' or \
                   (opts['fp_to_int'] == 'if_none' and
                    ret_arr is None):
                test_arr = self.cast_to_integer(arr,
                                                opts['prefer_int_type'])
                if test_arr is not None and \
                   test_arr.dtype.itemsize < curr_size:
                    if allclose(arr, test_arr, rtol, atol):
                        ret_arr = test_arr
                        curr_size = ret_arr.dtype.itemsize
            if opts['fp_to_fp'] == 'always' or \
                   (opts['fp_to_fp'] == 'if_none' and
                    ret_arr is None):
                if dtk == 'c' and not opts['guarantee_fp_to_fp_precision']:
                    # Try casting to float
                    max_size = min([self._c2f_capable_sctype_sizes[dtt],
                                    curr_size - 1])
                    test_arr = self.cast_to_fp(arr,
                                               'f',
                                               max_size,
                                               opts['downcast_within_fp'])
                    if test_arr is not None:
                        ret_arr = test_arr
                        curr_size = ret_arr.dtype.itemsize
                if opts['fp_to_fp'] == 'always' or \
                       (opts['fp_to_fp'] == 'if_none' and
                        ret_arr is None):
                    # Cast float or complex to another of same type
                    if opts['guarantee_fp_to_fp_precision']:
                        sct = self._capable_sctypes[dtt]
                        sz = self._capable_sctype_sizes[dtt]
                        if sz < curr_size and sct is not None:
                            ret_arr = arr.astype(sct)
                            curr_size = sz
                    else:
                        max_size = min([self._capable_sctype_sizes[dtt],
                                        curr_size - 1])
                        test_arr = self.cast_to_fp(arr,
                                                   dtk,
                                                   max_size,
                                                   opts['downcast_within_fp'])
                        if test_arr is not None:
                            ret_arr = test_arr
                            curr_size = ret_arr.dtype.itemsize
        elif dtk in ('u', 'i'):
            if opts['int_to_int'] == 'always' or \
                   (opts['int_to_int'] == 'if_none' and
                    ret_arr is None):
                test_arr = self.cast_to_integer(arr,
                                                opts['prefer_int_type'])
                if test_arr is not None and \
                       test_arr.dtype.itemsize < curr_size:
                    ret_arr = test_arr
                    curr_size = ret_arr.dtype.itemsize
            if opts['int_to_fp'] == 'always' or \
                   (opts['int_to_fp'] == 'if_none' and
                    ret_arr is None):
                test_arr = self.cast_to_fp(arr,
                                           'f',
                                           curr_size-1,
                                           opts['downcast_within_fp'])
                if test_arr is not None:
                    ret_arr = test_arr
        else:
            raise TypeError, 'Do not recognize array kind %s' % dtk

        if ret_arr is not None:
            return ret_arr
        raise RecastError, 'Cannot recast array within tolerance'

    def recast_best_sctype(self, arr):
        ''' Recast array, return closest sctype to original

        Returns tuple of recast array and best sctype to contain
        original data before recasting
        '''
        sct = arr.dtype.type
        arr = self.recast(arr)
        if sct not in self.sctype_list:
            sct = self._capable_sctypes[sct]
            if sct is None:
                sct = arr.dtype.type
        return arr, sct
