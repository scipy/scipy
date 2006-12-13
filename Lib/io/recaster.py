# Author: Matthew Brett

"""
Recaster class for recasting numeric arrays
"""

from numpy import *

def sctype_attributes():
    ''' Return dictionary describing numpy scalar types '''
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

class Recaster(object):
    ''' Class to recast arrays to one of acceptable scalar types

    Initialization specifies acceptable types (ATs)

    Implements downcast and recast method - returns array that may be
    of different storage type to the input array, where the new type
    is one of the ATs. Downcast forces return array to be same size or
    smaller than the input.  recast method will return a larger type
    if no smaller type will contain the data without loss of
    precision.

    At its simplest, the downcast method can reject arrays that
    are not in the list of ATs.
    '''

    _sctype_attributes = sctype_attributes()

    def __init__(self, sctype_list=None,
                 downcast_fp_to_fp = True,
                 downcast_fp_to_int = True,
                 downcast_int_to_int = True,
                 upcast_int_to_fp = True,
                 upcast_fp_to_int = True,
                 sctype_tols=None):
        ''' Set types for which we are attempting to downcast

        Input
        sctype_list  - list of acceptable scalar types
                     If None defaults to all system types
        downcast_fp_to_fp - if True, tries to downcast floats and complex
                       to smaller size of same type
        downcast_fp_to_int - if True, tries to downcast floats and complex
                       to integers
        downcast_int_to_int - if True, tries to downcast integers to
                       smaller of same type
        upcast_int_to_fp - if True, tries to upcast integers that could not
                       be downcast to floating point type
        upcast_fp_to_int - if True, tries to upcast floating point arrays
                       that cannot be downcast, to integers
        sctype_tols  - dictionary key datatype, values rtol, tol
                     to specify tolerances for checking near equality in
                     downcasting

        Note that tolerance values for integers are used for upcasting
        integers to floats
        '''
        if sctype_list is None:
            sctype_list = self._sctype_attributes.keys()
        self.sctype_list = sctype_list
        # Casting options
        self.sctype_tols = self.default_sctype_tols()
        self.downcast_fp_to_fp = downcast_fp_to_fp
        self.downcast_fp_to_int = downcast_fp_to_int
        self.downcast_int_to_int = downcast_int_to_int
        self.upcast_int_to_fp = upcast_int_to_fp
        self.upcast_fp_to_int = upcast_fp_to_int
        # Tolerances
        if sctype_tols is not None:
            self.sctype_tols.update(sctype_tols)
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
        # Cache capable types list
        self._capable_sctypes = {}
        for k in self._sctype_attributes:
            self._capable_sctypes[k] = self.get_capable_sctype(k)
    
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

    def capable_sctype(self, sct):
        ''' Return smallest type containing sct type without precision loss

        Value pulled fron dictionary cached from init - see
        get_capable_sctype method for algorithm
        '''
        try:
            return self._capable_sctypes[sct]
        except KeyError:
            return None
        
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

    def all_close(self, arr1, arr2):
        ''' True if arr1 arr2 close with tols for arr1 '''
        tols = self.sctype_tols[arr1.dtype.type]
        return allclose(arr1, arr2,
                        rtol=tols['rtol'],
                        atol=tols['atol'])

    def arr_if_valid(self, arr):
        ''' Returns array if of valid sctype, None otherwise '''
        if arr.dtype.type not in self.sctype_list:
            return None
        return arr

    def smallest_of_kind(self, arr, kind=None, max_size=None):
        ''' Return arr maybe downcast to same kind, smaller storage

        Inputs
        arr         - array to possibly downcast
        kind        - kind of array to downcast within
                      (if None (default) use arr.dtype.kind)
        max_size    - maximum size of sctype to return (in bytes)
                      (if None, set to arr.dtype.itemsize-1)
        If arr cannot be downcast within given tolerances, then:
        return arr if arr is in list of acceptable types, otherwise
        return None
        '''
        dtp = arr.dtype
        if kind is None:
            kind = dtp.kind
        if max_size is None:
            max_size = dtp.itemsize-1
        sctypes = self.sized_sctypes[kind]
        sctypes = [t[0] for i, t in enumerate(sctypes) if t[1] <= max_size]
        tols = self.sctype_tols[dtp.type]
        rtol, atol = tols['rtol'], tols['atol']
        ret_arr = arr
        for T in sctypes:
            test_arr = arr.astype(T)
            if allclose(test_arr, arr, rtol, atol):
                ret_arr = test_arr
            else:
                break
        return self.arr_if_valid(ret_arr)
        
    def smallest_int_sctype(self, mx, mn):
        ''' Return integer type with smallest storage containing mx and mn

        Inputs
        mx      - maximum value
        mn      - minumum value

        Returns None if no integer can contain this range
        '''
        sct = None
        for T, tsz in self.ints_sized_sctypes:
            t_dict = self._sctype_attributes[T]
            if t_dict['max'] >= mx and t_dict['min'] <= mn:
                if sct is None or tsz < sz:
                    sct = T
                    sz = tsz
        return sct

    def cast_to_integer(self, arr):
        ''' Casts arr to smallest integer containing range

        Returns None if range of arr cannot be contained in acceptable
        integer types
        '''
        mx = amax(arr)
        mn = amin(arr)
        idt = self.smallest_int_sctype(mx, mn)
        if idt is not None:
            return arr.astype(idt)
        return None

    def downcast_or_none(self, arr):
        ''' Downcast array to smaller or same type
        
        If cannot find smaller type within tolerance,
        return array if is already valid type, otherwise None
        '''
        dtp = arr.dtype
        dtk = dtp.kind
        dti = dtp.itemsize
        if dtk in ('c', 'f'):
            if self.downcast_fp_to_int:
                test_arr = self.cast_to_integer(arr)
                if test_arr is not None:
                    if self.all_close(arr, test_arr):
                        return test_arr
            if self.downcast_fp_to_fp:
                if dtk == 'c':
                    # Try downcasting to float
                    max_size = ceil(dti / 2.0)
                    test_arr = self.smallest_of_kind(arr, 'f', max_size)
                    if test_arr is not None:
                        return test_arr
                test_arr = self.smallest_of_kind(arr)
                if test_arr is not None:
                    return test_arr
        elif dtk in ('u', 'i'):
            if self.downcast_int_to_int:
                test_arr = self.cast_to_integer(arr)
                if test_arr is not None:
                    if test_arr.dtype.itemsize <= dti:
                        return test_arr
        else:
            raise TypeError, 'Do not recognize array kind %s' % dtk
        return self.arr_if_valid(arr)
        

    def recast_or_none(self, arr):
        ''' Recast array to type in type list
        
        If cannot find smaller type within tolerance, by downcasting,
        and array not of valid type already, then try larger
        types.  If none of these return an array within tolerance,
        return None
        '''
        test_arr = self.downcast_or_none(arr)
        if test_arr is not None:
            return test_arr
        # Could not downcast, arr dtype not in known list
        dtp = arr.dtype
        dtk = dtp.kind
        sct = dtp.type
        if dtk in ('c', 'f'):
            # Try upcast to larger dtype of same kind
            udt = self.capable_sctype[sct]
            if udt is not None:
                return arr.astype(udt)
            # Try casting to an integer
            if self.upcast_fp_to_int:
                test_arr = self.cast_to_integer(arr)
                if test_arr is not None:
                    if self.all_close(arr, test_arr):
                        return test_arr
        else: # integer types
            # try casting to any possible integer type
            test_arr = self.cast_to_integer(arr)
            if test_arr is not None:
                return test_arr
            # Can try casting integers to floats
            if self.upcast_int_to_fp:
                flts = self._sized_sctypes['f']
                if flts:
                    flt_arr = arr.astype(flts[0])
                    if self.all_close(arr, flt_arr):
                        if self.downcast_fp_to_fp:
                            max_size = flt_arr.dtype.itemsize - 1
                            test_arr = self.smallest_of_kind(arr, 'f', max_size)
                            if test_arr is not None:
                                return test_arr
                        return flt_arr
        return None

    def downcast(self, arr):
        ret = self.downcast_or_none(arr)
        if ret is None:
            raise TypeError, 'Cannot downcast array within tolerance'
        return ret
        
    def recast(self, arr):
        ret = self.recast_or_none(arr)
        if ret is None:
            raise TypeError, 'Cannot recast array within tolerance'
        return ret

    def recast_best_sctype(self, arr):
        ''' Recast array, return closest sctype to original

        Returns tuple of recast array and best sctype to contain
        original data before recasting
        '''
        sct = arr.dtype.type
        arr = self.recast(arr)
        if sct not in self.sctype_list:
            sct = self.capable_sctype[sct]
            if sct is None:
                sct = arr.dtype.type
        return arr, sct
