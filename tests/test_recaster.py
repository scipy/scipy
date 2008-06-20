# Test recasting module

from numpy.testing import *
import numpy as N

set_package_path()
from io.recaster import basecaster, intcaster, floatcaster, complexcaster, \
     rejectcastercollection, precisioncastercollection,\
     energeticcastercollection, fastcastercollection, \
     numerictypeinfo, RecastError
restore_path()

try:  # Python 2.3 support
    from sets import Set as set
except:
    pass


class test_numerictypeinfo(NumpyTestCase):
    def test_init(self):
        types = []
        for k in ('int', 'float', 'uint', 'complex'):
            types += N.sctypes[k]
        # All system numeric types by default
        nt = numerictypeinfo()
        assert nt.types == set(types)
        # Can include bools as uint type
        nt = numerictypeinfo(bool_as_uint=True)
        types.append(N.bool)
        assert nt.types == set(types)
        # Bools rejected unless flag is set
        self.assertRaises(TypeError, numerictypes,[N.bool])
        # Can restrict by kind
        nt = numerictypeinfo(kinds=['int', 'uint'])
        types = []
        for k in ('int', 'uint'):
            types += N.sctypes[k]
        assert nt.types == set(types)
        # By list of types
        type_list = [N.uint16, N.int16, N.float64]
        nt = numertypeinfo(types=type_list)
        assert nt.types == set(type_list)
        # And by intersection of two
        nt = numerictypeinfo(kinds=['int', 'uint'], types=type_list)
        assert nt.types == set([N.uint16, N.int16])
        # Reject non-numeric
        self.assertRaises(TypeError, numerictypeinfo, [N.void])

    def test_info(self):
        nt = numerictypeinfo()
        it_expected = {'kind': 'i',
                       'size': 1,
                       'min': -128,
                       'max': 127}
        assert nt.info(N.int8) == it_expected
        F = N.finfo(N.dtype(N.float64))
        ft_expected = {'kind': 'f',
                       'size', 8,
                       'min': F.min,
                       'max': F.max}
        assert nt.info(N.float64) == ft_expected

    def test_of_kind(self):
        # Can select sublist of specified kind
        nt = numerictypes()
        allints = set(N.sctypes['uint'] + N.sctypes['uint'] + [N.bool])
        assert set(nt.of_kind(['int', 'uint')) == uints

    def test_by_size(self):
        nt = numerictypes()
        t_by_size = nt.by_size()
        csz = N.inf
        for t in t_by_size:
            sz = N.dtype(t).itemsize
            assert sz <= csz
            csz = sz
        nt = numerictypes([N.int8, N.int16])            
        assert nt.by_size() == [N.int16, N.int8]
        # uints are lower (appear smaller) in search order
        # because larger range for same size means smaller
        nt = numerictypes([N.uint16, N.int16])
        assert nt.by_size() == [N.int16, N.uint16]
        # bools therefore tend to be higher in order
        nt = numerictypes([N.uint8, N.int8, N.bool], bool_as_uint=True)
        assert nt.by_size() == [N.bool, N.int8, N.uint8]

    def test_smallest_precise(self):
        nt = numerictypes(kinds=['complex'])
        largest_complex = nt.by_size()[0]
        nt = numerictypes(kinds=['float'])
        largest_float = nt.by_size()[0]        
        nt = numerictypes(types=[N.complex128, N.float64, N.int16])
        tests = [[N.int8, N.int16],
                 [N.bool, N.int16],
                 [N.uint8, N.int16],
                 [N.int16, N.int16],
                 [N.uint16, None],
                 [N.int32, None],
                 [N.float32, N.float64],
                 [N.float64, N.float64],
                 [N.float128, None],
                 [N.complex64, N.complex128],
                 [N.complex128, N.complex128],
                 [largest_complex, None]]
        for inp, outp in tests:
            assert nt.smallest_precise(inp) == outp
        # No floats, then float goes to complex
        nt = numerictypes(types=[N.complex128])
        assert nt.smallest_precise(N.float64) == N.complex128
        assert nt.smallest_precise(N.float32) == N.complex128
        assert nt.smallest_precise(largest_float) == None


class CastTestCase(NumpyTestCase):
    ''' Define helper function for running caster tests '''
    def run_rig(self, caster, test_list):
        for inpval, inptype, outval, outtype in test_list:
            inparr = N.array(inpval, dtype=inptype)
            if outval is None:
                self.assertRaises(RecastError, caster.do, inparr)
                continue
            res = caster.do(inparr)
            assert res.dtype == N.dtype(outtype)
            assert res == outval
    

class test_basecaster(CastTestCase):
    def test_init(self):
        bc = basecaster()
        assert bc.ordered_types is None
        bc = basecaster([N.float64])
        assert bc.ordered_types == [N.float64]        
        bc = basecaster([N.float64], tols={})
        bc = basecaster([N.int8, N.float64])
        assert bc.ordered_types = [N.float64, N.int8]
        # Reject non-numeric
        self.assertRaises(TypeError, intcaster, [N.void])
       
    def test_do(self):
        bc = basecaster()
        self.assertRaises(NotImplementedError, bc.do, 1)


class test_intcaster(CastTestCase):
    def test_init(self):
        # Default leads to None in typeinfo
        ic = intcaster()
        assert ic.typeinfo is None
        # Ordering
        ic = intcaster([N.int8, N.int16])
        assert bc.typeinfo.types == set([N.int16, N.int8])
        assert ic.tols is None
        # Test passed named args
        F = N.finfo(N.dtype(N.float64))
        tols = {N.int8: {'rtol': F.eps, 'atol': F.tiny}}
        ic = intcaster([N.int8], tols=tols)
        assert ic.tols == tols
        # Reject non-integer
        self.assertRaises(TypeError, intcaster, [N.float64])
        self.assertRaises(TypeError, intcaster, [N.complex128])
        # Accept bool
        ic = intcaster([N.bool])
        
    def test_int_int(self):
        ic = intcaster()
        self.assertRaises(RecastError, ic.do, 1)
        # Default caster tries to find smallest type with larger type range
        # If not, falls back to largest type containing data range
        ic = intcaster([N.int8])
        # Tests: input value, input dtype, output value, output dtype
        # Where None in output value indicates should raise error
        tests = [ # Inevitable shrink
            [1, N.int32, 1, N.int8],
            [128, N.int32, None, None],
            [-129, N.int32, None, None],
            ]
        self.run_rig(ic, tests)
        ic = intcaster([N.int32, N.int8])
        tests = [ # Up by default, down when 

    def test_float_int(self):


class test_floatcaster(CastTestCase):
    def test_init(self):
        pass

    def test_float_float(self):
        pass

    def test_int_float(self):
        pass

    def test_complex_float(self):
        pass


class test_complexcaster(CastTestCase):
    def test_init(self):
        pass

    def test_complex_complex(self):
        pass

    def test_complex_float(self):
        pass

    def test_complex_int(self):
        ''' Necessary?  Maybe only do complex->float->int '''
        pass
    
class test_rejectecastercollection(NumpyTestCase):
    def test_init(self):
        pass

    def test_do(self):
        pass
    
class test_precisioncastercollection(NumpyTestCase):
    def test_init(self):
        pass

    def test_do(self):
        pass
    
class test_energeticcastercollection(NumpyTestCase):
    def test_init(self):
        pass

    def test_do(self):
        pass
    
class test_fastcastercollection(NumpyTestCase):
    def test_init(self):
        pass

    def test_do(self):
        pass


class test_smallcastercollection(NumpyTestCase):
    def test_init(self):
        pass

    def test_do(self):
        pass
