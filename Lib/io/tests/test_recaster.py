from numpy.testing import *
import numpy as N

set_package_path()
from io.recaster import sctype_attributes, Recaster
restore_path()

try:  # Python 2.3 support
    from sets import Set as set
except:
    pass

class test_recaster(ScipyTestCase):
    def setUp(self):
        self.valid_types = [N.int32, N.complex128, N.float64]
        self.recaster = Recaster(self.valid_types,
                                 recast_options='smallest')
    
    def test_init(self):
        # Setting sctype_list
        R = Recaster()
        assert set(R.sctype_list) == set(sctype_attributes().keys()), \
                               'Default recaster should include all system types'
        T = N.float32
        R = Recaster([T])
        assert R.sctype_list == [T], 'Scalar type list not correctly set'
        # Setting tolerances
        tols = self.recaster.default_sctype_tols()
        assert tols == self.recaster.sctype_tols, 'Unexpected tols dictionary'
        F = N.finfo(T)
        R = Recaster(sctype_tols={T: {'rtol': F.eps*2, 'atol': F.tiny*2, 'silly': 'silly text'}})
        assert tols != R.sctype_tols, 'Tols dictionary not set correctly'
        assert R.sctype_tols[T]['rtol'] == F.eps*2, 'Rtol not correctly set'
        assert R.sctype_tols[T]['atol'] == F.tiny*2, 'Atol not correctly set'
        # Options
        # Sctype size lists
        # Integer sizes
        # Cabable types
        
    def test_cast_to_fp(self):
        R = self.recaster
        value = 1
        # Define expected type output from fp recast of value
        inp_outp = (
            (N.complex128, N.complex128),
            (N.complex64, N.complex128),
            )
        for inp, outp in inp_outp:
            arr = N.array(value, dtype=inp)
            rtol = R.sctype_tols[inp]['rtol']
            atol = R.sctype_tols[inp]['atol']
            kind = N.dtype(inp).kind
            arr = R.cast_to_fp(arr, rtol, atol, kind)
            if outp is None:
                assert arr is None, 'Expected None from type %s' % inp
            assert arr.dtype.type is outp, 'Expected output type %s from input %s' % (inp, outp)
            
    def test_smallest_int_sctype(self):
        # Smallest int sctype with testing recaster
        params = sctype_attributes()
        mmax = params[N.int32]['max']
        mmin = params[N.int32]['min']        
        for kind in ('int', 'uint'):
            for T in N.sctypes[kind]:
                mx = params[T]['max']
                mn = params[T]['min']
                rt = self.recaster.smallest_int_sctype(mx, mn)
                if mx <= mmax and mn >= mmin:
                    assert rt == N.int32, 'Expected int32 type'
                else:
                    assert rt is None, 'Expected None, got %s for %s' % (T, rt)
                   
        # Smallest int sctype with full recaster
        RF = Recaster()
        test_triples = [(N.uint8, 0, 255),
                      (N.int8, -128, 0),
                      (N.uint16, 0, params[N.uint16]['max']),                      
                      (N.int16, params[N.int16]['min'], 0),
                      (N.uint32, 0, params[N.uint32]['max']),
                      (N.int32, params[N.int32]['min'], 0),
                      (N.uint64, 0, params[N.uint64]['max']),
                      (N.int64, params[N.int64]['min'], 0)]
        for T, mn, mx in test_triples:
            rt = RF.smallest_int_sctype(mx, mn)
            assert N.dtype(rt) == N.dtype(T), \
                   'Expected %s, got %s type' % (T, rt)
        
    def test_recasts(self):
        value = 100
        R = self.recaster
        for T in (N.complex128, N.complex64,
                  N.float64, N.uint64):
            B = R.recast(N.array(value, T))
            assert B is not None, 'Got None for %s' % T
            Bt = B.dtype.type
            assert Bt == N.int32, 'Input %s, output %s' % (T, Bt)
        
