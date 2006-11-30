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
        self.recaster = Recaster([N.int32, N.complex64, N.float32])
    
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
        r, a = R.tols_from_sctype(T)
        assert r == F.eps*2, 'Rtol not correctly set'
        assert a == F.tiny*2, 'Atol not correctly set'
        # Sctype size lists
        # Integer sizes
        # Cabable types
        
    def test_methods(self):
        A = N.array(1, N.float64)
        B = A.astype(N.float32)
        # smallest from sctypes
        C = self.recaster.smallest_from_sctypes(A, [N.float32])
        # smaller same kind
        C = self.recaster.smallest_same_kind(A)
        assert C.dtype == N.dtype(N.float32), 'Dtype was not downcast'
