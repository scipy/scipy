import numpy as N
from numpy.testing import *

from scipy.io.recaster import sctype_attributes, Recaster, RecastError

try:  # Python 2.3 support
    from sets import Set as set
except:
    pass

class TestRecaster(TestCase):

    def test_init(self):
        # Setting sctype_list
        R = Recaster()
        assert set(R.sctype_list) == set(sctype_attributes().keys()), \
                               'Default recaster should include all system types'
        T = N.float32
        R = Recaster([T])
        assert R.sctype_list == [T], 'Scalar type list not correctly set'
        # Setting tolerances
        R = Recaster()
        tols = R.default_sctype_tols()
        assert tols == R.sctype_tols, 'Unexpected tols dictionary'
        F = N.finfo(T)
        R = Recaster(sctype_tols={T: {
            'rtol': F.eps*2,
            'atol': F.tiny*2,
            'silly': 'silly text'}})
        assert R.sctype_tols[T]['rtol'] == F.eps*2, \
               'Rtol not correctly set'
        assert R.sctype_tols[T]['atol'] == F.tiny*2, \
               'Atol not correctly set'
        T = N.complex128
        F = N.finfo(T)
        assert R.sctype_tols[T]['rtol'] == F.eps, \
               'Rtol defaults not correctly set'
        assert R.sctype_tols[T]['atol'] == F.tiny, \
               'Atol defaults not correctly set'
        # Options
        # Sctype size lists
        # Integer sizes
        # Cabable types

    def test_cast_to_fp(self):
        R = Recaster()
        # Define expected type output from fp recast of value
        sta = sctype_attributes()
        inp_outp = (
            (1, N.complex128, 'c', sta[N.complex128]['size'], 0, N.complex128),
            (1, N.complex128, 'c', sta[N.complex128]['size'], 1, N.complex64),
            (1, N.complex128, 'c', sta[N.complex64]['size'], 0, N.complex64),
            (1, N.complex128, 'f', sta[N.float64]['size'], 0, N.float64),
            (1.0+1j, N.complex128, 'f', sta[N.complex128]['size'], 0, None),
            (1, N.float64, 'f', sta[N.float64]['size'], 0, N.float64),
            (1, N.float64, 'f', sta[N.float64]['size'], 1, N.float32),
            (1, N.float64, 'f', sta[N.float32]['size'], 0, N.float32),
            (1, N.float64, 'c', sta[N.complex128]['size'], 0, N.complex128),
            (1, N.float64, 'c', sta[N.complex128]['size'], 1, N.complex64),
            (1, N.int32, 'f', sta[N.float64]['size'], 0, N.float64),
            (1, N.int32, 'f', sta[N.float64]['size'], 1, N.float32),
            (1, N.float64, 'f', 0, 0, None),
            )
        for value, inp, kind, max_size, continue_down, outp in inp_outp:
            arr = N.array(value, dtype=inp)
            arr = R.cast_to_fp(arr, kind, max_size, continue_down)
            if outp is None:
                assert arr is None, \
                       'Expected None from type %s, got %s' \
                       % (inp, arr.dtype.type)
                continue
            assert arr is not None, \
                   'Expected %s from %s, got None' % (outp, inp)
            dtt = arr.dtype.type
            assert dtt is outp, \
                   'Expected %s from %s, got %s' % (outp, inp, dtt)

    def test_smallest_int_sctype(self):
        # Smallest int sctype with full recaster
        params = sctype_attributes()
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
        # Smallest int sctype with restricted recaster
        mmax = params[N.int32]['max']
        mmin = params[N.int32]['min']
        RR = Recaster([N.int32])
        for kind in ('int', 'uint'):
            for T in N.sctypes[kind]:
                mx = params[T]['max']
                mn = params[T]['min']
                rt = RR.smallest_int_sctype(mx, mn)
                if mx <= mmax and mn >= mmin:
                    assert rt == N.int32, \
                           'Expected int32 type, got %s' % rt
                else:
                    assert rt is None, \
                           'Expected None, got %s for %s' % (T, rt)
        # Test preferred int flag
        mx = 1000
        mn = 0
        rt = RF.smallest_int_sctype(mx, mn)
        assert rt == N.int16, 'Expected int16, got %s' % rt
        rt = RF.smallest_int_sctype(mx, mn, 'i')
        assert rt == N.int16, 'Expected int16, got %s' % rt
        rt = RF.smallest_int_sctype(mx, mn, prefer='u')
        assert rt == N.uint16, 'Expected uint16, got %s' % rt

    def test_recasts(self):
        valid_types = [N.int32, N.complex128, N.float64]
        # Test smallest
        R = Recaster(valid_types, recast_options='smallest')
        inp_outp = (
            (1, N.complex128, N.int32),
            (1, N.complex64, N.int32),
            (1.0+1j, N.complex128, N.complex128),
            (1.0+1j, N.complex64, N.complex128),
            (1, N.float64, N.int32),
            (1, N.float32, N.int32),
            (1.1, N.float64, N.float64),
            (-1e12, N.int64, N.float64),
            )
        self.run_io_recasts(R, inp_outp)
        # Test only_if_none
        R = Recaster(valid_types, recast_options='only_if_none')
        inp_outp = (
            (1, N.complex128, N.complex128),
            (1, N.complex64, N.int32),
            (1.0+1j, N.complex128, N.complex128),
            (1.0+1j, N.complex64, N.complex128),
            (1, N.float64, N.float64),
            (1, N.float32, N.int32),
            (1.1, N.float64, N.float64),
            (-1e12, N.int64, N.float64),
            )
        self.run_io_recasts(R, inp_outp)
        # Test preserve_precision
        R = Recaster(valid_types, recast_options='preserve_precision')
        inp_outp = (
            (1, N.complex128, N.complex128),
            (1, N.complex64, N.complex128),
            (1.0+1j, N.complex128, N.complex128),
            (1.0+1j, N.complex64, N.complex128),
            (1, N.float64, N.float64),
            (1, N.float32, N.float64),
            (1.1, N.float64, N.float64),
            (-1e12, N.int64, None),
            )
        self.run_io_recasts(R, inp_outp)

    def run_io_recasts(self, R, inp_outp):
        ''' Runs sets of value, input, output tests '''
        for value, inp, outp in inp_outp:
            arr = N.array(value, inp)
            if outp is None:
                self.assertRaises(RecastError, R.recast, arr)
                continue
            arr = R.recast(N.array(value, inp))
            assert arr is not None, \
                   'Expected %s from %s, got None' % (outp, inp)
            dtt = arr.dtype.type
            assert dtt is outp, \
                   'Expected %s from %s, got %s' % (outp, inp, dtt)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
