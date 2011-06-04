import numpy as np
from numpy.testing import TestCase, assert_array_almost_equal, dec, \
                          assert_equal, assert_

from common import FUNCS_TP, FLAPACK_IS_EMPTY, CLAPACK_IS_EMPTY, FUNCS_FLAPACK, \
                   FUNCS_CLAPACK, PREC

SYEV_ARG = np.array([[1,2,3],[2,2,3],[3,3,6]])
SYEV_REF = np.array([-0.6699243371851365, 0.4876938861533345,
                     9.182230451031804])

class TestEsv(TestCase):
    def _test_base(self, func, lang):
        tp = FUNCS_TP[func]
        a = SYEV_ARG.astype(tp)
        if lang == 'C':
            f = FUNCS_CLAPACK[func]
        elif lang == 'F':
            f = FUNCS_FLAPACK[func]
        else:
            raise ValueError("Lang %s ??" % lang)

        w, v, info = f(a)

        assert_(not info, msg=repr(info))
        assert_array_almost_equal(w, SYEV_REF, decimal=PREC[tp])
        for i in range(3):
            assert_array_almost_equal(np.dot(a,v[:,i]), w[i]*v[:,i],
                                      decimal=PREC[tp])

    def _test_base_irange(self, func, irange, lang):
        tp = FUNCS_TP[func]
        a = SYEV_ARG.astype(tp)
        if lang == 'C':
            f = FUNCS_CLAPACK[func]
        elif lang == 'F':
            f = FUNCS_FLAPACK[func]
        else:
            raise ValueError("Lang %s ??" % lang)

        w, v, info = f(a, irange=irange)
        rslice = slice(irange[0], irange[1]+1)
        m = irange[1] - irange[0] + 1
        assert_(not info, msg=repr(info))

        assert_equal(len(w),m)
        assert_array_almost_equal(w, SYEV_REF[rslice], decimal=PREC[tp])

        for i in range(m):
            assert_array_almost_equal(np.dot(a,v[:,i]), w[i]*v[:,i],
                                      decimal=PREC[tp])

    def _test_base_vrange(self, func, vrange, lang):
        tp = FUNCS_TP[func]
        a = SYEV_ARG.astype(tp)
        ew = [value for value in SYEV_REF if vrange[0] < value <= vrange[1]]

        if lang == 'C':
            f = FUNCS_CLAPACK[func]
        elif lang == 'F':
            f = FUNCS_FLAPACK[func]
        else:
            raise ValueError("Lang %s ??" % lang)

        w, v, info = f(a, vrange=vrange)
        assert_(not info, msg=repr(info))

        assert_array_almost_equal(w, ew, decimal=PREC[tp])

        for i in range(len(w)):
            assert_array_almost_equal(np.dot(a,v[:,i]), w[i]*v[:,i],
                                      decimal=PREC[tp])

    def _test_syevr_ranges(self, func, lang):
        for irange in ([0, 2], [0, 1], [1, 1], [1, 2]):
            self._test_base_irange(func, irange, lang)

        for vrange in ([-1, 10], [-1, 1], [0, 1], [1, 10]):
            self._test_base_vrange(func, vrange, lang)

    # Flapack tests
    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_ssyev(self):
        self._test_base('ssyev', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dsyev(self):
        self._test_base('dsyev', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_ssyevr(self):
        self._test_base('ssyevr', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dsyevr(self):
        self._test_base('dsyevr', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_ssyevr_ranges(self):
        self._test_syevr_ranges('ssyevr', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dsyevr_ranges(self):
        self._test_syevr_ranges('dsyevr', 'F')

    # Clapack tests
    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["ssyev"],
                "Clapack empty, skip clapack test")
    def test_clapack_ssyev(self):
        self._test_base('ssyev', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dsyev"],
                "Clapack empty, skip clapack test")
    def test_clapack_dsyev(self):
        self._test_base('dsyev', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["ssyevr"],
                "Clapack empty, skip clapack test")
    def test_clapack_ssyevr(self):
        self._test_base('ssyevr', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dsyevr"],
                "Clapack empty, skip clapack test")
    def test_clapack_dsyevr(self):
        self._test_base('dsyevr', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["ssyevr"],
                "Clapack empty, skip clapack test")
    def test_clapack_ssyevr_ranges(self):
        self._test_syevr_ranges('ssyevr', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dsyevr"],
                "Clapack empty, skip clapack test")
    def test_clapack_dsyevr_ranges(self):
        self._test_syevr_ranges('dsyevr', 'C')

if __name__=="__main__":
    run_module_suite()
