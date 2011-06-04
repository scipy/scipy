import numpy as np
from numpy.testing import TestCase, assert_array_almost_equal, dec, \
                          assert_equal, assert_

from common import FUNCS_TP, FLAPACK_IS_EMPTY, CLAPACK_IS_EMPTY, FUNCS_FLAPACK, \
                   FUNCS_CLAPACK, PREC

A = np.array([[1,2,3],[2,2,3],[3,3,6]])
B = np.array([[10,-1,1],[-1,8,-2],[1,-2,6]])

class TestSygv(TestCase):
    def _test_base(self, func, lang, itype):
        tp = FUNCS_TP[func]
        a = A.astype(tp)
        b = B.astype(tp)
        if lang == 'C':
            f = FUNCS_CLAPACK[func]
        elif lang == 'F':
            f = FUNCS_FLAPACK[func]
        else:
            raise ValueError("Lang %s ??" % lang)

        w, v, info = f(a, b, itype=itype)

        assert_(not info, msg=repr(info))
        for i in range(3):
            if itype == 1:
                assert_array_almost_equal(np.dot(a,v[:,i]), w[i]*np.dot(b,v[:,i]),
                                          decimal=PREC[tp])
            elif itype == 2:
                assert_array_almost_equal(np.dot(a,np.dot(b,v[:,i])), w[i]*v[:,i],
                                          decimal=PREC[tp])
            elif itype == 3:
                assert_array_almost_equal(np.dot(b,np.dot(a,v[:,i])),
                                          w[i]*v[:,i], decimal=PREC[tp] - 1)
            else:
                raise ValueError(itype)

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_ssygv_1(self):
        self._test_base('ssygv', 'F', 1)

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_ssygv_2(self):
        self._test_base('ssygv', 'F', 2)

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_ssygv_3(self):
        self._test_base('ssygv', 'F', 3)

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dsygv_1(self):
        self._test_base('dsygv', 'F', 1)

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dsygv_2(self):
        self._test_base('dsygv', 'F', 2)

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dsygv_3(self):
        self._test_base('dsygv', 'F', 3)

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["ssygv"],
                "Clapack empty, skip flapack test")
    def test_clapack_ssygv_1(self):
        self._test_base('ssygv', 'C', 1)

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["ssygv"],
                "Clapack empty, skip flapack test")
    def test_clapack_ssygv_2(self):
        self._test_base('ssygv', 'C', 2)

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["ssygv"],
                "Clapack empty, skip flapack test")
    def test_clapack_ssygv_3(self):
        self._test_base('ssygv', 'C', 3)

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dsygv"],
                "Clapack empty, skip flapack test")
    def test_clapack_dsygv_1(self):
        self._test_base('dsygv', 'C', 1)

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dsygv"],
                "Clapack empty, skip flapack test")
    def test_clapack_dsygv_2(self):
        self._test_base('dsygv', 'C', 2)

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dsygv"],
                "Clapack empty, skip flapack test")
    def test_clapack_dsygv_3(self):
        self._test_base('dsygv', 'C', 3)

if __name__=="__main__":
    run_module_suite()
