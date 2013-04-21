#!/usr/bin/env python
#
# Created by: Pearu Peterson, September 2002
#
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import *

from common import FUNCS_TP, FUNCS_CLAPACK, FUNCS_FLAPACK, FLAPACK_IS_EMPTY, \
                   CLAPACK_IS_EMPTY


class TestLapack(TestCase):
    def _test_gebal_base(self, func, lang):
        tp = FUNCS_TP[func]

        a = np.array([[1,2,3],[4,5,6],[7,8,9]]).astype(tp)
        a1 = np.array([[1,0,0,3e-4],
                       [4,0,0,2e-3],
                       [7,1,0,0],
                       [0,1,0,0]]).astype(tp)

        if lang == 'C':
            f = FUNCS_CLAPACK[func]
        elif lang == 'F':
            f = FUNCS_FLAPACK[func]
        else:
            raise ValueError("Lang %s ??" % lang)

        ba, lo, hi, pivscale, info = f(a)
        assert_(not info, msg=repr(info))
        assert_array_almost_equal(ba, a)
        assert_equal((lo,hi), (0, len(a[0])-1))
        assert_array_almost_equal(pivscale, np.ones(len(a)))

        ba, lo, hi, pivscale, info = f(a1,permute=1,scale=1)
        assert_(not info, msg=repr(info))

    def _test_gehrd_base(self, func, lang):
        tp = FUNCS_TP[func]

        a = np.array([[-149, -50,-154],
             [537, 180, 546],
             [-27, -9, -25]]).astype(tp)

        if lang == 'C':
            f = FUNCS_CLAPACK[func]
        elif lang == 'F':
            f = FUNCS_FLAPACK[func]
        else:
            raise ValueError("Lang %s ??" % lang)

        ht, tau, info = f(a)
        assert_(not info, msg=repr(info))

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_sgebal(self):
        self._test_gebal_base('sgebal', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip flapack test")
    def test_dgebal(self):
        self._test_gebal_base('dgebal', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip clapack test")
    def test_sgehrd(self):
        self._test_gehrd_base('sgehrd', 'F')

    @dec.skipif(FLAPACK_IS_EMPTY, "Flapack empty, skip clapack test")
    def test_dgehrd(self):
        self._test_gehrd_base('dgehrd', 'F')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["sgebal"],
                "Clapack empty, skip flapack test")
    def test_clapack_sgebal(self):
        self._test_gebal_base('sgebal', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dgebal"],
                "Clapack empty, skip flapack test")
    def test_clapack_dgebal(self):
        self._test_gebal_base('dgebal', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["sgehrd"],
                "Clapack empty, skip flapack test")
    def test_clapack_sgehrd(self):
        self._test_gehrd_base('sgehrd', 'C')

    @dec.skipif(CLAPACK_IS_EMPTY or not FUNCS_CLAPACK["dgehrd"],
                "Clapack empty, skip flapack test")
    def test_clapack_dgehrd(self):
        self._test_gehrd_base('dgehrd', 'C')
