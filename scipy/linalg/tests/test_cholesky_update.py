import itertools

import numpy as np
from numpy.testing import assert_allclose
from pytest import raises as assert_raises
import scipy.linalg
from scipy.linalg._cholesky_update import cholesky_update


def assert_upper_tri(a, rtol=None, atol=None):
    if rtol is None:
        rtol = 10.0 ** -(np.finfo(a.dtype).precision-2)
    if atol is None:
        atol = 2*np.finfo(a.dtype).eps
    mask = np.tri(a.shape[0], a.shape[1], -1, np.bool_)
    assert_allclose(a[mask], 0.0, rtol=rtol, atol=atol)


def assert_lower_tri(a, rtol=None, atol=None):
    if rtol is None:
        rtol = 10.0 ** -(np.finfo(a.dtype).precision-2)
    if atol is None:
        atol = 2*np.finfo(a.dtype).eps
    mask = np.triu(np.ones((a.shape[0], a.shape[1]), dtype=np.bool_), 1)
    assert_allclose(a[mask], 0.0, rtol=rtol, atol=atol)


def assert_allclose_tri(exp, act, lower, rtol, atol):
    if lower:
        assert_allclose(np.tril(exp, -1), np.tril(act, -1),
                        rtol=rtol, atol=atol)
    else:
        assert_allclose(np.triu(exp, 1), np.triu(act, 1),
                        rtol=rtol, atol=atol)


def check_cholesky(a, r, lower, rtol, atol):
    if lower:
        assert_lower_tri(r, rtol=rtol, atol=atol)
        a1 = np.dot(r, r.T.conj())
    else:
        assert_upper_tri(r, rtol=rtol, atol=atol)
        a1 = np.dot(r.T.conj(), r)
    assert_allclose(a1, a, rtol=rtol, atol=atol)
    assert_allclose(scipy.linalg.cholesky(a, lower=lower), r, rtol=rtol, atol=atol)


def make_strided(arrs):
    strides = [(3, 7), (2, 2), (3, 4), (4, 2), (5, 4), (2, 3), (2, 1), (4, 5)]
    kmax = len(strides)
    k = 0
    ret = []
    for a in arrs:
        if a.ndim == 1:
            s = strides[k % kmax]
            k += 1
            base = np.zeros(s[0]*a.shape[0]+s[1], a.dtype)
            view = base[s[1]::s[0]]
            view[...] = a
        elif a.ndim == 2:
            s = strides[k % kmax]
            t = strides[(k+1) % kmax]
            k += 2
            base = np.zeros((s[0]*a.shape[0]+s[1], t[0]*a.shape[1]+t[1]),
                            a.dtype)
            view = base[s[1]::s[0], t[1]::t[0]]
            view[...] = a
        else:
            raise ValueError('make_strided only works for ndim = 1 or'
                             ' 2 arrays')
        ret.append(view)
    return ret


def negate_strides(arrs):
    ret = []
    for a in arrs:
        b = np.zeros_like(a)
        if b.ndim == 2:
            b = b[::-1, ::-1]
        elif b.ndim == 1:
            b = b[::-1]
        else:
            raise ValueError('negate_strides only works for ndim = 1 or'
                             ' 2 arrays')
        b[...] = a
        ret.append(b)
    return ret


def nonitemsize_strides(arrs):
    out = []
    for a in arrs:
        a_dtype = a.dtype
        b = np.zeros(a.shape, [('a', a_dtype), ('junk', 'S1')])
        c = b.getfield(a_dtype)
        c[...] = a
        out.append(c)
    return out


def make_nonnative(arrs):
    return [a.astype(a.dtype.newbyteorder()) for a in arrs]


class BaseCholdeltas:
    def setup_method(self):
        self.rtol = 10.0 ** -(np.finfo(self.dtype).precision)
        self.atol = 1 * np.finfo(self.dtype).eps

    def generate(self, type, which='1d'):
        np.random.seed(12)
        shape = {'sqr': (3, 3), '1x1': (1, 1)}[type]
        a = np.random.random(shape)
        if np.iscomplexobj(self.dtype.type(1)):
            b = np.random.random(shape)
            a = a + 1j * b
        a = a.astype(self.dtype)

        # Make symmetric and PD
        a = a + a.T.conj() + a.shape[0] * np.eye(a.shape[0], dtype=self.dtype)
        l = scipy.linalg.cholesky(a, lower=True)
        r = scipy.linalg.cholesky(a, lower=False)

        # Generate the update vector
        if which == '1d':
            u = np.random.random(a.shape[0])
        elif which == 'col':
            u = np.random.random((a.shape[0], 1))
        elif which == 'p':
            # u = np.random.random((2, a.shape[0]))
            u = np.random.random((a.shape[0], 2))
        else:
            raise ValueError('which should be either "1d" or "col" or "p"')
        if np.iscomplexobj(self.dtype.type(1)):
            b = np.random.random(u.shape)
            u = u + 1j * b
        u = u.astype(self.dtype)
        return a, u, l, r


def get_updated(a, u, downdate):
    if downdate:
        if u.ndim > 1:
            return a - u @ u.T.conjugate()
        return a - np.outer(u, u.conj())
    else:
        if u.ndim > 1:
            return a + u @ u.T.conjugate()
        return a + np.outer(u, u.conj())


class BaseCholUpdate(BaseCholdeltas):
    def base_update1_test(self, a, u, l, r, lower: bool, downdate: bool):
        a1 = get_updated(a, u, downdate)
        if lower:
            d = cholesky_update(l, u, downdate=downdate, lower=True,
                                overwrite_rz=False)
        else:
            d = cholesky_update(r, u, downdate=downdate, lower=False,
                                overwrite_rz=False)
        check_cholesky(a1, d, lower=lower, rtol=self.rtol, atol=self.atol)

        # Check with the same matrix passed to fn to ensure no overwrite
        if lower:
            a2 = np.dot(l, l.T.conj())
        else:
            a2 = np.dot(r.T.conj(), r)
        a2 = get_updated(a2, u, downdate)
        check_cholesky(a2, d, lower=lower, rtol=self.rtol, atol=self.atol)

    def test_u1_l_dd_sqr(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=True)

    def test_u1_l_dd_1x1(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=True)

    def test_u1_u_dd_sqr(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=True)

    def test_u1_u_dd_1x1(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=True)

    def test_u1_l_ud_sqr(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=False)

    def test_u1_l_ud_1x1(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=False)

    def test_u1_u_ud_sqr(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=False)

    def test_u1_u_ud_1x1(self):
        for which in ['1d', 'col']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=False)

    def test_up_l_dd_sqr(self):
        for which in ['p']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=True)

    def test_up_l_dd_1x1(self):
        for which in ['p']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=True)

    def test_up_u_dd_sqr(self):
        for which in ['p']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=True)

    def test_up_u_dd_1x1(self):
        for which in ['p']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=True)

    def test_up_l_ud_sqr(self):
        for which in ['p']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=False)

    def test_up_l_ud_1x1(self):
        for which in ['p']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=True, downdate=False)

    def test_up_u_ud_sqr(self):
        for which in ['p']:
            a, u, l, r = self.generate('sqr', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=False)

    def test_up_u_ud_1x1(self):
        for which in ['p']:
            a, u, l, r = self.generate('1x1', which=which)
            self.base_update1_test(a, u, l, r, lower=False, downdate=False)

    def test_list_input(self):
        a, u, l, r = self.generate('sqr', '1d')

        d = cholesky_update(l, list(u), downdate=False, lower=True,
                            overwrite_rz=False)
        check_cholesky(get_updated(a, u, downdate=False), d, lower=True,
                       rtol=self.rtol, atol=self.atol)

        d = cholesky_update(list(l), u, downdate=False, lower=True,
                            overwrite_rz=False)
        check_cholesky(get_updated(a, u, downdate=False), d, lower=True,
                       rtol=self.rtol, atol=self.atol)

        d = cholesky_update(list(l), list(u), downdate=False, lower=True,
                            overwrite_rz=False)
        check_cholesky(get_updated(a, u, downdate=False), d, lower=True,
                       rtol=self.rtol, atol=self.atol)

    def test_check_nan(self):
        a, u0, l0, r = self.generate('sqr', '1d')

        l = l0.copy('F')
        l[1, 1] = np.nan
        assert_raises(ValueError, cholesky_update, l, u0, False, False)
        assert_raises(ValueError, cholesky_update, l, u0, False, True)
        assert_raises(ValueError, cholesky_update, l, u0, True, False)
        assert_raises(ValueError, cholesky_update, l, u0, True, True)

        u = u0.copy('F')
        u[1] = np.nan
        assert_raises(ValueError, cholesky_update, l0, u, False, False)
        assert_raises(ValueError, cholesky_update, l0, u, False, True)
        assert_raises(ValueError, cholesky_update, l0, u, True, False)
        assert_raises(ValueError, cholesky_update, l0, u, True, True)

    def test_check_nan_p(self):
        a, u0, l0, r = self.generate('sqr', 'p')

        l = l0.copy('F')
        l[1, 1] = np.nan
        assert_raises(ValueError, cholesky_update, l, u0, False, False)
        assert_raises(ValueError, cholesky_update, l, u0, False, True)
        assert_raises(ValueError, cholesky_update, l, u0, True, False)
        assert_raises(ValueError, cholesky_update, l, u0, True, True)

        u = u0.copy('F')
        u[0, 0] = np.nan
        assert_raises(ValueError, cholesky_update, l0, u, False, False)
        assert_raises(ValueError, cholesky_update, l0, u, False, True)
        assert_raises(ValueError, cholesky_update, l0, u, True, False)
        assert_raises(ValueError, cholesky_update, l0, u, True, True)

    def test_check_finite(self):
        a, u0, l0, r = self.generate('sqr', '1d')

        l = l0.copy('F')
        l[1, 1] = np.inf
        assert_raises(ValueError, cholesky_update, l, u0, False, False)
        assert_raises(ValueError, cholesky_update, l, u0, False, True)
        assert_raises(ValueError, cholesky_update, l, u0, True, False)
        assert_raises(ValueError, cholesky_update, l, u0, True, True)

        u = u0.copy('F')
        u[1] = np.inf
        assert_raises(ValueError, cholesky_update, l0, u, False, False)
        assert_raises(ValueError, cholesky_update, l0, u, False, True)
        assert_raises(ValueError, cholesky_update, l0, u, True, False)
        assert_raises(ValueError, cholesky_update, l0, u, True, True)

    def test_check_finite_p(self):
        a, u0, l0, r = self.generate('sqr', 'p')

        l = l0.copy('F')
        l[1, 1] = np.inf
        assert_raises(ValueError, cholesky_update, l, u0, False, False)
        assert_raises(ValueError, cholesky_update, l, u0, False, True)
        assert_raises(ValueError, cholesky_update, l, u0, True, False)
        assert_raises(ValueError, cholesky_update, l, u0, True, True)

        u = u0.copy('F')
        u[1, 0] = np.inf
        assert_raises(ValueError, cholesky_update, l0, u, False, False)
        assert_raises(ValueError, cholesky_update, l0, u, False, True)
        assert_raises(ValueError, cholesky_update, l0, u, True, False)
        assert_raises(ValueError, cholesky_update, l0, u, True, True)

    def base_non_simple_strides(self, adjust_strides, which, r_overwriteable,
                                u_overwriteable):
        for type in ['sqr', '1x1']:
            a, u0, l0, r = self.generate(type, which=which)
            us, ls = adjust_strides((u0, l0))
            a1 = get_updated(a, u0, downdate=False)

            # if overwriteable is True
            # For l and u we try first with overwrite=False, then with
            # overwrite=True and make sure that overwriting occurred
            # appropriately.
            u = u0.copy('F')
            l = l0.copy('F')
            d1 = cholesky_update(ls, u, False, lower=True, overwrite_rz=False)
            check_cholesky(a1, d1, lower=True, rtol=self.rtol, atol=self.atol)
            d1o = cholesky_update(ls, u, False, lower=True, overwrite_rz=True)
            check_cholesky(a1, d1o, lower=True, rtol=self.rtol, atol=self.atol)
            if r_overwriteable:
                assert_allclose(d1o, ls, rtol=self.rtol, atol=self.atol)
            elif type != '1x1':
                assert_raises(AssertionError, assert_allclose, d1o, ls,
                              rtol=self.rtol, atol=self.atol)
            if u_overwriteable and len(u0) != 1:
                # if u contains single element it won't be modified.
                assert_raises(AssertionError, assert_allclose, u0, u,
                              self.rtol, self.atol)

            u = u0.copy('F')
            l = l0.copy('F')
            d2 = cholesky_update(l, us, False, lower=True, overwrite_rz=False)
            check_cholesky(a1, d2, lower=True, rtol=self.rtol, atol=self.atol)
            d2o = cholesky_update(l, us, False, lower=True, overwrite_rz=True)
            check_cholesky(a1, d2o, lower=True, rtol=self.rtol, atol=self.atol)
            if r_overwriteable:
                assert_allclose(d2o, l, rtol=self.rtol, atol=self.atol)
            elif type != '1x1':
                assert_raises(AssertionError, assert_allclose, d2o, ls,
                              rtol=self.rtol, atol=self.atol)
            if u_overwriteable and len(u0) != 1:
                assert_raises(AssertionError, assert_allclose, u0, us,
                              rtol=self.rtol, atol=self.atol)

            # Now with both strided (need to be recreated since they were
            # overwritten above)
            u = u0.copy('F')
            l = l0.copy('F')
            # since some of these were consumed above
            us, ls = adjust_strides((u, l))
            d3 = cholesky_update(ls, us, False, lower=True, overwrite_rz=False)
            check_cholesky(a1, d3, lower=True, rtol=self.rtol, atol=self.atol)
            d3o = cholesky_update(ls, us, False, lower=True, overwrite_rz=True)
            check_cholesky(a1, d3o, lower=True, rtol=self.rtol, atol=self.atol)
            if r_overwriteable:
                assert_allclose(d3o, ls, rtol=self.rtol, atol=self.atol)
            elif type != '1x1':
                assert_raises(AssertionError, assert_allclose, d3o, ls,
                              rtol=self.rtol, atol=self.atol)
            if u_overwriteable and len(u0) != 1:
                assert_raises(AssertionError, assert_allclose, u0, us,
                              rtol=self.rtol, atol=self.atol)

    def test_non_unit_strides_1d(self):
        self.base_non_simple_strides(make_strided, '1d', True, True)

    def test_neg_strides_1d(self):
        self.base_non_simple_strides(negate_strides, '1d', False, False)

    def test_non_itemize_strides_1d(self):
        self.base_non_simple_strides(nonitemsize_strides, '1d', False, False)

    def test_non_native_byte_order_1d(self):
        self.base_non_simple_strides(make_nonnative, '1d', False, False)

    def test_non_unit_strides_col(self):
        self.base_non_simple_strides(make_strided, 'col', True, True)

    def test_neg_strides_col(self):
        self.base_non_simple_strides(negate_strides, 'col', False, False)

    def test_non_itemize_strides_col(self):
        self.base_non_simple_strides(nonitemsize_strides, 'col', False, False)

    def test_non_native_byte_order_col(self):
        self.base_non_simple_strides(make_nonnative, 'col', False, False)

    def test_non_unit_strides_p(self):
        self.base_non_simple_strides(make_strided, 'p', True, False)

    def test_neg_strides_p(self):
        self.base_non_simple_strides(negate_strides, 'p', False, False)

    def test_non_itemize_strides_p(self):
        self.base_non_simple_strides(nonitemsize_strides, 'p', False, False)

    def test_non_native_byte_order_p(self):
        self.base_non_simple_strides(make_nonnative, 'p', False, False)

    def base_overwrite_chol(self, which):
        for lower, downdate, type in itertools.product(
                [True, False], [True, False], ['sqr', '1x1']):
            a, u0, l0, r0 = self.generate(type, which=which)
            a1 = get_updated(a, u0, downdate)
            mat0 = l0 if lower else r0

            # Don't overwrite
            u = u0.copy('F')
            mat = mat0.copy('F')
            d1 = cholesky_update(mat, u, downdate, lower, overwrite_rz=False)
            check_cholesky(a1, d1, lower=lower, rtol=self.rtol, atol=self.atol)
            check_cholesky(a, mat, lower=lower, rtol=self.rtol, atol=self.atol)

            u = u0.copy('F')
            mat = mat0.copy('F')
            d2 = cholesky_update(mat, u, downdate, lower, overwrite_rz=True)
            check_cholesky(a1, d2, lower=lower, rtol=self.rtol, atol=self.atol)
            assert_allclose(d2, mat, rtol=self.rtol, atol=self.atol)
            # Check that the untouched triangle remains the same
            assert_allclose_tri(mat0, d2, lower=not lower, rtol=self.rtol,
                                atol=self.atol)
            if len(u) > 1:
                assert_raises(AssertionError, assert_allclose, u, u0,
                              rtol=self.rtol, atol=self.atol)

            u = u0.copy('C')
            mat = mat0.copy('C')
            d3 = cholesky_update(mat, u, downdate, lower, overwrite_rz=True)
            check_cholesky(a1, d3, lower=lower, rtol=self.rtol, atol=self.atol)
            assert_allclose(d3, mat, rtol=self.rtol, atol=self.atol)
            assert_allclose_tri(mat0, d3, lower=not lower, rtol=self.rtol,
                                atol=self.atol)
            if len(u) > 1:
                assert_raises(AssertionError, assert_allclose, u, u0,
                              rtol=self.rtol, atol=self.atol)

            u = u0.copy('C')
            mat = mat0.copy('F')
            d4 = cholesky_update(mat, u, downdate, lower, overwrite_rz=True)
            check_cholesky(a1, d4, lower=lower, rtol=self.rtol, atol=self.atol)
            assert_allclose(d4, mat, rtol=self.rtol, atol=self.atol)
            # Check that the untouched triangle remains the same
            assert_allclose_tri(mat0, d4, lower=not lower, rtol=self.rtol,
                                atol=self.atol)
            if len(u) > 1:
                assert_raises(AssertionError, assert_allclose, u, u0,
                              rtol=self.rtol, atol=self.atol)

            u = u0.copy('F')
            mat = mat0.copy('C')
            d5 = cholesky_update(mat, u, downdate, lower, overwrite_rz=True)
            check_cholesky(a1, d5, lower=lower, rtol=self.rtol, atol=self.atol)
            assert_allclose(d5, mat, rtol=self.rtol, atol=self.atol)
            # Check that the untouched triangle remains the same
            assert_allclose_tri(mat0, d5, lower=not lower, rtol=self.rtol,
                                atol=self.atol)
            if len(u) > 1:
                assert_raises(AssertionError, assert_allclose, u, u0,
                              rtol=self.rtol, atol=self.atol)

    def test_overwrite_1d(self):
        self.base_overwrite_chol('1d')

    def test_overwrite_col(self):
        self.base_overwrite_chol('col')

    def test_overwrite_p(self):
        self.base_overwrite_chol('p')

    """ Test edge cases which should error out """
    def test_empty_l(self):
        a, u, l, r = self.generate('sqr', which='col')
        assert_raises(ValueError, cholesky_update, np.array([]), u)

    def test_empty_u(self):
        a, u, l, r = self.generate('sqr', which='col')
        assert_raises(ValueError, cholesky_update, l, np.array([]))

    def test_unsupported_dtypes(self):
        dts = ['int8', 'int16', 'int32', 'int64',
               'uint8', 'uint16', 'uint32', 'uint64',
               'float16', 'longdouble', 'longcomplex',
               'bool']
        a, u0, l0, r = self.generate('sqr', which='col')
        for dtype in dts:
            l = l0.real.astype(dtype)
            u = u0.real.astype(dtype)
            assert_raises(ValueError, cholesky_update, l, u0)
            assert_raises(ValueError, cholesky_update, l, u0)
            assert_raises(ValueError, cholesky_update, l, u0)
            assert_raises(ValueError, cholesky_update, l, u0)

            assert_raises(ValueError, cholesky_update, l0, u)
            assert_raises(ValueError, cholesky_update, l0, u)
            assert_raises(ValueError, cholesky_update, l0, u)
            assert_raises(ValueError, cholesky_update, l0, u)

    def test_scalar(self):
        a, u0, l0, r = self.generate('1x1', which='col')
        assert_raises(ValueError, cholesky_update, l0[0, 0], u0)
        assert_raises(ValueError, cholesky_update, l0, u0[0, 0])

    def test_mismatched_shapes(self):
        a, u, l, r = self.generate('sqr', which='col')
        a0, u0, l0, r0 = self.generate('1x1', which='col')

        assert_raises(ValueError, cholesky_update, l, u0)
        assert_raises(ValueError, cholesky_update, l0, u)

    def test_l_nonsquare(self):
        a, u, l0, r = self.generate('sqr', which='col')
        l = np.concatenate((l0, l0[0, None]), axis=0)
        assert_raises(ValueError, cholesky_update, l, u)

        l = np.concatenate((l0, l0[:, 0][:, None]), axis=1)
        assert_raises(ValueError, cholesky_update, l, u)

    def test_3d_inputs(self):
        a, u0, l0, r = self.generate('sqr', which='col')
        l = l0[..., None]
        assert_raises(ValueError, cholesky_update, l, u0)
        u = u0[..., None]
        assert_raises(ValueError, cholesky_update, l0, u)

    def test_non_pd(self):
        a0, u0, l0, r0 = self.generate('sqr', which='col')
        a1 = a0 - (a0.shape[0] - 1.9) * np.eye(a0.shape[0], dtype=self.dtype)
        for lower in [True, False]:
            l1 = scipy.linalg.cholesky(a1, lower=lower)

            u1 = u0 * 5
            a2 = get_updated(a1, u1, downdate=True)
            assert_raises(np.linalg.LinAlgError, scipy.linalg.cholesky, a2,
                          lower=lower)
            assert_raises(np.linalg.LinAlgError, cholesky_update, l1, u1,
                          lower=lower, downdate=True)
            # update should never fail
            cholesky_update(l1, u1, lower=lower, downdate=False)

    def test_non_pd_p(self):
        a0, u0, l0, r0 = self.generate('sqr', which='p')
        a1 = a0 - (a0.shape[0] - 1.9) * np.eye(a0.shape[0], dtype=self.dtype)
        for lower in [False]:

            l1 = scipy.linalg.cholesky(a1, lower=lower)

            u1 = u0 * 1.3
            a2 = get_updated(a1, u1, downdate=True)
            assert_raises(np.linalg.LinAlgError, scipy.linalg.cholesky, a2,
                          lower=lower)
            assert_raises(np.linalg.LinAlgError, cholesky_update, l1, u1,
                          lower=lower, downdate=True)
            # update should never fail
            cholesky_update(l1, u1, lower=lower, downdate=False)


class TestCholUpdate_f(BaseCholUpdate):
    dtype = np.dtype('f')


class TestCholUpdate_d(BaseCholUpdate):
    dtype = np.dtype('d')


class TestCholUpdate_F(BaseCholUpdate):
    dtype = np.dtype('F')


class TestCholUpdate_D(BaseCholUpdate):
    dtype = np.dtype('D')
