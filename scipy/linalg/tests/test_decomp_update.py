import itertools
import numpy as np
from numpy.testing import assert_, assert_allclose, assert_raises
from scipy import linalg
import scipy.linalg._decomp_update as _decomp_update
from scipy.linalg._decomp_update import *

def assert_unitary(a, rtol=None, atol=None, assert_sqr=True):
    if rtol is None:
        rtol = 10.0 ** -(np.finfo(a.dtype).precision-2)
    if atol is None:
        atol = 2*np.finfo(a.dtype).eps

    if assert_sqr:
        assert_(a.shape[0] == a.shape[1], 'unitary matrices must be square')
    aTa = np.dot(a.T.conj(), a)
    assert_allclose(aTa, np.eye(a.shape[1]), rtol=rtol, atol=atol)

def assert_upper_tri(a, rtol=None, atol=None):
    if rtol is None:
        rtol = 10.0 ** -(np.finfo(a.dtype).precision-2)
    if atol is None:
        atol = 2*np.finfo(a.dtype).eps
    mask = np.tri(a.shape[0], a.shape[1], -1, np.bool_)
    assert_allclose(a[mask], 0.0, rtol=rtol, atol=atol)

def check_qr(q, r, a, rtol, atol, assert_sqr=True):
    assert_unitary(q, rtol, atol, assert_sqr)
    assert_upper_tri(r, rtol, atol)
    assert_allclose(q.dot(r), a, rtol=rtol, atol=atol)

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
            base = np.zeros((s[0]*a.shape[0]+s[1], t[0]*a.shape[1]+t[1]), a.dtype)
            view = base[s[1]::s[0], t[1]::t[0]]
            view[...] = a
        else:
            raise ValueError('make_strided only works for ndim = 1 or 2 arrays')
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
            raise ValueError('negate_strides only works for ndim = 1 or 2 arrays')
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
    out = []
    for a in arrs:
        out.append(a.astype(a.dtype.newbyteorder()))
    return out

class BaseQRdeltas(object):
    def __init__(self):
        self.rtol = 10.0 ** -(np.finfo(self.dtype).precision-2)
        self.atol = 5 * np.finfo(self.dtype).eps

    def generate(self, type, mode='full'):
        np.random.seed(29382)
        shape = {'sqr': (8, 8), 'tall': (12, 7), 'fat': (7, 12)}[type]
        a = np.random.random(shape)
        if np.iscomplexobj(self.dtype.type(1)):
            b = np.random.random(shape)
            a = a + 1j * b
        a = a.astype(self.dtype)
        q, r = linalg.qr(a, mode=mode)
        # numpy 1.5.1 np.triu can modify array dtype. So qr of 'f' arrays
        # will return a 'd' r and a 'f' q.
        if r.dtype != self.dtype:
            r = r.astype(self.dtype)
        return a, q, r

class BaseQRdelete(BaseQRdeltas):
    def test_sqr_1_row(self):
        a, q, r = self.generate('sqr')
        for row in range(r.shape[0]):
            q1, r1 = qr_delete(q, r, row, overwrite_qr=False)
            a1 = np.delete(a, row, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_sqr_p_row(self):
        a, q, r = self.generate('sqr')
        for ndel in range(2, 6):
            for row in range(a.shape[0]-ndel):
                q1, r1 = qr_delete(q, r, row, ndel, overwrite_qr=False)
                a1 = np.delete(a, slice(row, row+ndel), 0)
                check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_sqr_1_col(self):
        a, q, r = self.generate('sqr')
        for col in range(r.shape[1]):
            q1, r1 = qr_delete(q, r, col, which='col', overwrite_qr=False)
            a1 = np.delete(a, col, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_sqr_p_col(self):
        a, q, r = self.generate('sqr')
        for ndel in range(2, 6):
            for col in range(r.shape[1]-ndel):
                q1, r1 = qr_delete(q, r, col, ndel, which='col',
                                   overwrite_qr=False)
                a1 = np.delete(a, slice(col, col+ndel), 1)
                check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_1_row(self):
        a, q, r = self.generate('tall')
        for row in range(r.shape[0]):
            q1, r1 = qr_delete(q, r, row, overwrite_qr=False)
            a1 = np.delete(a, row, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_p_row(self):
        a, q, r = self.generate('tall')
        for ndel in range(2, 6):
            for row in range(a.shape[0]-ndel):
                q1, r1 = qr_delete(q, r, row, ndel, overwrite_qr=False)
                a1 = np.delete(a, slice(row, row+ndel), 0)
                check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_1_col(self):
        a, q, r = self.generate('tall')
        for col in range(r.shape[1]):
            q1, r1 = qr_delete(q, r, col, which='col', overwrite_qr=False)
            a1 = np.delete(a, col, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_p_col(self):
        a, q, r = self.generate('tall')
        for ndel in range(2, 6):
            for col in range(r.shape[1]-ndel):
                q1, r1 = qr_delete(q, r, col, ndel, which='col',
                                   overwrite_qr=False)
                a1 = np.delete(a, slice(col, col+ndel), 1)
                check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_1_row(self):
        a, q, r = self.generate('fat')
        for row in range(r.shape[0]):
            q1, r1 = qr_delete(q, r, row, overwrite_qr=False)
            a1 = np.delete(a, row, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_p_row(self):
        a, q, r = self.generate('fat')
        for ndel in range(2, 6):
            for row in range(a.shape[0]-ndel):
                q1, r1 = qr_delete(q, r, row, ndel, overwrite_qr=False)
                a1 = np.delete(a, slice(row, row+ndel), 0)
                check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_1_col(self):
        a, q, r = self.generate('fat')
        for col in range(r.shape[1]):
            q1, r1 = qr_delete(q, r, col, which='col', overwrite_qr=False)
            a1 = np.delete(a, col, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_p_col(self):
        a, q, r = self.generate('fat')
        for ndel in range(2, 6):
            for col in range(r.shape[1]-ndel):
                q1, r1 = qr_delete(q, r, col, ndel, which='col',
                                   overwrite_qr=False)
                a1 = np.delete(a, slice(col, col+ndel), 1)
                check_qr(q1, r1, a1, self.rtol, self.atol)

    # all row deletes and single column deletes should be able to
    # handle any non negative strides. (only row and column vector
    # operations are used.) p column delete require fortran ordered
    # Q and R and will make a copy as necessary.

    def base_non_simple_strides(self, adjust_strides, ks, p, which, overwriteable):
        if which == 'row':
            qind = (slice(p,None), slice(p,None))
            rind = (slice(p,None), slice(None))
        else:
            qind = (slice(None), slice(None))
            rind = (slice(None), slice(None,-p))

        for type, k in itertools.product(['sqr', 'tall', 'fat'], ks):
            a, q0, r0, = self.generate(type)
            qs, rs = adjust_strides((q0, r0))
            if p == 1:
                a1 = np.delete(a, k, 0 if which == 'row' else 1)
            else:
                s = slice(k,k+p)
                if k < 0:
                    s = slice(k, k+p+ (a.shape[0] if which == 'row' else a.shape[1]))
                a1 = np.delete(a, s, 0 if which == 'row' else 1)

            # for each variable, q, r we try with it strided and
            # overwrite=False. Then we try with overwrite=True, and make
            # sure that q and r are still overwritten.

            q = q0.copy('F'); r = r0.copy('F') 
            q1, r1 = qr_delete(qs, r, k, p, which, False)
            check_qr(q1, r1, a1, self.rtol, self.atol)
            q1o, r1o = qr_delete(qs, r, k, p, which, True)
            check_qr(q1o, r1o, a1, self.rtol, self.atol)
            if overwriteable:
                assert_allclose(q1o, qs[qind], rtol=self.rtol, atol=self.atol)
                assert_allclose(r1o, r[rind], rtol=self.rtol, atol=self.atol)

            q = q0.copy('F'); r = r0.copy('F')
            q2, r2 = qr_delete(q, rs, k, p, which, False)
            check_qr(q2, r2, a1, self.rtol, self.atol)
            q2o, r2o = qr_delete(q, rs, k, p, which, True)
            check_qr(q2o, r2o, a1, self.rtol, self.atol)
            if overwriteable:
                assert_allclose(q2o, q[qind], rtol=self.rtol, atol=self.atol)
                assert_allclose(r2o, rs[rind], rtol=self.rtol, atol=self.atol)

            q = q0.copy('F'); r = r0.copy('F')
            # since some of these were consumed above
            qs, rs = adjust_strides((q, r))
            q3, r3 = qr_delete(qs, rs, k, p, which, False)
            check_qr(q3, r3, a1, self.rtol, self.atol)
            q3o, r3o = qr_delete(qs, rs, k, p, which, True)
            check_qr(q3o, r3o, a1, self.rtol, self.atol)
            if overwriteable:
                assert_allclose(q2o, qs[qind], rtol=self.rtol, atol=self.atol)
                assert_allclose(r3o, rs[rind], rtol=self.rtol, atol=self.atol)

    def test_non_unit_strides_1_row(self):
        self.base_non_simple_strides(make_strided, [0], 1, 'row', True)

    def test_non_unit_strides_p_row(self):
        self.base_non_simple_strides(make_strided, [0], 3, 'row', True)

    def test_non_unit_strides_1_col(self):
        self.base_non_simple_strides(make_strided, [0], 1, 'col', True)

    def test_non_unit_strides_p_col(self):
        self.base_non_simple_strides(make_strided, [0], 3, 'col', False)

    def test_neg_strides_1_row(self):
        self.base_non_simple_strides(negate_strides, [0], 1, 'row', False)

    def test_neg_strides_p_row(self):
        self.base_non_simple_strides(negate_strides, [0], 3, 'row', False)

    def test_neg_strides_1_col(self):
        self.base_non_simple_strides(negate_strides, [0], 1, 'col', False)

    def test_neg_strides_p_col(self):
        self.base_non_simple_strides(negate_strides, [0], 3, 'col', False)

    def test_non_itemize_strides_1_row(self):
        self.base_non_simple_strides(nonitemsize_strides, [0], 1, 'row', False)

    def test_non_itemize_strides_p_row(self):
        self.base_non_simple_strides(nonitemsize_strides, [0], 3, 'row', False)

    def test_non_itemize_strides_1_col(self):
        self.base_non_simple_strides(nonitemsize_strides, [0], 1, 'col', False)

    def test_non_itemize_strides_p_col(self):
        self.base_non_simple_strides(nonitemsize_strides, [0], 3, 'col', False)

    def test_non_native_byte_order_1_row(self):
        self.base_non_simple_strides(make_nonnative, [0], 1, 'row', False)

    def test_non_native_byte_order_p_row(self):
        self.base_non_simple_strides(make_nonnative, [0], 3, 'row', False)

    def test_non_native_byte_order_1_col(self):
        self.base_non_simple_strides(make_nonnative, [0], 1, 'col', False)

    def test_non_native_byte_order_p_col(self):
        self.base_non_simple_strides(make_nonnative, [0], 3, 'col', False)

    def test_neg_k(self):
        a, q, r = self.generate('sqr')
        for k, p, w in itertools.product([-3, -7], [1, 3], ['row', 'col']):
            q1, r1 = qr_delete(q, r, k, p, w, overwrite_qr=False)
            if w == 'row':
                a1 = np.delete(a, slice(k+a.shape[0], k+p+a.shape[0]), 0)
            else:
                a1 = np.delete(a, slice(k+a.shape[0], k+p+a.shape[1]), 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def base_overwrite_qr(self, which, p, test_C, test_F):
        if which == 'row':
            qind = (slice(p,None), slice(p,None))
            rind = (slice(p,None), slice(None))
        else:
            qind = (slice(None), slice(None))
            rind = (slice(None), slice(None,-p))
        a, q0, r0 = self.generate('sqr')
        if p == 1:
            a1 = np.delete(a, 3, 0 if which == 'row' else 1)
        else:
            a1 = np.delete(a, slice(3, 3+p), 0 if which == 'row' else 1)

        # don't overwrite
        q = q0.copy('F'); r = r0.copy('F')
        q1, r1 = qr_delete(q, r, 3, p, which, False)
        check_qr(q1, r1, a1, self.rtol, self.atol)
        check_qr(q, r, a, self.rtol, self.atol)

        if test_F:
            q = q0.copy('F'); r = r0.copy('F')
            q2, r2 = qr_delete(q, r, 3, p, which, True)
            check_qr(q2, r2, a1, self.rtol, self.atol)
            # verify the overwriting
            assert_allclose(q2, q[qind], rtol=self.rtol, atol=self.atol)
            assert_allclose(r2, r[rind], rtol=self.rtol, atol=self.atol)

        if test_C:
            q = q0.copy('C'); r = r0.copy('C')
            q3, r3 = qr_delete(q, r, 3, p, which, True)
            check_qr(q3, r3, a1, self.rtol, self.atol)
            assert_allclose(q3, q[qind], rtol=self.rtol, atol=self.atol)
            assert_allclose(r3, r[rind], rtol=self.rtol, atol=self.atol)

    def test_overwrite_qr_1_row(self):
        # any positively strided q and r.
        self.base_overwrite_qr('row', 1, True, True)

    def test_overwrite_qr_1_col(self):
        # any positively strided q and r.
        self.base_overwrite_qr('col', 1, True, True)

    def test_overwrite_qr_p_row(self):
        # any positively strided q and r.
        self.base_overwrite_qr('row', 3, True, True)

    def test_overwrite_qr_p_col(self):
        # only F orderd q and r can be overwritten for cols
        self.base_overwrite_qr('col', 3, False, True)

    def test_economic_qr(self):
        a, q, r = self.generate('tall', mode='economic')
        # only test row delete, this logic is shared.
        assert_raises(ValueError, qr_delete, q, r, 0)

    def test_bad_which(self):
        a, q, r = self.generate('sqr')
        assert_raises(ValueError, qr_delete, q, r, 0, which='foo')

    def test_bad_k(self):
        a, q, r = self.generate('tall')
        assert_raises(ValueError, qr_delete, q, r, q.shape[0], 1)
        assert_raises(ValueError, qr_delete, q, r, -q.shape[0]-1, 1)
        assert_raises(ValueError, qr_delete, q, r, r.shape[0], 1, 'col')
        assert_raises(ValueError, qr_delete, q, r, -r.shape[0]-1, 1, 'col')

    def test_bad_p(self):
        a, q, r = self.generate('tall')
        # p must be positive
        assert_raises(ValueError, qr_delete, q, r, 0, -1)
        assert_raises(ValueError, qr_delete, q, r, 0, -1, 'col')

        # and nonzero
        assert_raises(ValueError, qr_delete, q, r, 0, 0)
        assert_raises(ValueError, qr_delete, q, r, 0, 0, 'col')

        # must have at least k+p rows or cols, depending.
        assert_raises(ValueError, qr_delete, q, r, 3, q.shape[0]-2)
        assert_raises(ValueError, qr_delete, q, r, 3, r.shape[1]-2, 'col')

    def test_empty_q(self):
        a, q, r = self.generate('tall')
        # same code path for 'row' and 'col'
        assert_raises(ValueError, qr_delete, np.array([]), r, 0, 1)

    def test_empty_r(self):
        a, q, r = self.generate('tall')
        # same code path for 'row' and 'col'
        assert_raises(ValueError, qr_delete, q, np.array([]), 0, 1)

    def test_mismatched_q_and_r(self):
        a, q, r = self.generate('tall')
        r = r[1:]
        assert_raises(ValueError, qr_delete, q, r, 0, 1)

    def test_integer_input(self):
        q = np.arange(16).reshape(4, 4)
        r = q.copy()  # doesn't matter
        assert_raises(ValueError, qr_delete, q, r, 0, 1)

class TestQRdelete_f(BaseQRdelete):
    dtype = np.dtype('f')

class TestQRdelete_F(BaseQRdelete):
    dtype = np.dtype('F')

class TestQRdelete_d(BaseQRdelete):
    dtype = np.dtype('d')

class TestQRdelete_D(BaseQRdelete):
    dtype = np.dtype('D')


