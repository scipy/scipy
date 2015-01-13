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
        self.atol = 10 * np.finfo(self.dtype).eps

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

    def test_economic_1_row(self):
        # this test always starts and ends with an economic decomp.
        a, q, r = self.generate('tall', 'economic')
        for row in range(r.shape[0]):
            q1, r1 = qr_delete(q, r, row, overwrite_qr=False)
            a1 = np.delete(a, row, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    # for economic row deletes 
    # eco - prow = eco
    # eco - prow = sqr
    # eco - prow = fat
    def base_economic_p_row_xxx(self, ndel):
        a, q, r = self.generate('tall', 'economic')
        for row in range(a.shape[0]-ndel):
            q1, r1 = qr_delete(q, r, row, ndel, overwrite_qr=False)
            a1 = np.delete(a, slice(row, row+ndel), 0)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def test_economic_p_row_economic(self):
        # (12, 7) - (3, 7) = (9,7) --> stays economic
        self.base_economic_p_row_xxx(3)

    def test_economic_p_row_sqr(self):
        # (12, 7) - (5, 7) = (7, 7) --> becomes square
        self.base_economic_p_row_xxx(5)

    def test_economic_p_row_fat(self):
        # (12, 7) - (7,7) = (5, 7) --> becomes fat
        self.base_economic_p_row_xxx(7)

    def test_economic_1_col(self):
        a, q, r = self.generate('tall', 'economic')
        for col in range(r.shape[1]):
            q1, r1 = qr_delete(q, r, col, which='col', overwrite_qr=False)
            a1 = np.delete(a, col, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def test_economic_p_col(self):
        a, q, r = self.generate('tall', 'economic')
        for ndel in range(2, 6):
            for col in range(r.shape[1]-ndel):
                q1, r1 = qr_delete(q, r, col, ndel, which='col',
                                   overwrite_qr=False)
                a1 = np.delete(a, slice(col, col+ndel), 1)
                check_qr(q1, r1, a1, self.rtol, self.atol, False)

    # all full qr, row deletes and single column deletes should be able to
    # handle any non negative strides. (only row and column vector
    # operations are used.) p column delete require fortran ordered
    # Q and R and will make a copy as necessary.  Economic qr row deletes 
    # requre a contigous q.

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
                    s = slice(k, k + p + (a.shape[0] if which == 'row' else a.shape[1]))
                a1 = np.delete(a, s, 0 if which == 'row' else 1)

            # for each variable, q, r we try with it strided and
            # overwrite=False. Then we try with overwrite=True, and make
            # sure that q and r are still overwritten.

            q = q0.copy('F')
            r = r0.copy('F') 
            q1, r1 = qr_delete(qs, r, k, p, which, False)
            check_qr(q1, r1, a1, self.rtol, self.atol)
            q1o, r1o = qr_delete(qs, r, k, p, which, True)
            check_qr(q1o, r1o, a1, self.rtol, self.atol)
            if overwriteable:
                assert_allclose(q1o, qs[qind], rtol=self.rtol, atol=self.atol)
                assert_allclose(r1o, r[rind], rtol=self.rtol, atol=self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            q2, r2 = qr_delete(q, rs, k, p, which, False)
            check_qr(q2, r2, a1, self.rtol, self.atol)
            q2o, r2o = qr_delete(q, rs, k, p, which, True)
            check_qr(q2o, r2o, a1, self.rtol, self.atol)
            if overwriteable:
                assert_allclose(q2o, q[qind], rtol=self.rtol, atol=self.atol)
                assert_allclose(r2o, rs[rind], rtol=self.rtol, atol=self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
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

    def base_overwrite_qr(self, which, p, test_C, test_F, mode='full'):
        assert_sqr = True if mode == 'full' else False
        if which == 'row':
            qind = (slice(p,None), slice(p,None))
            rind = (slice(p,None), slice(None))
        else:
            qind = (slice(None), slice(None))
            rind = (slice(None), slice(None,-p))
        a, q0, r0 = self.generate('sqr', mode)
        if p == 1:
            a1 = np.delete(a, 3, 0 if which == 'row' else 1)
        else:
            a1 = np.delete(a, slice(3, 3+p), 0 if which == 'row' else 1)

        # don't overwrite
        q = q0.copy('F')
        r = r0.copy('F')
        q1, r1 = qr_delete(q, r, 3, p, which, False)
        check_qr(q1, r1, a1, self.rtol, self.atol, assert_sqr)
        check_qr(q, r, a, self.rtol, self.atol, assert_sqr)

        if test_F:
            q = q0.copy('F')
            r = r0.copy('F')
            q2, r2 = qr_delete(q, r, 3, p, which, True)
            check_qr(q2, r2, a1, self.rtol, self.atol, assert_sqr)
            # verify the overwriting
            assert_allclose(q2, q[qind], rtol=self.rtol, atol=self.atol)
            assert_allclose(r2, r[rind], rtol=self.rtol, atol=self.atol)

        if test_C:
            q = q0.copy('C')
            r = r0.copy('C')
            q3, r3 = qr_delete(q, r, 3, p, which, True)
            check_qr(q3, r3, a1, self.rtol, self.atol, assert_sqr)
            assert_allclose(q3, q[qind], rtol=self.rtol, atol=self.atol)
            assert_allclose(r3, r[rind], rtol=self.rtol, atol=self.atol)

    def test_overwrite_qr_1_row(self):
        # any positively strided q and r.
        self.base_overwrite_qr('row', 1, True, True)

    def test_overwrite_economic_qr_1_row(self):
        # Any contiguous q and positively strided r.
        self.base_overwrite_qr('row', 1, True, True, 'economic')

    def test_overwrite_qr_1_col(self):
        # any positively strided q and r.
        # full and eco share code paths
        self.base_overwrite_qr('col', 1, True, True)

    def test_overwrite_qr_p_row(self):
        # any positively strided q and r.
        self.base_overwrite_qr('row', 3, True, True)

    def test_overwrite_economic_qr_p_row(self):
        # any contiguous q and positively strided r
        self.base_overwrite_qr('row', 3, True, True, 'economic')

    def test_overwrite_qr_p_col(self):
        # only F orderd q and r can be overwritten for cols
        # full and eco share code paths 
        self.base_overwrite_qr('col', 3, False, True)

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

    def test_check_finite(self):
        a0, q0, r0 = self.generate('tall')

        q = q0.copy('F')
        q[1,1] = np.nan
        assert_raises(ValueError, qr_delete, q, r0, 0, 1, 'row')
        assert_raises(ValueError, qr_delete, q, r0, 0, 3, 'row')
        assert_raises(ValueError, qr_delete, q, r0, 0, 1, 'col')
        assert_raises(ValueError, qr_delete, q, r0, 0, 3, 'col')

        r = r0.copy('F')
        r[1,1] = np.nan
        assert_raises(ValueError, qr_delete, q0, r, 0, 1, 'row')
        assert_raises(ValueError, qr_delete, q0, r, 0, 3, 'row')
        assert_raises(ValueError, qr_delete, q0, r, 0, 1, 'col')
        assert_raises(ValueError, qr_delete, q0, r, 0, 3, 'col')

class TestQRdelete_f(BaseQRdelete):
    dtype = np.dtype('f')

class TestQRdelete_F(BaseQRdelete):
    dtype = np.dtype('F')

class TestQRdelete_d(BaseQRdelete):
    dtype = np.dtype('d')

class TestQRdelete_D(BaseQRdelete):
    dtype = np.dtype('D')

class BaseQRinsert(BaseQRdeltas):
    def generate(self, type, mode='full', which='row', p=1):
        a, q, r = super(BaseQRinsert, self).generate(type, mode)

        assert_(p > 0)

        # super call set the seed...
        if which == 'row':
            u = np.random.random((p, a.shape[1]))
        elif which == 'col':
            u = np.random.random((a.shape[0], p))
        else:
            ValueError('which should be either "row" or "col"')

        u = np.squeeze(u)

        if np.iscomplexobj(self.dtype.type(1)):
            b = np.random.random(u.shape)
            u = u + 1j * b

        u = u.astype(self.dtype)
        return a, q, r, u

    def test_sqr_1_row(self):
        a, q, r, u = self.generate('sqr', which='row')
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row)
            a1 = np.insert(a, row, u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)
    
    def test_sqr_p_row(self):
        # sqr + rows --> fat always
        a, q, r, u = self.generate('sqr', which='row', p=3)
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row)
            a1 = np.insert(a, row*np.ones(3, np.intp), u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_sqr_1_col(self):
        a, q, r, u = self.generate('sqr', which='col')
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col, u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_sqr_p_col(self):
        # sqr + cols --> fat always
        a, q, r, u = self.generate('sqr', which='col', p=3)
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col*np.ones(3, np.intp), u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_1_row(self):
        a, q, r, u = self.generate('tall', which='row')
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row)
            a1 = np.insert(a, row, u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_p_row(self):
        # tall + rows --> tall always
        a, q, r, u = self.generate('tall', which='row', p=3)
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row)
            a1 = np.insert(a, row*np.ones(3, np.intp), u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_1_col(self):
        a, q, r, u = self.generate('tall', which='col')
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col, u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    # for column adds to tall matrices there are three cases to test
    # tall + pcol --> tall
    # tall + pcol --> sqr
    # tall + pcol --> fat
    def base_tall_p_col_xxx(self, p):
        a, q, r, u = self.generate('tall', which='col', p=p)
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col*np.ones(p, np.intp), u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_p_col_tall(self):
        # 12x7 + 12x3 = 12x10 --> stays tall
        self.base_tall_p_col_xxx(3)

    def test_tall_p_col_sqr(self):
        # 12x7 + 12x5 = 12x12 --> becomes sqr
        self.base_tall_p_col_xxx(5)

    def test_tall_p_col_fat(self):
        # 12x7 + 12x7 = 12x14 --> becomes fat
        self.base_tall_p_col_xxx(7)

    def test_fat_1_row(self):
        a, q, r, u = self.generate('fat', which='row')
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row)
            a1 = np.insert(a, row, u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    # for row adds to fat matrices there are three cases to test
    # fat + prow --> fat
    # fat + prow --> sqr
    # fat + prow --> tall
    def base_fat_p_row_xxx(self, p):
        a, q, r, u = self.generate('fat', which='row', p=p)
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row)
            a1 = np.insert(a, row*np.ones(p, np.intp), u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol)
    
    def test_fat_p_row_fat(self):
        # 7x12 + 3x12 = 10x12 --> stays fat
        self.base_fat_p_row_xxx(3)

    def test_fat_p_row_sqr(self):
        # 7x12 + 5x12 = 12x12 --> becomes sqr
        self.base_fat_p_row_xxx(5)

    def test_fat_p_row_tall(self):
        # 7x12 + 7x12 = 14x12 --> becomes tall
        self.base_fat_p_row_xxx(7)

    def test_fat_1_col(self):
        a, q, r, u = self.generate('fat', which='col')
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col, u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_p_col(self):
        # fat + cols --> fat always
        a, q, r, u = self.generate('fat', which='col', p=3)
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col*np.ones(3, np.intp), u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_economic_1_row(self):
        a, q, r, u = self.generate('tall', 'economic', 'row')
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row, overwrite_qru=False)
            a1 = np.insert(a, row, u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def test_economic_p_row(self):
        # tall + rows --> tall always
        a, q, r, u = self.generate('tall', 'economic', 'row', 3)
        for row in range(r.shape[0]):
            q1, r1 = qr_insert(q, r, u, row, overwrite_qru=False)
            a1 = np.insert(a, row*np.ones(3, np.intp), u, 0)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def test_economic_1_col(self):
        a, q, r, u = self.generate('tall', 'economic',  which='col')
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u.copy(), col, 'col', False)
            a1 = np.insert(a, col, u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    # for column adds to economic matrices there are three cases to test
    # eco + pcol --> eco
    # eco + pcol --> sqr
    # eco + pcol --> fat
    def base_economic_p_col_xxx(self, p):
        a, q, r, u = self.generate('tall', 'economic', which='col', p=p)
        for col in range(r.shape[1]):
            q1, r1 = qr_insert(q, r, u, col, 'col', False)
            a1 = np.insert(a, col*np.ones(p, np.intp), u, 1)
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def test_economic_p_col_eco(self):
        # 12x7 + 12x3 = 12x10 --> stays eco
        self.base_economic_p_col_xxx(3)

    def test_economic_p_col_sqr(self):
        # 12x7 + 12x5 = 12x12 --> becomes sqr
        self.base_economic_p_col_xxx(5)

    def test_economic_p_col_fat(self):
        # 12x7 + 12x7 = 12x14 --> becomes fat
        self.base_economic_p_col_xxx(7)

    def base_non_simple_strides(self, adjust_strides, k, p, which):
        for type in ['sqr', 'tall', 'fat']:
            a, q0, r0, u0 = self.generate(type, which=which, p=p)
            qs, rs, us = adjust_strides((q0, r0, u0))
            if p == 1:
                ai = np.insert(a, k, u0, 0 if which == 'row' else 1)
            else:
                ai = np.insert(a, k*np.ones(p, np.intp),
                        u0 if which == 'row' else u0, 
                        0 if which == 'row' else 1)

            # for each variable, q, r, u we try with it strided and
            # overwrite=False. Then we try with overwrite=True. Nothing
            # is checked to see if it can be overwritten, since only
            # F ordered Q can be overwritten when adding columns.

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            q1, r1 = qr_insert(qs, r, u, k, which, False)
            check_qr(q1, r1, ai, self.rtol, self.atol)
            q1o, r1o = qr_insert(qs, r, u, k, which, True)
            check_qr(q1o, r1o, ai, self.rtol, self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            q2, r2 = qr_insert(q, rs, u, k, which, False)
            check_qr(q2, r2, ai, self.rtol, self.atol)
            q2o, r2o = qr_insert(q, rs, u, k, which, True)
            check_qr(q2o, r2o, ai, self.rtol, self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            q3, r3 = qr_insert(q, r, us, k, which, False)
            check_qr(q3, r3, ai, self.rtol, self.atol)
            q3o, r3o = qr_insert(q, r, us, k, which, True)
            check_qr(q3o, r3o, ai, self.rtol, self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            # since some of these were consumed above
            qs, rs, us = adjust_strides((q, r, u))
            q5, r5 = qr_insert(qs, rs, us, k, which, False)
            check_qr(q5, r5, ai, self.rtol, self.atol)
            q5o, r5o = qr_insert(qs, rs, us, k, which, True)
            check_qr(q5o, r5o, ai, self.rtol, self.atol)
    
    def test_non_unit_strides_1_row(self):
        self.base_non_simple_strides(make_strided, 0, 1, 'row')

    def test_non_unit_strides_p_row(self):
        self.base_non_simple_strides(make_strided, 0, 3, 'row')

    def test_non_unit_strides_1_col(self):
        self.base_non_simple_strides(make_strided, 0, 1, 'col')

    def test_non_unit_strides_p_col(self):
        self.base_non_simple_strides(make_strided, 0, 3, 'col')

    def test_neg_strides_1_row(self):
        self.base_non_simple_strides(negate_strides, 0, 1, 'row')

    def test_neg_strides_p_row(self):
        self.base_non_simple_strides(negate_strides, 0, 3, 'row')

    def test_neg_strides_1_col(self):
        self.base_non_simple_strides(negate_strides, 0, 1, 'col')

    def test_neg_strides_p_col(self):
        self.base_non_simple_strides(negate_strides, 0, 3, 'col')

    def test_non_itemsize_strides_1_row(self):
        self.base_non_simple_strides(nonitemsize_strides, 0, 1, 'row')

    def test_non_itemsize_strides_p_row(self):
        self.base_non_simple_strides(nonitemsize_strides, 0, 3, 'row')

    def test_non_itemsize_strides_1_col(self):
        self.base_non_simple_strides(nonitemsize_strides, 0, 1, 'col')

    def test_non_itemsize_strides_p_col(self):
        self.base_non_simple_strides(nonitemsize_strides, 0, 3, 'col')

    def test_non_native_byte_order_1_row(self):
        self.base_non_simple_strides(make_nonnative, 0, 1, 'row')

    def test_non_native_byte_order_p_row(self):
        self.base_non_simple_strides(make_nonnative, 0, 3, 'row')

    def test_non_native_byte_order_1_col(self):
        self.base_non_simple_strides(make_nonnative, 0, 1, 'col')

    def test_non_native_byte_order_p_col(self):
        self.base_non_simple_strides(make_nonnative, 0, 3, 'col')

    def test_overwrite_qu_rank_1(self):
        # when inserting rows, the size of both Q and R change, so only
        # column inserts can overwrite q. Only complex column inserts 
        # with C ordered Q overwrite u. Any contiguous Q is overwritten
        # when inserting 1 column
        a, q0, r, u, = self.generate('sqr', which='col', p=1)
        q = q0.copy('C')
        u0 = u.copy()
        # don't overwrite
        q1, r1 = qr_insert(q, r, u, 0, 'col', False)
        a1 = np.insert(a, 0, u0, 1)
        check_qr(q1, r1, a1, self.rtol, self.atol)
        check_qr(q, r, a, self.rtol, self.atol)

        # try overwriting
        q2, r2 = qr_insert(q, r, u, 0, 'col', True)
        check_qr(q2, r2, a1, self.rtol, self.atol)
        # verify the overwriting
        assert_allclose(q2, q, rtol=self.rtol, atol=self.atol)
        assert_allclose(u, u0.conj(), self.rtol, self.atol)

        # now try with a fortran ordered Q
        qF = q0.copy('F')
        u1 = u0.copy()
        q3, r3 = qr_insert(qF, r, u1, 0, 'col', False)
        check_qr(q3, r3, a1, self.rtol, self.atol)
        check_qr(qF, r, a, self.rtol, self.atol)

        # try overwriting
        q4, r4 = qr_insert(qF, r, u1, 0, 'col', True)
        check_qr(q4, r4, a1, self.rtol, self.atol)
        assert_allclose(q4, qF, rtol=self.rtol, atol=self.atol)

    def test_overwrite_qu_rank_p(self):
        # when inserting rows, the size of both Q and R change, so only
        # column inserts can potentially overwrite Q.  In practice, only
        # F ordered Q are overwritten with a rank p update.
        a, q0, r, u, = self.generate('sqr', which='col', p=3)
        q = q0.copy('F')
        a1 = np.insert(a, np.zeros(3, np.intp), u, 1)

        # don't overwrite
        q1, r1 = qr_insert(q, r, u, 0, 'col', False)
        check_qr(q1, r1, a1, self.rtol, self.atol)
        check_qr(q, r, a, self.rtol, self.atol)

        # try overwriting
        q2, r2 = qr_insert(q, r, u, 0, 'col', True)
        check_qr(q2, r2, a1, self.rtol, self.atol)
        assert_allclose(q2, q, rtol=self.rtol, atol=self.atol)

    def test_empty_inputs(self):
        a, q, r, u = self.generate('sqr', which='row')
        assert_raises(ValueError, qr_insert, np.array([]), r, u, 0, 'row')
        assert_raises(ValueError, qr_insert, q, np.array([]), u, 0, 'row')
        assert_raises(ValueError, qr_insert, q, r, np.array([]), 0, 'row')
        assert_raises(ValueError, qr_insert, np.array([]), r, u, 0, 'col')
        assert_raises(ValueError, qr_insert, q, np.array([]), u, 0, 'col')
        assert_raises(ValueError, qr_insert, q, r, np.array([]), 0, 'col')

    def test_mismatched_shapes(self):
        a, q, r, u = self.generate('tall', which='row')
        assert_raises(ValueError, qr_insert, q, r[1:], u, 0, 'row')
        assert_raises(ValueError, qr_insert, q[:-2], r, u, 0, 'row')
        assert_raises(ValueError, qr_insert, q, r, u[1:], 0, 'row')
        assert_raises(ValueError, qr_insert, q, r[1:], u, 0, 'col')
        assert_raises(ValueError, qr_insert, q[:-2], r, u, 0, 'col')
        assert_raises(ValueError, qr_insert, q, r, u[1:], 0, 'col')

    def test_integer_input(self):
        q = np.arange(16).reshape(4, 4)
        r = q.copy()  # doesn't matter
        u = q[:, 0].copy()
        assert_raises(ValueError, qr_insert, q, r, u, 0, 'row')
        assert_raises(ValueError, qr_insert, q, r, u, 0, 'col')
    
    def test_check_finite(self):
        a0, q0, r0, u0 = self.generate('sqr', which='row', p=3)

        q = q0.copy('F')
        q[1,1] = np.nan
        assert_raises(ValueError, qr_insert, q, r0, u0[:,0], 0, 'row')
        assert_raises(ValueError, qr_insert, q, r0, u0, 0, 'row')
        assert_raises(ValueError, qr_insert, q, r0, u0[:,0], 0, 'col')
        assert_raises(ValueError, qr_insert, q, r0, u0, 0, 'col')

        r = r0.copy('F')
        r[1,1] = np.nan
        assert_raises(ValueError, qr_insert, q0, r, u0[:,0], 0, 'row')
        assert_raises(ValueError, qr_insert, q0, r, u0, 0, 'row')
        assert_raises(ValueError, qr_insert, q0, r, u0[:,0], 0, 'col')
        assert_raises(ValueError, qr_insert, q0, r, u0, 0, 'col')

        u = u0.copy('F')
        u[0,0] = np.nan
        assert_raises(ValueError, qr_insert, q0, r0, u[:,0], 0, 'row')
        assert_raises(ValueError, qr_insert, q0, r0, u, 0, 'row')
        assert_raises(ValueError, qr_insert, q0, r0, u[:,0], 0, 'col')
        assert_raises(ValueError, qr_insert, q0, r0, u, 0, 'col')

class TestQRinsert_f(BaseQRinsert):
    dtype = np.dtype('f')

class TestQRinsert_F(BaseQRinsert):
    dtype = np.dtype('F')

class TestQRinsert_d(BaseQRinsert):
    dtype = np.dtype('d')

class TestQRinsert_D(BaseQRinsert):
    dtype = np.dtype('D')

class BaseQRupdate(BaseQRdeltas):
    def generate(self, type, mode='full', p=1):
        a, q, r = super(BaseQRupdate, self).generate(type, mode)

        # super call set the seed...
        u = np.random.random((q.shape[0], p))
        v = np.random.random((r.shape[1], p))

        u = np.squeeze(u)
        v = np.squeeze(v)

        if np.iscomplexobj(self.dtype.type(1)):
            b = np.random.random(u.shape)
            u = u + 1j * b

            c = np.random.random(v.shape)
            v = v + 1j * c

        u = u.astype(self.dtype)
        v = v.astype(self.dtype)
        return a, q, r, u, v

    def test_sqr_rank_1(self):
        a, q, r, u, v = self.generate('sqr')
        q1, r1 = qr_update(q, r, u, v, False)
        a1 = a + np.outer(u, v.conj())
        check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_sqr_rank_p(self):
        # test ndim = 2, rank 1 updates here too
        for p in [1, 2, 3, 5]:
            a, q, r, u, v = self.generate('sqr', p=p)
            if p == 1:
                u = u.reshape(u.size, 1)
                v = v.reshape(v.size, 1)
            q1, r1 = qr_update(q, r, u, v, False)
            a1 = a + np.dot(u, v.T.conj())
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_rank_1(self):
        a, q, r, u, v = self.generate('tall')
        q1, r1 = qr_update(q, r, u, v, False)
        a1 = a + np.outer(u, v.conj())
        check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_tall_rank_p(self):
        for p in [1, 2, 3, 5]:
            a, q, r, u, v = self.generate('tall', p=p)
            if p == 1:
                u = u.reshape(u.size, 1)
                v = v.reshape(v.size, 1)
            q1, r1 = qr_update(q, r, u, v, False)
            a1 = a + np.dot(u, v.T.conj())
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_rank_1(self):
        a, q, r, u, v = self.generate('fat')
        q1, r1 = qr_update(q, r, u, v, False)
        a1 = a + np.outer(u, v.conj())
        check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_fat_rank_p(self):
        for p in [1, 2, 3, 5]:
            a, q, r, u, v = self.generate('fat', p=p)
            if p == 1:
                u = u.reshape(u.size, 1)
                v = v.reshape(v.size, 1)
            q1, r1 = qr_update(q, r, u, v, False)
            a1 = a + np.dot(u, v.T.conj())
            check_qr(q1, r1, a1, self.rtol, self.atol)

    def test_economic_rank_1(self):
        a, q, r, u, v = self.generate('tall', 'economic')
        q1, r1 = qr_update(q, r, u, v, False)
        a1 = a + np.outer(u, v.conj())
        check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def test_economic_rank_p(self):
        for p in [1, 2, 3, 5]:
            a, q, r, u, v = self.generate('tall', 'economic', p)
            if p == 1:
                u = u.reshape(u.size, 1)
                v = v.reshape(v.size, 1)
            q1, r1 = qr_update(q, r, u, v, False)
            a1 = a + np.dot(u, v.T.conj())
            check_qr(q1, r1, a1, self.rtol, self.atol, False)

    def base_non_simple_strides(self, adjust_strides, mode, p, overwriteable):
        assert_sqr = False if mode == 'economic' else True
        for type in ['sqr', 'tall', 'fat']:
            a, q0, r0, u0, v0 = self.generate(type, mode, p)
            qs, rs, us, vs = adjust_strides((q0, r0, u0, v0))
            if p == 1:
                aup = a + np.outer(u0, v0.conj())
            else:
                aup = a + np.dot(u0, v0.T.conj())

            # for each variable, q, r, u, v we try with it strided and
            # overwrite=False. Then we try with overwrite=True, and make
            # sure that if p == 1, r and v are still overwritten.
            # a strided q and u must always be copied.

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            v = v0.copy('C')
            q1, r1 = qr_update(qs, r, u, v, False)
            check_qr(q1, r1, aup, self.rtol, self.atol, assert_sqr)
            q1o, r1o = qr_update(qs, r, u, v, True)
            check_qr(q1o, r1o, aup, self.rtol, self.atol, assert_sqr)
            if overwriteable:
                assert_allclose(r1o, r, rtol=self.rtol, atol=self.atol)
                assert_allclose(v, v0.conj(), rtol=self.rtol, atol=self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            v = v0.copy('C')
            q2, r2 = qr_update(q, rs, u, v, False)
            check_qr(q2, r2, aup, self.rtol, self.atol, assert_sqr)
            q2o, r2o = qr_update(q, rs, u, v, True)
            check_qr(q2o, r2o, aup, self.rtol, self.atol, assert_sqr)
            if overwriteable:
                assert_allclose(r2o, rs, rtol=self.rtol, atol=self.atol)
                assert_allclose(v, v0.conj(), rtol=self.rtol, atol=self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            v = v0.copy('C')
            q3, r3 = qr_update(q, r, us, v, False)
            check_qr(q3, r3, aup, self.rtol, self.atol, assert_sqr)
            q3o, r3o = qr_update(q, r, us, v, True)
            check_qr(q3o, r3o, aup, self.rtol, self.atol, assert_sqr)
            if overwriteable:
                assert_allclose(r3o, r, rtol=self.rtol, atol=self.atol)
                assert_allclose(v, v0.conj(), rtol=self.rtol, atol=self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            v = v0.copy('C')
            q4, r4 = qr_update(q, r, u, vs, False)
            check_qr(q4, r4, aup, self.rtol, self.atol, assert_sqr)
            q4o, r4o = qr_update(q, r, u, vs, True)
            check_qr(q4o, r4o, aup, self.rtol, self.atol, assert_sqr)
            if overwriteable:
                assert_allclose(r4o, r, rtol=self.rtol, atol=self.atol)
                assert_allclose(vs, v0.conj(), rtol=self.rtol, atol=self.atol)

            q = q0.copy('F')
            r = r0.copy('F')
            u = u0.copy('F')
            v = v0.copy('C')
            # since some of these were consumed above
            qs, rs, us, vs = adjust_strides((q, r, u, v))
            q5, r5 = qr_update(qs, rs, us, vs, False)
            check_qr(q5, r5, aup, self.rtol, self.atol, assert_sqr)
            q5o, r5o = qr_update(qs, rs, us, vs, True)
            check_qr(q5o, r5o, aup, self.rtol, self.atol, assert_sqr)
            if overwriteable:
                assert_allclose(r5o, rs, rtol=self.rtol, atol=self.atol)
                assert_allclose(vs, v0.conj(), rtol=self.rtol, atol=self.atol)

    def test_non_unit_strides_rank_1(self):
        self.base_non_simple_strides(make_strided, 'full', 1, True)

    def test_non_unit_strides_economic_rank_1(self):
        self.base_non_simple_strides(make_strided, 'economic', 1, True)

    def test_non_unit_strides_rank_p(self):
        self.base_non_simple_strides(make_strided, 'full', 3, False)

    def test_non_unit_strides_economic_rank_p(self):
        self.base_non_simple_strides(make_strided, 'economic', 3, False)

    def test_neg_strides_rank_1(self):
        self.base_non_simple_strides(negate_strides, 'full', 1, False)

    def test_neg_strides_economic_rank_1(self):
        self.base_non_simple_strides(negate_strides, 'economic', 1, False)

    def test_neg_strides_rank_p(self):
        self.base_non_simple_strides(negate_strides, 'full', 3, False)

    def test_neg_strides_economic_rank_p(self):
        self.base_non_simple_strides(negate_strides, 'economic', 3, False)
     
    def test_non_itemsize_strides_rank_1(self):
        self.base_non_simple_strides(nonitemsize_strides, 'full', 1, False)

    def test_non_itemsize_strides_economic_rank_1(self):
        self.base_non_simple_strides(nonitemsize_strides, 'economic', 1, False)

    def test_non_itemsize_strides_rank_p(self):
        self.base_non_simple_strides(nonitemsize_strides, 'full', 3, False)

    def test_non_itemsize_strides_economic_rank_p(self):
        self.base_non_simple_strides(nonitemsize_strides, 'economic', 3, False)

    def test_non_native_byte_order_rank_1(self):
        self.base_non_simple_strides(make_nonnative, 'full', 1, False)

    def test_non_native_byte_order_economic_rank_1(self):
        self.base_non_simple_strides(make_nonnative, 'economic', 1, False)

    def test_non_native_byte_order_rank_p(self):
        self.base_non_simple_strides(make_nonnative, 'full', 3, False)

    def test_non_native_byte_order_economic_rank_p(self):
        self.base_non_simple_strides(make_nonnative, 'economic', 3, False)

    def test_overwrite_qruv_rank_1(self):
        # Any positive strided q, r, u, and v can be overwritten for a rank 1 
        # update, only checking C and F contiguous.
        a, q0, r0, u0, v0 = self.generate('sqr')
        a1 = a + np.outer(u0, v0.conj())
        q = q0.copy('F')
        r = r0.copy('F')
        u = u0.copy('F')
        v = v0.copy('F')

        # don't overwrite
        q1, r1 = qr_update(q, r, u, v, False)
        check_qr(q1, r1, a1, self.rtol, self.atol)
        check_qr(q, r, a, self.rtol, self.atol)

        q2, r2 = qr_update(q, r, u, v, True)
        check_qr(q2, r2, a1, self.rtol, self.atol)
        # verify the overwriting, no good way to check u and v.
        assert_allclose(q2, q, rtol=self.rtol, atol=self.atol)
        assert_allclose(r2, r, rtol=self.rtol, atol=self.atol)

        q = q0.copy('C')
        r = r0.copy('C')
        u = u0.copy('C')
        v = v0.copy('C')
        q3, r3 = qr_update(q, r, u, v, True)
        check_qr(q3, r3, a1, self.rtol, self.atol)
        assert_allclose(q3, q, rtol=self.rtol, atol=self.atol)
        assert_allclose(r3, r, rtol=self.rtol, atol=self.atol)

    def test_overwrite_qruv_rank_1_economic(self):
        # updating economic decompositions can overwrite any contigous r,
        # and positively strided r and u. V is only ever read.
        # only checking C and F contiguous.
        a, q0, r0, u0, v0 = self.generate('tall', 'economic')
        a1 = a + np.outer(u0, v0.conj())
        q = q0.copy('F')
        r = r0.copy('F')
        u = u0.copy('F')
        v = v0.copy('F')

        # don't overwrite
        q1, r1 = qr_update(q, r, u, v, False)
        check_qr(q1, r1, a1, self.rtol, self.atol, False)
        check_qr(q, r, a, self.rtol, self.atol, False)

        q2, r2 = qr_update(q, r, u, v, True)
        check_qr(q2, r2, a1, self.rtol, self.atol, False)
        # verify the overwriting, no good way to check u and v.
        assert_allclose(q2, q, rtol=self.rtol, atol=self.atol)
        assert_allclose(r2, r, rtol=self.rtol, atol=self.atol)

        q = q0.copy('C')
        r = r0.copy('C')
        u = u0.copy('C')
        v = v0.copy('C')
        q3, r3 = qr_update(q, r, u, v, True)
        check_qr(q3, r3, a1, self.rtol, self.atol, False)
        assert_allclose(q3, q, rtol=self.rtol, atol=self.atol)
        assert_allclose(r3, r, rtol=self.rtol, atol=self.atol)

    def test_overwrite_qruv_rank_p(self):
        # for rank p updates, q r must be F contiguous, v must be C (v.T --> F)
        # and u can be C or F, but is only overwritten if Q is C and complex
        a, q0, r0, u0, v0 = self.generate('sqr', p=3)
        a1 = a + np.dot(u0, v0.T.conj())
        q = q0.copy('F')
        r = r0.copy('F')
        u = u0.copy('F')
        v = v0.copy('C')

        # don't overwrite
        q1, r1 = qr_update(q, r, u, v, False)
        check_qr(q1, r1, a1, self.rtol, self.atol)
        check_qr(q, r, a, self.rtol, self.atol)

        q2, r2 = qr_update(q, r, u, v, True)
        check_qr(q2, r2, a1, self.rtol, self.atol)
        # verify the overwriting, no good way to check u and v.
        assert_allclose(q2, q, rtol=self.rtol, atol=self.atol)
        assert_allclose(r2, r, rtol=self.rtol, atol=self.atol)

    def test_empty_inputs(self):
        a, q, r, u, v = self.generate('tall')
        assert_raises(ValueError, qr_update, np.array([]), r, u, v)
        assert_raises(ValueError, qr_update, q, np.array([]), u, v)
        assert_raises(ValueError, qr_update, q, r, np.array([]), v)
        assert_raises(ValueError, qr_update, q, r, u, np.array([]))

    def test_mismatched_shapes(self):
        a, q, r, u, v = self.generate('tall')
        assert_raises(ValueError, qr_update, q, r[1:], u, v)
        assert_raises(ValueError, qr_update, q[:-2], r, u, v)
        assert_raises(ValueError, qr_update, q, r, u[1:], v)
        assert_raises(ValueError, qr_update, q, r, u, v[1:])

    def test_integer_input(self):
        q = np.arange(16).reshape(4, 4)
        r = q.copy()  # doesn't matter
        u = q[:, 0].copy()
        v = r[0, :].copy()
        assert_raises(ValueError, qr_update, q, r, u, v)

    def test_check_finite(self):
        a0, q0, r0, u0, v0 = self.generate('tall', p=3)

        q = q0.copy('F')
        q[1,1] = np.nan
        assert_raises(ValueError, qr_update, q, r0, u0[:,0], v0[:,0])
        assert_raises(ValueError, qr_update, q, r0, u0, v0)

        r = r0.copy('F')
        r[1,1] = np.nan
        assert_raises(ValueError, qr_update, q0, r, u0[:,0], v0[:,0])
        assert_raises(ValueError, qr_update, q0, r, u0, v0)

        u = u0.copy('F')
        u[0,0] = np.nan
        assert_raises(ValueError, qr_update, q0, r0, u[:,0], v0[:,0])
        assert_raises(ValueError, qr_update, q0, r0, u, v0)

        v = v0.copy('F')
        v[0,0] = np.nan
        assert_raises(ValueError, qr_update, q0, r0, u[:,0], v[:,0])
        assert_raises(ValueError, qr_update, q0, r0, u, v)

    def test_economic_check_finite(self):
        a0, q0, r0, u0, v0 = self.generate('tall', mode='economic', p=3)

        q = q0.copy('F')
        q[1,1] = np.nan
        assert_raises(ValueError, qr_update, q, r0, u0[:,0], v0[:,0])
        assert_raises(ValueError, qr_update, q, r0, u0, v0)

        r = r0.copy('F')
        r[1,1] = np.nan
        assert_raises(ValueError, qr_update, q0, r, u0[:,0], v0[:,0])
        assert_raises(ValueError, qr_update, q0, r, u0, v0)

        u = u0.copy('F')
        u[0,0] = np.nan
        assert_raises(ValueError, qr_update, q0, r0, u[:,0], v0[:,0])
        assert_raises(ValueError, qr_update, q0, r0, u, v0)

        v = v0.copy('F')
        v[0,0] = np.nan
        assert_raises(ValueError, qr_update, q0, r0, u[:,0], v[:,0])
        assert_raises(ValueError, qr_update, q0, r0, u, v)

class TestQRupdate_f(BaseQRupdate):
    dtype = np.dtype('f')

class TestQRupdate_F(BaseQRupdate):
    dtype = np.dtype('F')

class TestQRupdate_d(BaseQRupdate):
    dtype = np.dtype('d')

class TestQRupdate_D(BaseQRupdate):
    dtype = np.dtype('D')

def test_form_qTu():
    # we want to ensure that all of the code paths through this function are tested,
    # Most of them should be hit with the rest of test suite, but explicit tests make
    # clear precisely what is being tested.
    #
    # this function expects that Q is either C or F contiguous. But U have any positive strides.
    #
    # some of these test are duplicates, since contig 1d arrays are both C and F. 
    # should we test a Mx1 2d array?

    q_order = ['F', 'C']
    u_order = ['F', 'C', 'A']  # here A means is not F not C
    u_shape = [1, 3]
    dtype = ['f', 'd', 'F', 'D']

    for qo, uo, us, d in itertools.product(q_order, u_order, u_shape, dtype):
        yield check_form_qTu, qo, uo, us, d
    
def check_form_qTu(q_order, u_order, u_shape, dtype):
    np.random.seed(47)
    q_shape = (8,8)
    if u_shape == 1:
        u_shape = (8,)
    else:
        u_shape = (8, u_shape)
    dtype = np.dtype(dtype)

    if dtype.char in 'fd':
        q = np.random.random(q_shape)
        u = np.random.random(u_shape)
    elif dtype.char in 'FD':
        q = np.random.random(q_shape) + 1j*np.random.random(q_shape)
        u = np.random.random(u_shape) + 1j*np.random.random(u_shape)
    else:
        ValueError("form_qTu doesn't support this dtype")

    q = np.require(q, dtype, q_order)
    if u_order != 'A':
        u = np.require(u, dtype, u_order)
    else:
        u, = make_strided((u.astype(dtype),))
    
    rtol = 10.0 ** -(np.finfo(dtype).precision-2)
    atol = 2*np.finfo(dtype).eps

    expected = np.dot(q.T.conj(), u)
    res = _decomp_update._form_qTu(q, u)
    assert_allclose(res, expected, rtol=rtol, atol=atol)
 
