from scipy.sparse import (
    bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dia_matrix, lil_matrix, dok_matrix,
    bsr_array, coo_array, csc_array, csr_array, dia_array, lil_array, dok_array,
    is_bsr, is_coo, is_csc, is_csr, is_dia, is_lil, is_dok,
    is_bsr_matrix, is_coo_matrix, is_csc_matrix, is_csr_matrix,
    is_dia_matrix, is_lil_matrix, is_dok_matrix,
    is_bsr_array, is_coo_array, is_csc_array, is_csr_array,
    is_dia_array, is_lil_array, is_dok_array,
    issparse, isspmatrix, issparray,
)


MAT_CLASSES = [bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dia_matrix, lil_matrix, dok_matrix]
ARR_CLASSES = [bsr_array, coo_array, csc_array, csr_array, dia_array, lil_array, dok_array]
IS_CLASSES = [is_bsr, is_coo, is_csc, is_csr, is_dia, is_lil, is_dok]
IS_MAT_CLASSES = [is_bsr_matrix, is_coo_matrix, is_csc_matrix, is_csr_matrix,
                  is_dia_matrix, is_lil_matrix, is_dok_matrix]
IS_ARR_CLASSES = [is_bsr_array, is_coo_array, is_csc_array, is_csr_array,
                  is_dia_array, is_lil_array, is_dok_array]


def test_ismatrix():
    for cls, is_xxx, is_xxx_matrix in zip(MAT_CLASSES, IS_CLASSES, IS_MAT_CLASSES):
        M = cls((2, 3))
        assert issparse(M)
        #assert isspmatrix(M)
        assert not issparray(M)
        assert is_xxx(M)
        assert is_xxx_matrix(M)
        for is_, is_m, is_a in zip(IS_CLASSES, IS_MAT_CLASSES, IS_ARR_CLASSES):
            if is_ != is_xxx:
                assert not is_(M)
                assert not is_m(M)
            assert not is_a(M)


def test_isarray():
    for cls, is_xxx, is_xxx_array in zip(ARR_CLASSES, IS_CLASSES, IS_ARR_CLASSES):
        A = cls((2, 3))
        assert issparse(A)
        #assert not isspmatrix(A)
        assert issparray(A)
        assert is_xxx(A)
        assert is_xxx_array(A)
        for is_, is_m, is_a in zip(IS_CLASSES, IS_MAT_CLASSES, IS_ARR_CLASSES):
            if is_ != is_xxx:
                assert not is_(A)
                assert not is_a(A)
            assert not is_m(A)
