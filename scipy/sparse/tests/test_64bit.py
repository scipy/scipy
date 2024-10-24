""" Test functions involving 64bit or 32bit indexing """
import pytest
import numpy as np
from scipy.sparse import (
    bsr_array, coo_array, csc_array, csr_array, dia_array,
    bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dia_matrix,
)

# rename to avoid pytest collecting them in this module
from .test_base import (
    TestBSR as _TestBSR,
    TestCOO as _TestCOO,
    TestCSC as _TestCSC,
    TestCSR as _TestCSR,
    TestDIA as _TestDIA,
    TestDOK as _TestDOK,
    TestLIL as _TestLIL,
    TestBSRMatrix as _TestBSRMatrix,
    TestCOOMatrix as _TestCOOMatrix,
    TestCSCMatrix as _TestCSCMatrix,
    TestCSRMatrix as _TestCSRMatrix,
    TestDIAMatrix as _TestDIAMatrix,
    TestDOKMatrix as _TestDOKMatrix,
    TestLILMatrix as _TestLILMatrix,
    with_64bit_maxval_limit,
)

counter_container = [0]
zero_get_index_dtype_calls = set()

# name : reason not tested here
SKIP_TESTS = {
    'test_expm': 'expm for 64-bit indices not available',
    'test_inv': 'linsolve for 64-bit indices not available',
    'test_solve': 'linsolve for 64-bit indices not available',
    'test_scalar_idx_dtype': 'test implemented in base class',
    'test_large_dimensions_reshape': 'test actually requires 64-bit to work',
    'test_constructor_smallcol': 'test verifies int32 indexes',
    'test_constructor_largecol': 'test verifies int64 indexes',
    'test_tocoo_tocsr_tocsc_gh19245': 'test verifies int32 indexes',
}
# See the bottom of the module for IGNORE_TESTS


def cases_64bit(sp_api, ignore=True):
    """Yield all tests for all formats that use get_index_dtype

    This is more than testing get_index_dtype. It allows checking whether upcasting
    or downcasting the index dtypes affects test results. The approach used here
    does not try to figure out which tests might fail due to 32/64-bit issues.
    We just run them all.
    So, each test method in that uses cases_64bit reruns most of the test suite!
    """
    IGNORE = IGNORE_TESTS if ignore else {}
    if sp_api == "sparray":
        TEST_CLASSES = [_TestBSR, _TestCOO, _TestCSC, _TestCSR, _TestDIA]
    elif sp_api == "sparray-extra":
        # lil/dok->other conversion operations use get_index_dtype
        # so we include lil & dok test suite even though they do not
        # use get_index_dtype within the class. That means many of
        # these tests are superfluous, but it's hard to pick which
        TEST_CLASSES = [_TestDOK, _TestLIL]
    elif sp_api == "spmatrix":
        TEST_CLASSES = [_TestBSRMatrix, _TestCOOMatrix, _TestCSCMatrix,
                        _TestCSRMatrix, _TestDIAMatrix]
    elif sp_api == "spmatrix-extra":
        # lil/dok->other conversion operations use get_index_dtype
        TEST_CLASSES = [_TestDOKMatrix, _TestLILMatrix]
    else:
        raise ValueError(f"parameter {sp_api=} is not valid")

    for cls in TEST_CLASSES:
        for method_name in sorted(dir(cls)):
            name = f'{cls.__name__}-{method_name}'
            method = getattr(cls, method_name)
            if (method_name.startswith('test_') and
                    not getattr(method, 'slow', False)):
                marks = []

                if name in IGNORE:
                    continue

                msg = SKIP_TESTS.get(method_name)
                if msg:
                    marks.append(pytest.mark.skip(reason=msg))

                markers = getattr(method, 'pytestmark', [])
                for mark in markers:
                    if mark.name in ('skipif', 'skip', 'xfail', 'xslow'):
                        marks.append(mark)

                yield pytest.param(cls, method_name, marks=marks)


class RunAll64Bit:
    def _check_resiliency(self, cls, method_name, **kw):
        # Resiliency test, to check that sparse matrices deal reasonably
        # with varying index data types.

        @with_64bit_maxval_limit(**kw)
        def check(cls, method_name):
            instance = cls()
            if hasattr(instance, 'setup_method'):
                instance.setup_method()
            try:
                getattr(instance, method_name)()
            finally:
                if hasattr(instance, 'teardown_method'):
                    instance.teardown_method()

        check(cls, method_name)

    def _check_resiliency_and_list_of_ignores(self, cls, method_name):
        @with_64bit_maxval_limit(random=counter_container)
        def check_and_collect(cls, method_name):
            instance = cls()
            if hasattr(instance, 'setup_method'):
                instance.setup_method()
            try:
                counter_container[0] = 0
                getattr(instance, method_name)()
                if counter_container[0] == 0:
                    zero_calls.add(name)
            finally:
                if hasattr(instance, 'teardown_method'):
                    instance.teardown_method()

        name = f'{cls.__name__}-{method_name}'
        zero_calls = zero_get_index_dtype_calls

        check_and_collect(cls, method_name)

        if name in zero_calls:
            assert name in IGNORE_TESTS


class Test64BitArray(RunAll64Bit):
    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray", ignore=False))
    def test_resiliency_random(self, cls, method_name):
        self._check_resiliency_and_list_of_ignores(cls, method_name)


class Test64BitMatrix(RunAll64Bit):
    # assert_32bit=True only for spmatrix cuz sparray does not check index content
    @pytest.mark.fail_slow(5)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_no_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, assert_32bit=True)


#@pytest.mark.slow
class Test64BitMatrixSameAsArray(RunAll64Bit):
    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix", ignore=False))
    def test_resiliency_random(self, cls, method_name):
        # Resiliency check that sparse deals with varying index data types.
        self._check_resiliency_and_list_of_ignores(cls, method_name)


@pytest.mark.xslow
class Test64BitArrayExtra(RunAll64Bit):
    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra", False))
    def test_resiliency_random(self, cls, method_name):
        # Resiliency check that sparse deals with varying index data types.
        self._check_resiliency_and_list_of_ignores(cls, method_name)


@pytest.mark.xslow
class Test64BitMatrixExtra(RunAll64Bit):
    # assert_32bit=True only for spmatrix cuz sparray does not check index content
    @pytest.mark.fail_slow(5)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_no_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, assert_32bit=True)

    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra", False))
    def test_resiliency_random(self, cls, method_name):
        # Resiliency check that sparse deals with varying index data types.
        self._check_resiliency_and_list_of_ignores(cls, method_name)


@pytest.mark.xslow
def test_report_zero_get_index_dtype_calls():
    print(zero_get_index_dtype_calls)
    print('IGNORE_TESTS = {')
    for name in sorted(zero_get_index_dtype_calls):
        print(f"    '{name}': 'get_index_dtype not called',")
    print('}')

    # these < unless you run entire module with xslow and slow. Then ==
    # Test order matters here. Must be after all Test64Bit___.test_resiliency_random
    # Usually put last in the order (see conftest.py)
    assert zero_get_index_dtype_calls == IGNORE_TESTS.keys(), 'Found New Tests'


class Test64BitTools:
    # classes that use get_index_dtype
    MAT_CLASSES = [
        bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dia_matrix,
        bsr_array, coo_array, csc_array, csr_array, dia_array,
    ]

    def _compare_index_dtype(self, m, dtype):
        dtype = np.dtype(dtype)
        if m.format in ['csc', 'csr', 'bsr']:
            return (m.indices.dtype == dtype) and (m.indptr.dtype == dtype)
        elif m.format == 'coo':
            return (m.row.dtype == dtype) and (m.col.dtype == dtype)
        elif m.format == 'dia':
            return (m.offsets.dtype == dtype)
        else:
            raise ValueError(f"matrix {m!r} has no integer indices")

    def test_decorator_maxval_limit(self):
        # Test that the with_64bit_maxval_limit decorator works

        @with_64bit_maxval_limit(maxval_limit=10)
        def check(mat_cls):
            m = mat_cls(np.random.rand(10, 1))
            assert self._compare_index_dtype(m, np.int32)
            m = mat_cls(np.random.rand(11, 1))
            assert self._compare_index_dtype(m, np.int64)

        for mat_cls in self.MAT_CLASSES:
            check(mat_cls)

    def test_decorator_maxval_random(self):
        # Test that the with_64bit_maxval_limit decorator works (2)

        @with_64bit_maxval_limit(random=True)
        def check(mat_cls):
            seen_32 = False
            seen_64 = False
            for k in range(100):
                m = mat_cls(np.random.rand(9, 9))
                seen_32 = seen_32 or self._compare_index_dtype(m, np.int32)
                seen_64 = seen_64 or self._compare_index_dtype(m, np.int64)
                if seen_32 and seen_64:
                    break
            else:
                raise AssertionError("both 32 and 64 bit indices not seen")

        for mat_cls in self.MAT_CLASSES:
            check(mat_cls)

    def test_downcast_intp(self):
        # Check that bincount and ufunc.reduceat intp downcasts are
        # dealt with. The point here is to trigger points in the code
        # that can fail on 32-bit systems when using 64-bit indices,
        # due to use of functions that only work with intp-size indices.

        @with_64bit_maxval_limit(fixed_dtype=np.int64, downcast_maxval=1)
        def check_limited(csc_container, csr_container, coo_container):
            # These involve indices larger than `downcast_maxval`
            a = csc_container([[1, 2], [3, 4], [5, 6]])
            pytest.raises(AssertionError, a.count_nonzero, axis=1)
            pytest.raises(AssertionError, a.sum, axis=0)

            a = csr_container([[1, 2, 3], [3, 4, 6]])
            pytest.raises(AssertionError, a.count_nonzero, axis=0)
            pytest.raises(AssertionError, a.sum, axis=1)

            a = coo_container([[1, 2, 3], [3, 4, 5]])
            pytest.raises(AssertionError, a.count_nonzero, axis=0)
            a.has_canonical_format = False
            pytest.raises(AssertionError, a.sum_duplicates)

        @with_64bit_maxval_limit(fixed_dtype=np.int64)
        def check_unlimited(csc_container, csr_container, coo_container):
            # These involve indices smaller than `downcast_maxval`
            a = csc_container([[1, 2], [3, 4], [5, 6]])
            a.count_nonzero(axis=1)
            a.sum(axis=0)

            a = csr_container([[1, 2, 3], [3, 4, 6]])
            a.count_nonzero(axis=0)
            a.sum(axis=1)

            a = coo_container([[1, 2, 3], [3, 4, 5]])
            a.count_nonzero(axis=0)
            a.has_canonical_format = False
            a.sum_duplicates()

        check_limited(csc_array, csr_array, coo_array)
        check_unlimited(csc_array, csr_array, coo_array)
        check_limited(csc_matrix, csr_matrix, coo_matrix)
        check_unlimited(csc_matrix, csr_matrix, coo_matrix)


# We could track includes rather than ignores. But this way no new tests get skipped
IGNORE_TESTS = {
#}
#BYPASS = {
    'TestBSR-test_invalid_shapes': 'get_index_dtype not called',
    'TestBSR-test_rmul_scalar_type_error': 'get_index_dtype not called',
    'TestBSRMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestCOO-test_add_dense': 'get_index_dtype not called',
    'TestCOO-test_invalid_shapes': 'get_index_dtype not called',
    'TestCOO-test_radd': 'get_index_dtype not called',
    'TestCOOMatrix-test_add_dense': 'get_index_dtype not called',
    'TestCOOMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestCOOMatrix-test_radd': 'get_index_dtype not called',
    'TestCSC-test_add_dense': 'get_index_dtype not called',
    'TestCSC-test_invalid_shapes': 'get_index_dtype not called',
    'TestCSC-test_radd': 'get_index_dtype not called',
    'TestCSC-test_rmul_scalar_type_error': 'get_index_dtype not called',
    'TestCSCMatrix-test_add_dense': 'get_index_dtype not called',
    'TestCSCMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestCSCMatrix-test_radd': 'get_index_dtype not called',
    'TestCSR-test_add_dense': 'get_index_dtype not called',
    'TestCSR-test_invalid_shapes': 'get_index_dtype not called',
    'TestCSR-test_radd': 'get_index_dtype not called',
    'TestCSR-test_rmul_scalar_type_error': 'get_index_dtype not called',
    'TestCSRMatrix-test_add_dense': 'get_index_dtype not called',
    'TestCSRMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestCSRMatrix-test_radd': 'get_index_dtype not called',
    'TestDIA-test_invalid_shapes': 'get_index_dtype not called',
    'TestDIA-test_setdiag_dtype': 'get_index_dtype not called',
    'TestDIAMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestDIAMatrix-test_setdiag_dtype': 'get_index_dtype not called',
    'TestDOK-test_dtype_preservation_empty_index': 'get_index_dtype not called',
    'TestDOK-test_fancy_assignment_dtypes': 'get_index_dtype not called',
    'TestDOK-test_invalid_shapes': 'get_index_dtype not called',
    'TestDOK-test_negative_index_assignment': 'get_index_dtype not called',
    'TestDOK-test_scalar_assign_2': 'get_index_dtype not called',
    'TestDOK-test_ticket1160': 'get_index_dtype not called',
    'TestDOKMatrix-test_dtype_preservation': 'get_index_dtype not called',
    'TestDOKMatrix-test_dtype_preservation_empty_index': 'get_index_dtype not called',
    'TestDOKMatrix-test_dtype_preservation_empty_slice': 'get_index_dtype not called',
    'TestDOKMatrix-test_fancy_assignment_dtypes': 'get_index_dtype not called',
    'TestDOKMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestDOKMatrix-test_negative_index_assignment': 'get_index_dtype not called',
    'TestDOKMatrix-test_rmul_scalar_type_error': 'get_index_dtype not called',
    'TestDOKMatrix-test_scalar_assign_2': 'get_index_dtype not called',
    'TestDOKMatrix-test_ticket1160': 'get_index_dtype not called',
    'TestLIL-test_copy': 'get_index_dtype not called',
    'TestLIL-test_dtype_preservation_empty_index': 'get_index_dtype not called',
    'TestLIL-test_empty': 'get_index_dtype not called',
    'TestLIL-test_fancy_indexing_multidim_set': 'get_index_dtype not called',
    'TestLIL-test_fancy_indexing_set': 'get_index_dtype not called',
    'TestLIL-test_idiv_scalar': 'get_index_dtype not called',
    'TestLIL-test_imul_scalar': 'get_index_dtype not called',
    'TestLIL-test_index_scalar_assign': 'get_index_dtype not called',
    'TestLIL-test_invalid_shapes': 'get_index_dtype not called',
    'TestLIL-test_negative_index_assignment': 'get_index_dtype not called',
    'TestLIL-test_pickle': 'get_index_dtype not called',
    'TestLIL-test_scalar_mul': 'get_index_dtype not called',
    'TestLIL-test_setdiag_comprehensive': 'get_index_dtype not called',
    'TestLIL-test_setelement': 'get_index_dtype not called',
    'TestLIL-test_slice_assign_2': 'get_index_dtype not called',
    'TestLIL-test_slice_scalar_assign': 'get_index_dtype not called',
    'TestLIL-test_truediv_scalar': 'get_index_dtype not called',
    'TestLILMatrix-test_copy': 'get_index_dtype not called',
    'TestLILMatrix-test_dtype_preservation': 'get_index_dtype not called',
    'TestLILMatrix-test_dtype_preservation_empty_index': 'get_index_dtype not called',
    'TestLILMatrix-test_dtype_preservation_empty_slice': 'get_index_dtype not called',
    'TestLILMatrix-test_empty': 'get_index_dtype not called',
    'TestLILMatrix-test_fancy_indexing_multidim_set': 'get_index_dtype not called',
    'TestLILMatrix-test_fancy_indexing_set': 'get_index_dtype not called',
    'TestLILMatrix-test_idiv_scalar': 'get_index_dtype not called',
    'TestLILMatrix-test_imul_scalar': 'get_index_dtype not called',
    'TestLILMatrix-test_index_scalar_assign': 'get_index_dtype not called',
    'TestLILMatrix-test_invalid_shapes': 'get_index_dtype not called',
    'TestLILMatrix-test_negative_index_assignment': 'get_index_dtype not called',
    'TestLILMatrix-test_pickle': 'get_index_dtype not called',
    'TestLILMatrix-test_scalar_mul': 'get_index_dtype not called',
    'TestLILMatrix-test_setdiag_comprehensive': 'get_index_dtype not called',
    'TestLILMatrix-test_setelement': 'get_index_dtype not called',
    'TestLILMatrix-test_slice_assign_2': 'get_index_dtype not called',
    'TestLILMatrix-test_slice_scalar_assign': 'get_index_dtype not called',
    'TestLILMatrix-test_truediv_scalar': 'get_index_dtype not called',
}
