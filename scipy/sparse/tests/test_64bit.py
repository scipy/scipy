""" Test functions involving 64bit or 32bit indexing """
import pytest
import numpy as np
import scipy
from scipy._lib.decorator import decorator
from scipy.sparse._sputils import get_index_dtype

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


#Todo: Revisit 64bit tests: avoid rerun of all tests for each version of get_index_dtype
def cases_64bit(sp_api):
    """Yield all tests for all formats that use get_index_dtype

    This is more than testing get_index_dtype. It allows checking whether upcasting
    or downcasting the index dtypes affects test results. The approach used here
    does not try to figure out which tests might fail due to 32/64-bit issues.
    We just run them all.
    So, each test method in that uses cases_64bit reruns most of the test suite!
    """
    IGNORE = IGNORE_TESTS
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
    elif sp_api == "all":
        TEST_CLASSES = [_TestBSR, _TestCOO, _TestCSC, _TestCSR, _TestDIA, _TestDOK, _TestLIL,
                        _TestBSRMatrix, _TestCOOMatrix, _TestCSCMatrix, _TestCSRMatrix,
                        _TestDIAMatrix, _TestDOKMatrix, _TestLILMatrix]
        IGNORE = {}
    else:
        raise ValueError(f"parameter {sp_api=} is not valid")

    for cls in TEST_CLASSES:
        for method_name in sorted(dir(cls)):
            name = f'{cls.__name__}-{method_name}'
            method = getattr(cls, method_name)
            if (method_name.startswith('test_') and
                    not getattr(method, 'slow', False)):
                marks = []

                if IGNORE.get(name):
                    continue

                msg = SKIP_TESTS.get(method_name)
                if msg:
                    marks.append(pytest.mark.skip(reason=msg))

                markers = getattr(method, 'pytestmark', [])
                for mark in markers:
                    if mark.name in ('skipif', 'skip', 'xfail', 'xslow'):
                        marks.append(mark)

                yield pytest.param(cls, method_name, marks=marks)


# Testing both spmatrices and sparrays for 64bit index dtype handling is
# expensive and double-checks the same code (e.g. _coobase)
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

    @pytest.mark.xslow
    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_random(self, cls, method_name):
        # bsr_array.eliminate_zeros relies on csr_array constructor
        # not making copies of index arrays --- this is not
        # necessarily true when we pick the index data type randomly
        self._check_resiliency(cls, method_name, random=True)


class Test64BitMatrix(RunAll64Bit):
    # assert_32bit=True only for spmatrix cuz sparray does not check index content
    @pytest.mark.fail_slow(5)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_no_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, assert_32bit=True)


@pytest.mark.slow
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
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_random(self, cls, method_name):
        # bsr_array.eliminate_zeros relies on csr_array constructor
        # not making copies of index arrays --- this is not
        # necessarily true when we pick the index data type randomly
        self._check_resiliency(cls, method_name, random=True)


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
   @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
   def test_resiliency_random(self, cls, method_name):
       # bsr_array.eliminate_zeros relies on csr_array constructor
       # not making copies of index arrays --- this is not
       # necessarily true when we pick the index data type randomly
       self._check_resiliency(cls, method_name, random=True)


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
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_random(self, cls, method_name):
        # bsr_array.eliminate_zeros relies on csr_array constructor
        # not making copies of index arrays --- this is not
        # necessarily true when we pick the index data type randomly
        self._check_resiliency(cls, method_name, random=True)


# Start tests of the accuracy of IGNORE_TESTS
# It should include all test methods that dont call get_index_dtype

counter_container = [0]
zero_get_index_dtype_calls = set()

# mimic with_64bit_maxval_limit but needs access to counter_container
def with_index_dtype_tracking():
    def new_get_index_dtype(arrays=(), maxval=None, check_contents=False):
        counter_container[0] = counter_container[0] + 1
        return get_index_dtype(arrays, maxval, check_contents)

    @decorator
    def deco(func, *a, **kw):
        backup = []
        modules = [scipy.sparse._bsr, scipy.sparse._coo, scipy.sparse._csc,
                   scipy.sparse._csr, scipy.sparse._dia, scipy.sparse._dok,
                   scipy.sparse._lil, scipy.sparse._sputils,
                   scipy.sparse._compressed, scipy.sparse._construct]
        try:
            for mod in modules:
                backup.append((mod, 'get_index_dtype',
                               getattr(mod, 'get_index_dtype', None)))
                setattr(mod, 'get_index_dtype', new_get_index_dtype)
            return func(*a, **kw)
        finally:
            for mod, name, oldfunc in backup:
                if oldfunc is not None:
                    setattr(mod, name, oldfunc)

    return deco

@pytest.mark.xslow
@pytest.mark.parametrize('cls,method_name', cases_64bit("all"))
def test_methods_to_ignore(cls, method_name):

    @with_index_dtype_tracking()
    def collect(cls, method_name):
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

    # move to local namespace
    zero_calls = zero_get_index_dtype_calls

    name = f'{cls.__name__}-{method_name}'
    collect(cls, method_name)
    print(name, counter_container)

    if name in zero_calls:
        assert name in IGNORE_TESTS

@pytest.mark.xslow
def test_report_zero_get_index_dtype_calls():
    print('IGNORE_TESTS = {')
    for name in sorted(zero_get_index_dtype_calls):
        print(f"    '{name}': 'get_index_dtype not called',")
    print('}')
    assert zero_get_index_dtype_calls == IGNORE_TESTS.keys(), 'Found New Tests'

# The following features are missing, so skip the tests:
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
# Of the 1069 test methods in TestCSR, 722 do not create any get_index_dtype calls
# We could list includes rather than skip. But this way no new tests get skipped
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

##################################3
#IGNORE_TESTS = {
#    'test_ellipsis_fancy_slicing': "get_index_dtype not called",
#    'test_unary_ufunc_overrides': "get_index_dtype not called",
#    'test_truediv_scalar': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_slicing_2': "get_index_dtype not called",
#    'test_slice_scalar_assign': "get_index_dtype not called",
#    'test_slice_assign_2': "get_index_dtype not called",
#    'test_setelement': "get_index_dtype not called",
#    'test_setdiag_comprehensive': "get_index_dtype not called",
#    'test_set_slice': "get_index_dtype not called",
#    'test_scalar_mul': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_ne': "get_index_dtype not called",
#    'test_mul_scalar': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_lt': "get_index_dtype not called",
#    'test_lil_multiply_removal': "get_index_dtype not called",
#    'test_lil_iteration': "get_index_dtype not called",
#    'test_le': "get_index_dtype not called",
#    'test_iterator': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_index_scalar_assign': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_gt': "get_index_dtype not called",
#    'test_getnnz_axis': "get_index_dtype not called",
#    'test_get_slices': "get_index_dtype not called",
#    'test_ge': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_set': "get_index_dtype not called",
#    'test_fancy_indexing_seq_assign': "get_index_dtype not called",
#    'test_fancy_indexing_multidim_set': "get_index_dtype not called",
#    'test_fancy_indexing_boolean': "get_index_dtype not called",
#    'test_fancy_indexing': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_fancy_assign_empty': "get_index_dtype not called",
#    'test_eq': "get_index_dtype not called",
#    'test_empty': "get_index_dtype not called",
#    'test_ellipsis_slicing': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dtype_preservation_empty_index': "get_index_dtype not called",
#    'test_dtype_preservation': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_assign_1d_slice': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_add0': "get_index_dtype not called",
#    'test_add': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#    'test_unary_ufunc_overrides': "get_index_dtype not called",
#    'test_trace': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_tobsr': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_ticket1160': "get_index_dtype not called",
#    'test_sum_out': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_sum': "get_index_dtype not called",
#    'test_sub_dense': "get_index_dtype not called",
#    'test_str_maxprint': "get_index_dtype not called",
#    'test_str': "get_index_dtype not called",
#    'test_star_vs_at_sign_for_sparray_and_spmatrix': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_size_zero_conversions': "get_index_dtype not called",
#    'test_scalar_assign_2': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_round': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_resize': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_real': "get_index_dtype not called",
#    'test_radd': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_nonzero': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_ne': "get_index_dtype not called",
#    'test_nanminmax': "get_index_dtype not called",
#    'test_multiple_ellipsis_slicing': "get_index_dtype not called",
#    'test_mul_scalar': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_minmax_axis': "get_index_dtype not called",
#    'test_minmax': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmul': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_lt': "get_index_dtype not called",
#    'test_le': "get_index_dtype not called",
#    'test_iterator': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_inplace_success': "get_index_dtype not called",
#    'test_inplace_dense': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_gt': "get_index_dtype not called",
#    'test_getrow': "get_index_dtype not called",
#    'test_getnnz_axis': "get_index_dtype not called",
#    'test_getelement': "get_index_dtype not called",
#    'test_getcol': "get_index_dtype not called",
#    'test_get_horiz_slice': "get_index_dtype not called",
#    'test_ge': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_seq_assign': "get_index_dtype not called",
#    'test_fancy_indexing_ndarray': "get_index_dtype not called",
#    'test_fancy_indexing_empty': "get_index_dtype not called",
#    'test_fancy_assignment_dtypes': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_eq': "get_index_dtype not called",
#    'test_empty_arithmetic': "get_index_dtype not called",
#    'test_empty': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dtype_preservation_empty_index': "get_index_dtype not called",
#    'test_dtype_preservation': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bool_rollover': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_asfptype': "get_index_dtype not called",
#    'test_argmax': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_add0': "get_index_dtype not called",
#    'test_add': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#    'test_trace': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_sum_out': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_sum': "get_index_dtype not called",
#    'test_sub_dense': "get_index_dtype not called",
#    'test_str': "get_index_dtype not called",
#    'test_sparsity_modifying_assignment': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_slicing_3': "get_index_dtype not called",
#    'test_slicing_2': "get_index_dtype not called",
#    'test_slice_scalar_assign': "get_index_dtype not called",
#    'test_slice_assignment': "get_index_dtype not called",
#    'test_slice_assign_2': "get_index_dtype not called",
#    'test_size_zero_conversions': "get_index_dtype not called",
#    'test_setelement': "get_index_dtype not called",
#    'test_setdiag_dtype': "get_index_dtype not called",
#    'test_set_slice': "get_index_dtype not called",
#    'test_sequence_assignment': "get_index_dtype not called",
#    'test_self_self_assignment': "get_index_dtype not called",
#    'test_scalar_assign_2': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_round': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_radd': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_non_unit_stride_2d_indexing': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_ne': "get_index_dtype not called",
#    'test_nanminmax': "get_index_dtype not called",
#    'test_multiple_ellipsis_slicing': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_minmax_axis': "get_index_dtype not called",
#    'test_minmax': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_lt': "get_index_dtype not called",
#    'test_le': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_inplace_success': "get_index_dtype not called",
#    'test_index_scalar_assign': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_gt': "get_index_dtype not called",
#    'test_getrow': "get_index_dtype not called",
#    'test_getelement': "get_index_dtype not called",
#    'test_getcol': "get_index_dtype not called",
#    'test_get_vert_slice': "get_index_dtype not called",
#    'test_get_slices': "get_index_dtype not called",
#    'test_get_horiz_slice': "get_index_dtype not called",
#    'test_ge': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_sparse_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_set': "get_index_dtype not called",
#    'test_fancy_indexing_seq_assign': "get_index_dtype not called",
#    'test_fancy_indexing_regression_3087': "get_index_dtype not called",
#    'test_fancy_indexing_randomized': "get_index_dtype not called",
#    'test_fancy_indexing_ndarray': "get_index_dtype not called",
#    'test_fancy_indexing_multidim_set': "get_index_dtype not called",
#    'test_fancy_indexing_empty': "get_index_dtype not called",
#    'test_fancy_indexing_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_2d_assign': "get_index_dtype not called",
#    'test_fancy_indexing': "get_index_dtype not called",
#    'test_fancy_assignment_dtypes': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_fancy_assign_empty': "get_index_dtype not called",
#    'test_eq': "get_index_dtype not called",
#    'test_empty_arithmetic': "get_index_dtype not called",
#    'test_ellipsis_slicing': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dtype_preservation_empty_index': "get_index_dtype not called",
#    'test_dtype_preservation': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_diagonal': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_assign_empty': "get_index_dtype not called",
#    'test_assign_1d_slice': "get_index_dtype not called",
#    'test_asfptype': "get_index_dtype not called",
#    'test_argmax': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_add0': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#    'test_ufuncs': "get_index_dtype not called",
#    'test_transpose': "get_index_dtype not called",
#    'test_trace': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_tobsr': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_sum_out': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_sum': "get_index_dtype not called",
#    'test_sub_dense': "get_index_dtype not called",
#    'test_str_maxprint': "get_index_dtype not called",
#    'test_str': "get_index_dtype not called",
#    'test_star_vs_at_sign_for_sparray_and_spmatrix': "get_index_dtype not called",
#    'test_sparse_format_conversions': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_size_zero_matrix_arithmetic': "get_index_dtype not called",
#    'test_size_zero_conversions': "get_index_dtype not called",
#    'test_scalar_assign_2': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_round': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_rmatvec': "get_index_dtype not called",
#    'test_resize': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_real': "get_index_dtype not called",
#    'test_radd': "get_index_dtype not called",
#    'test_pow': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_nonzero': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_neg': "get_index_dtype not called",
#    'test_multiple_ellipsis_slicing': "get_index_dtype not called",
#    'test_mul_scalar': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_mu': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmul': "get_index_dtype not called",
#    'test_matmat_sparse': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_iterator': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_inplace_dense': "get_index_dtype not called",
#    'test_index_scalar_assign': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_getnnz_axis': "get_index_dtype not called",
#    'test_getelement': "get_index_dtype not called",
#    'test_getcol': "get_index_dtype not called",
#    'test_get_horiz_slice': "get_index_dtype not called",
#    'test_from_sparse': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_sparse_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_seq_assign': "get_index_dtype not called",
#    'test_fancy_indexing_randomized': "get_index_dtype not called",
#    'test_fancy_indexing_broadcast': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_fancy_assign_empty': "get_index_dtype not called",
#    'test_empty_arithmetic': "get_index_dtype not called",
#    'test_empty': "get_index_dtype not called",
#    'test_elementwise_power': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_diagonal': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_constructor3': "get_index_dtype not called",
#    'test_constructor1_base': "get_index_dtype not called",
#    'test_constructor1': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bool_rollover': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_assign_empty': "get_index_dtype not called",
#    'test_assign_1d_slice': "get_index_dtype not called",
#    'test_asfptype': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_add': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#    'test_transpose': "get_index_dtype not called",
#    'test_trace': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_sum_out': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_sum': "get_index_dtype not called",
#    'test_str': "get_index_dtype not called",
#    'test_star_vs_at_sign_for_sparray_and_spmatrix': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_size_zero_conversions': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_round': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_rmatvec': "get_index_dtype not called",
#    'test_resize': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_real': "get_index_dtype not called",
#    'test_radd': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_nonzero': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_multiple_ellipsis_slicing': "get_index_dtype not called",
#    'test_mul_scalar': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_minmax': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_inplace_success': "get_index_dtype not called",
#    'test_inplace_dense': "get_index_dtype not called",
#    'test_index_scalar_assign': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_getnnz_axis': "get_index_dtype not called",
#    'test_getelement': "get_index_dtype not called",
#    'test_getcol': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_sparse_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_seq_assign': "get_index_dtype not called",
#    'test_fancy_indexing_randomized': "get_index_dtype not called",
#    'test_fancy_indexing_broadcast': "get_index_dtype not called",
#    'test_fancy_assignment_dtypes': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_fancy_assign_empty': "get_index_dtype not called",
#    'test_empty_arithmetic': "get_index_dtype not called",
#    'test_empty': "get_index_dtype not called",
#    'test_ellipsis_slicing': "get_index_dtype not called",
#    'test_eliminate_zeros': "get_index_dtype not called",
#    'test_elementwise_power': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dtype_preservation': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_constructor6': "get_index_dtype not called",
#    'test_constructor5': "get_index_dtype not called",
#    'test_constructor4': "get_index_dtype not called",
#    'test_constructor1': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_assign_empty': "get_index_dtype not called",
#    'test_assign_1d_slice': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#    'test_unary_ufunc_overrides': "get_index_dtype not called",
#    'test_transpose': "get_index_dtype not called",
#    'test_trace': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_tobsr': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_sum_out': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_sub_dense': "get_index_dtype not called",
#    'test_sub': "get_index_dtype not called",
#    'test_str_maxprint': "get_index_dtype not called",
#    'test_str': "get_index_dtype not called",
#    'test_star_vs_at_sign_for_sparray_and_spmatrix': "get_index_dtype not called",
#    'test_sparsity_modifying_assignment': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_slicing_3': "get_index_dtype not called",
#    'test_slicing_2': "get_index_dtype not called",
#    'test_slice_scalar_assign': "get_index_dtype not called",
#    'test_slice_assignment': "get_index_dtype not called",
#    'test_slice_assign_2': "get_index_dtype not called",
#    'test_size_zero_matrix_arithmetic': "get_index_dtype not called",
#    'test_size_zero_conversions': "get_index_dtype not called",
#    'test_setelement': "get_index_dtype not called",
#    'test_setdiag_comprehensive': "get_index_dtype not called",
#    'test_setdiag': "get_index_dtype not called",
#    'test_set_slice': "get_index_dtype not called",
#    'test_sequence_assignment': "get_index_dtype not called",
#    'test_self_self_assignment': "get_index_dtype not called",
#    'test_scalar_assign_2': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_round': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_resize': "get_index_dtype not called",
#    'test_reshape_copy': "get_index_dtype not called",
#    'test_reshape': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_real': "get_index_dtype not called",
#    'test_radd': "get_index_dtype not called",
#    'test_pow': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_nonzero': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_non_unit_stride_2d_indexing': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_neg': "get_index_dtype not called",
#    'test_ne': "get_index_dtype not called",
#    'test_nanminmax': "get_index_dtype not called",
#    'test_multiple_ellipsis_slicing': "get_index_dtype not called",
#    'test_mul_scalar': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_mu': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_minmax_axis': "get_index_dtype not called",
#    'test_minmax': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmul': "get_index_dtype not called",
#    'test_matmat_sparse': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_lt': "get_index_dtype not called",
#    'test_le': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_inplace_dense': "get_index_dtype not called",
#    'test_index_scalar_assign': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_gt': "get_index_dtype not called",
#    'test_getrow': "get_index_dtype not called",
#    'test_getnnz_axis': "get_index_dtype not called",
#    'test_getelement': "get_index_dtype not called",
#    'test_getcol': "get_index_dtype not called",
#    'test_get_vert_slice': "get_index_dtype not called",
#    'test_get_slices': "get_index_dtype not called",
#    'test_get_horiz_slice': "get_index_dtype not called",
#    'test_ge': "get_index_dtype not called",
#    'test_from_sparse': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_sparse_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_set': "get_index_dtype not called",
#    'test_fancy_indexing_seq_assign': "get_index_dtype not called",
#    'test_fancy_indexing_randomized': "get_index_dtype not called",
#    'test_fancy_indexing_ndarray': "get_index_dtype not called",
#    'test_fancy_indexing_multidim_set': "get_index_dtype not called",
#    'test_fancy_indexing_empty': "get_index_dtype not called",
#    'test_fancy_indexing_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_2d_assign': "get_index_dtype not called",
#    'test_fancy_indexing': "get_index_dtype not called",
#    'test_fancy_assignment_dtypes': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_fancy_assign_empty': "get_index_dtype not called",
#    'test_eq': "get_index_dtype not called",
#    'test_empty_arithmetic': "get_index_dtype not called",
#    'test_empty': "get_index_dtype not called",
#    'test_ellipsis_slicing': "get_index_dtype not called",
#    'test_elementwise_power': "get_index_dtype not called",
#    'test_elementwise_multiply_broadcast': "get_index_dtype not called",
#    'test_elementwise_divide': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dtype_preservation_empty_index': "get_index_dtype not called",
#    'test_dtype_preservation': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_diagonal': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_constructor3': "get_index_dtype not called",
#    'test_constructor1_base': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bool_rollover': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_assign_empty': "get_index_dtype not called",
#    'test_assign_1d_slice': "get_index_dtype not called",
#    'test_argmax': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_add0': "get_index_dtype not called",
#    'test_add': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#    'test_unary_ufunc_overrides': "get_index_dtype not called",
#    'test_todense': "get_index_dtype not called",
#    'test_toarray': "get_index_dtype not called",
#    'test_sum_invalid_params': "get_index_dtype not called",
#    'test_sum_dtype': "get_index_dtype not called",
#    'test_sum': "get_index_dtype not called",
#    'test_sub_dense': "get_index_dtype not called",
#    'test_str_maxprint': "get_index_dtype not called",
#    'test_str': "get_index_dtype not called",
#    'test_star_vs_at_sign_for_sparray_and_spmatrix': "get_index_dtype not called",
#    'test_sparsity_modifying_assignment': "get_index_dtype not called",
#    'test_small_multiplication': "get_index_dtype not called",
#    'test_slicing_3': "get_index_dtype not called",
#    'test_slicing_2': "get_index_dtype not called",
#    'test_slice_scalar_assign': "get_index_dtype not called",
#    'test_slice_assignment': "get_index_dtype not called",
#    'test_slice_assign_2': "get_index_dtype not called",
#    'test_setelement': "get_index_dtype not called",
#    'test_set_slice': "get_index_dtype not called",
#    'test_sequence_assignment': "get_index_dtype not called",
#    'test_self_self_assignment': "get_index_dtype not called",
#    'test_scalar_assign_2': "get_index_dtype not called",
#    'test_rsub': "get_index_dtype not called",
#    'test_rmul_scalar_type_error': "get_index_dtype not called",
#    'test_rmul_scalar': "get_index_dtype not called",
#    'test_rmatvec': "get_index_dtype not called",
#    'test_resize': "get_index_dtype not called",
#    'test_repr': "get_index_dtype not called",
#    'test_real': "get_index_dtype not called",
#    'test_radd': "get_index_dtype not called",
#    'test_pickle': "get_index_dtype not called",
#    'test_numpy_sum': "get_index_dtype not called",
#    'test_numpy_nonzero': "get_index_dtype not called",
#    'test_numpy_minmax': "get_index_dtype not called",
#    'test_numpy_mean': "get_index_dtype not called",
#    'test_nonzero': "get_index_dtype not called",
#    'test_negative_index_assignment': "get_index_dtype not called",
#    'test_neg': "get_index_dtype not called",
#    'test_nanminmax': "get_index_dtype not called",
#    'test_multiple_ellipsis_slicing': "get_index_dtype not called",
#    'test_mul_scalar': "get_index_dtype not called",
#    'test_mul_custom_type': "get_index_dtype not called",
#    'test_missized_masking': "get_index_dtype not called",
#    'test_minmax_invalid_params': "get_index_dtype not called",
#    'test_mean_out': "get_index_dtype not called",
#    'test_mean_invalid_params': "get_index_dtype not called",
#    'test_mean_dtype': "get_index_dtype not called",
#    'test_mean': "get_index_dtype not called",
#    'test_matvec': "get_index_dtype not called",
#    'test_matmat_dense': "get_index_dtype not called",
#    'test_invalid_shapes': "get_index_dtype not called",
#    'test_inplace_success': "get_index_dtype not called",
#    'test_inplace_dense': "get_index_dtype not called",
#    'test_index_scalar_assign': "get_index_dtype not called",
#    'test_imul_scalar': "get_index_dtype not called",
#    'test_imag': "get_index_dtype not called",
#    'test_idiv_scalar': "get_index_dtype not called",
#    'test_getrow': "get_index_dtype not called",
#    'test_getnnz_axis': "get_index_dtype not called",
#    'test_get_vert_slice': "get_index_dtype not called",
#    'test_get_slices': "get_index_dtype not called",
#    'test_get_horiz_slice': "get_index_dtype not called",
#    'test_from_matrix': "get_index_dtype not called",
#    'test_from_list': "get_index_dtype not called",
#    'test_from_array': "get_index_dtype not called",
#    'test_fancy_indexing_sparse_boolean': "get_index_dtype not called",
#    'test_fancy_indexing_set': "get_index_dtype not called",
#    'test_fancy_indexing_regression_3087': "get_index_dtype not called",
#    'test_fancy_indexing_ndarray': "get_index_dtype not called",
#    'test_fancy_indexing_multidim_set': "get_index_dtype not called",
#    'test_fancy_indexing_empty': "get_index_dtype not called",
#    'test_fancy_indexing_boolean': "get_index_dtype not called",
#    'test_fancy_indexing': "get_index_dtype not called",
#    'test_fancy_assignment_dtypes': "get_index_dtype not called",
#    'test_fancy_assign_slice': "get_index_dtype not called",
#    'test_fancy_assign_ndarray': "get_index_dtype not called",
#    'test_fancy_assign_list': "get_index_dtype not called",
#    'test_fancy_assign_empty': "get_index_dtype not called",
#    'test_empty_arithmetic': "get_index_dtype not called",
#    'test_empty': "get_index_dtype not called",
#    'test_ellipsis_slicing': "get_index_dtype not called",
#    'test_dtype_preservation_empty_slice': "get_index_dtype not called",
#    'test_dtype_preservation_empty_index': "get_index_dtype not called",
#    'test_dtype_preservation': "get_index_dtype not called",
#    'test_dot_scalar': "get_index_dtype not called",
#    'test_default_dtype': "get_index_dtype not called",
#    'test_count_nonzero': "get_index_dtype not called",
#    'test_copy': "get_index_dtype not called",
#    'test_constructor4': "get_index_dtype not called",
#    'test_constructor1_base': "get_index_dtype not called",
#    'test_comparisons_custom_type': "get_index_dtype not called",
#    'test_bsr_matvec': "get_index_dtype not called",
#    'test_bool_rollover': "get_index_dtype not called",
#    'test_bool': "get_index_dtype not called",
#    'test_binop_custom_type_with_shape': "get_index_dtype not called",
#    'test_binop_custom_type': "get_index_dtype not called",
#    'test_bad_index_assign': "get_index_dtype not called",
#    'test_bad_index': "get_index_dtype not called",
#    'test_astype_immutable': "get_index_dtype not called",
#    'test_astype': "get_index_dtype not called",
#    'test_assign_empty': "get_index_dtype not called",
#    'test_assign_1d_slice': "get_index_dtype not called",
#    'test_asfptype': "get_index_dtype not called",
#    'test_add_dense': "get_index_dtype not called",
#    'test_add0': "get_index_dtype not called",
#    'test_add': "get_index_dtype not called",
#    'test_abs': "get_index_dtype not called",
#}
