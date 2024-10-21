""" Test functions involving 64bit or 32bit indexing """
import pytest
import numpy as np

from .test_base import (
    TestBSR,
    TestCOO,
    TestCSC,
    TestCSR,
    TestDIA,
    TestDOK,
    TestLIL,
    TestBSRMatrix,
    TestCOOMatrix,
    TestCSCMatrix,
    TestCSRMatrix,
    TestDIAMatrix,
    TestDOKMatrix,
    TestLILMatrix,
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
    if sp_api == "sparray":
        TEST_CLASSES = [TestBSR, TestCOO, TestCSC, TestCSR, TestDIA]
    elif sp_api == "sparray-extra":
        # lil/dok->other conversion operations use get_index_dtype
        # so we include lil & dok test suite even though they do not
        # use get_index_dtype within the class. That means many of
        # these tests are superfluous, but it's hard to pick which
        TEST_CLASSES = [TestDOK, TestLIL]
    elif sp_api == "spmatrix":
        TEST_CLASSES = [TestBSRMatrix, TestCOOMatrix, TestCSCMatrix,
                        TestCSRMatrix, TestDIAMatrix,]
    elif sp_api == "spmatrix-extra":
        # lil/dok->other conversion operations use get_index_dtype
        TEST_CLASSES = [TestDOKMatrix, TestLILMatrix]
    else:
        raise ValueError(f"parameter {sp_api=} is not valid")

    for cls in TEST_CLASSES:
        for method_name in sorted(dir(cls)):
            method = getattr(cls, method_name)
            if (method_name.startswith('test_') and
                    not getattr(method, 'slow', False)):
                marks = []

                msg = SKIP_TESTS.get(method_name)
                if bool(msg):
                    marks += [pytest.mark.skip(reason=msg)]

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
