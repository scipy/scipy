"""
Module containing private utility functions
===========================================

The ``scipy._lib`` namespace is empty (for now). Tests for all
utilities in submodules of ``_lib`` can be run with::

    from scipy import _lib
    _lib.test()

"""
from scipy._lib._testutils import PytestTester
# Windows AVX512 fix - minimal implementation
def is_windows_non_avx512():
    """Check if we're on Windows with a non-AVX512 CPU."""
    import sys
    return sys.platform == 'win32'

def skip_on_windows_non_avx512(reason="Test hangs on Windows with non-AVX512 CPUs"):
    """Decorator to skip tests on Windows systems with non-AVX512 CPUs."""
    def decorator(func):
        import pytest
        return pytest.mark.skipif(
            is_windows_non_avx512(),
            reason=reason
        )(func)
    return decorator

test = PytestTester(__name__)
del PytestTester
