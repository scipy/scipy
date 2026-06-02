cdef public int foo() noexcept:
    return 1

cdef public int foo_deprecated() noexcept:
    return 1

from scipy._lib.deprecation import deprecate_cython_api
import scipy._lib._test_deprecation_def as mod
deprecate_cython_api(mod, "foo_deprecated", new_name="foo",
                     message="Deprecated in Scipy 42.0.0")
del deprecate_cython_api, mod
