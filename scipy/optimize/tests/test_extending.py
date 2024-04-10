import os
import platform

import pytest

from scipy._lib._testutils import _test_cython_extension, cython


@pytest.mark.skipif(platform.machine() in ["wasm32", "wasm64"],
                    reason="Can't start subprocess")
@pytest.mark.skipif(cython is None, reason="requires cython")
def test_cython(tmp_path):
    srcdir = os.path.dirname(os.path.dirname(__file__))
    extensions, extensions_cpp = _test_cython_extension(tmp_path, srcdir)
    # actually test the cython c-extensions
    # From docstring for scipy.optimize.cython_optimize module
    x = extensions.brentq_example()
    assert x == 0.6999942848231314
    x = extensions_cpp.brentq_example()
    assert x == 0.6999942848231314
