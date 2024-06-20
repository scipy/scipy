import pytest

from scipy.signal._support_alternative_backends import array_special_func_map
from scipy.conftest import array_api_compatible
from scipy import signal
from scipy._lib._array_api import xp_assert_close, is_jax, get_array_subpackage_func
from scipy._lib.array_api_compat import numpy as np

try:
    import array_api_strict
    HAVE_ARRAY_API_STRICT = True
except ImportError:
    HAVE_ARRAY_API_STRICT = False


@pytest.mark.fail_slow(5)
@array_api_compatible
# @pytest.mark.skip_xp_backends('numpy', reasons=['skip while debugging'])
# @pytest.mark.usefixtures("skip_xp_backends")
# `reversed` is for developer convenience: test new function first = less waiting
@pytest.mark.parametrize('f_name_n_args', reversed(array_special_func_map.items()))
@pytest.mark.parametrize('dtype', ['float32', 'float64'])
@pytest.mark.parametrize('shapes', [[100, 150]])
def test_support_alternative_backends(xp, f_name_n_args, dtype, shapes):
    f_name, n_args = f_name_n_args
    shapes = shapes[:n_args]
    f = getattr(signal, f_name)

    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)

    rng = np.random.default_rng(984254252920492019)
    args_np = [rng.standard_normal(size=shape, dtype=dtype_np) for shape in shapes]
    args_xp = [xp.asarray(arg[()], dtype=dtype_xp) for arg in args_np]

    res = f(*args_xp)
    ref = xp.asarray(f(*args_np), dtype=dtype_xp)

    eps = np.finfo(dtype_np).eps
    xp_assert_close(res, ref, atol=10*eps)
