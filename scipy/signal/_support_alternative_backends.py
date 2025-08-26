from scipy._lib._array_api import SCIPY_ARRAY_API

from ._signal_api import *   # noqa: F403
from . import _signal_api
from . import _delegators
__all__ = _signal_api.__all__


MODULE_NAME = 'signal'

# jax.scipy.signal has only partial coverage of scipy.signal, so we keep the list
# of functions we can delegate to JAX
# https://jax.readthedocs.io/en/latest/jax.scipy.html
JAX_SIGNAL_FUNCS = [
    'fftconvolve', 'convolve', 'convolve2d', 'correlate', 'correlate2d',
    'csd', 'detrend', 'istft', 'welch'
]

# some cupyx.scipy.signal functions are incompatible with their scipy counterparts
CUPY_BLACKLIST = [
    'lfilter_zi', 'sosfilt_zi', 'get_window', 'besselap', 'envelope', 'remez'
]

# freqz_sos is a sosfreqz rename, and cupy does not have the new name yet (in v13.x)
CUPY_RENAMES = {'freqz_sos': 'sosfreqz'}


def drop_xp_wrapper(func):
    # We could do this only if there actually is an `xp=` argument.
    def new_func(*args, **kwargs):
        kwargs.pop("xp", None)
        return func(*args, **kwargs)
    return new_func


class _JaxFunctions:
    # Small hack to wrap functions to drop `xp=`, we need it to look like
    # a module...
    @classmethod
    def __getattr__(cls, func_name):
        import jax.scipy.signal
        return drop_xp_wrapper(getattr(jax.scipy.signal, func_name))


class JaxBackend:
    name = "jax"
    # A class, just for convenience (may need to change)
    primary_types = ["~jax:Array"]  # allow subclasses otherwise need _jax.ArrayImpl?
    secondary_types = []
    requires_opt_in = False

    functions = {
        f"scipy.signal:{func}": {
            "function":
                f"scipy.signal._support_alternative_backends:"
                f"JaxBackend.JaxFunctions.{func}"}
        for func in JAX_SIGNAL_FUNCS
    }
    JaxFunctions = _JaxFunctions()


class _CupyFunctions:
    # Small hack to wrap functions to drop `xp=`, we need it to look like
    # a module...
    @classmethod
    def __getattr__(cls, func_name):
        import cupyx.signal
        return drop_xp_wrapper(getattr(cupyx.signal, func_name))


class CupyBackend:
    name = "cupy"
    # A class, just for convenience (may need to change)
    primary_types = ["cupy:ndarray"]
    secondary_types = []
    requires_opt_in = False
    # See JaxBackend about the `__getattr__` dance here.
    functions = {
        f"scipy.signal:{func}": {
            "function":
                f"scipy.signal._support_alternative_backends:"
                f"CupyBackend.CupyFunctions.{CUPY_RENAMES.get(func, func)}"}
        for func in _signal_api.__all__ if func not in CUPY_BLACKLIST
    }
    CupyFunctions = _CupyFunctions()


if SCIPY_ARRAY_API:
    from spatch.backend_system import BackendSystem

    _bs = BackendSystem(
        "scipy_backends",
        "_SCIPY_INTERNAL_BACKENDS",  # spatch env-var prefix
        default_primary_types=["numpy:ndarray"],
        backends=[JaxBackend],
    )

    # expose backend opts to scipy.signal to allow playing with it
    vars()["backend_opts"] =  _bs.backend_opts
    __all__ = __all__ + ["backend_opts"]

# ### decorate ###
for obj_name in _signal_api.__all__:
    bare_obj = getattr(_signal_api, obj_name)
    delegator = getattr(_delegators, obj_name + "_signature", None)

    if SCIPY_ARRAY_API and delegator is not None:
        f = _bs.dispatchable(delegator, module="scipy.signal")(bare_obj)
    else:
        f = bare_obj

    # add the decorated function to the namespace, to be imported in __init__.py
    vars()[obj_name] = f
