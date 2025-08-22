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


class JaxBackend:
    name = "jax"
    # A class, just for convenience (may need to change)
    primary_types = ["jax:Array"]  # presumably needs more!
    secondary_types = []
    requires_opt_in = False
    # TODO: We are not dropping `xp=` argument here.  Would need a helper (or teach spatch)
    functions = {
        f"scipy.signal:{func}": "jax.scipy.signal:{func}" for func in JAX_SIGNAL_FUNCS
    }


class CupyBackend:
    name = "cupy"
    # A class, just for convenience (may need to change)
    primary_types = ["cupy:ndarray"]
    secondary_types = []
    requires_opt_in = False
    # TODO: We are not dropping `xp=` argument here.  Would need a helper (or teach spatch)
    functions = {
        f"scipy.signal:{func}": "cupyx.signal:{CUPY_RENAMES.get(func.__name__, func.__name__)}"
        for func in _signal_api.__all__ if func not in CUPY_BLACKLIST
    }


if SCIPY_ARRAY_API:
    from spatch.backend_system import BackendSystem

    _bs = BackendSystem(
        None,  # don't load entry-points for now (no 3rd party backends)
        "_SCIPY_INTERNAL_BACKENDS",  # spatch env-var prefix
        default_primary_types=["numpy:ndarray"],
        backends=[JaxBackend],
    )


# ### decorate ###
for obj_name in _signal_api.__all__:
    bare_obj = getattr(_signal_api, obj_name)
    delegator = getattr(_delegators, obj_name + "_signature", None)

    if SCIPY_ARRAY_API and delegator is not None:
        f = _bs.dispatchable(delegator)(bare_obj)
    else:
        f = bare_obj

    # add the decorated function to the namespace, to be imported in __init__.py
    vars()[obj_name] = f
