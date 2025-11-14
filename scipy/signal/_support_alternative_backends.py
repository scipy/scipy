import functools
from scipy._lib._array_api import (
    is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API, xp_capabilities
)

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
    'lfilter_zi', 'sosfilt_zi', 'get_window', 'besselap', 'envelope', 'remez', 'bessel'
]

# freqz_sos is a sosfreqz rename, and cupy does not have the new name yet (in v13.x)
CUPY_RENAMES = {'freqz_sos': 'sosfreqz'}


def delegate_xp(delegator, module_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwds):
            try:
                xp = delegator(*args, **kwds)
            except TypeError:
                # object arrays
                import numpy as np
                xp = np

            # try delegating to a cupyx/jax namesake
            if is_cupy(xp) and func.__name__ not in CUPY_BLACKLIST:
                func_name = CUPY_RENAMES.get(func.__name__, func.__name__)

                # https://github.com/cupy/cupy/issues/8336
                import importlib
                cupyx_module = importlib.import_module(f"cupyx.scipy.{module_name}")
                cupyx_func = getattr(cupyx_module, func_name)
                kwds.pop('xp', None)
                return cupyx_func(*args, **kwds)
            elif is_jax(xp) and func.__name__ in JAX_SIGNAL_FUNCS:
                spx = scipy_namespace_for(xp)
                jax_module = getattr(spx, module_name)
                jax_func = getattr(jax_module, func.__name__)
                kwds.pop('xp', None)
                return jax_func(*args, **kwds)
            else:
                # the original function
                return func(*args, **kwds)
        return wrapper
    return inner


# Although most of these functions currently exist in CuPy and some in JAX,
# there are no alternative backend tests for any of them in the current
# test suite. Each will be documented as np_only until tests are added.
untested = {
    "argrelextrema",
    "argrelmax",
    "argrelmin",
    "band_stop_obj",
    "check_COLA",
    "check_NOLA",
    "chirp",
    "coherence",
    "csd",
    "czt_points",
    "dbode",
    "dfreqresp",
    "dlsim",
    "dstep",
    "find_peaks",
    "find_peaks_cwt",
    "findfreqs",
    "freqresp",
    "gausspulse",
    "istft",
    "lombscargle",
    "lsim",
    "max_len_seq",
    "peak_prominences",
    "peak_widths",
    "periodogram",
    "place_pols",
    "sawtooth",
    "sepfir2d",
    "spectrogram",
    "square",
    "ss2tf",
    "ss2zpk",
    "step",
    "stft",
    "sweep_poly",
    "symiirorder1",
    "symiirorder2",
    "tf2ss",
    "unit_impulse",
    "welch",
    "zoom_fft",
    "zpk2ss",
}


def get_default_capabilities(func_name, delegator):
    if delegator is None or func_name in untested:
        return xp_capabilities(np_only=True)
    return xp_capabilities()


capabilities_overrides = {
    "bessel": xp_capabilities(cpu_only=True, jax_jit=False, allow_dask_compute=True),
    "besselap": xp_capabilities(np_only=True),
    "bilinear": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                jax_jit=False, allow_dask_compute=True),
    "bilinear_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                    jax_jit=False, allow_dask_compute=True),
    "butter": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False, 
                              allow_dask_compute=True),
    "buttord": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               jax_jit=False, allow_dask_compute=True),
    "cheby1": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "cheb1ord": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                jax_jit=False, allow_dask_compute=True),
    "cheby2": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "cheb2ord": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                jax_jit=False, allow_dask_compute=True),
    "cont2discrete": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "convolve": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                 allow_dask_compute=True),
    "convolve2d": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                 allow_dask_compute=True),
    "correlate": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                 allow_dask_compute=True),
    "correlate2d": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                   allow_dask_compute=True),
    "correlation_lags": xp_capabilities(out_of_scope=True),
    "cspline1d": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                 jax_jit=False, allow_dask_compute=True),
    "cspline1d_eval": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                      jax_jit=False, allow_dask_compute=True),
    "czt": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "deconvolve": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                  allow_dask_compute=True,
                                  skip_backends=[("jax.numpy", "item assignment")]),
    "decimate": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "detrend": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                               allow_dask_compute=True),
    "dlti": xp_capabilities(np_only=True,
                            reason="works in CuPy but delegation isn't set up yet"),
    "ellip": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False, 
                             allow_dask_compute=True),
    "ellipord": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                jax_jit=False, allow_dask_compute=True),
    "firls": xp_capabilities(cpu_only=True, allow_dask_compute=True, jax_jit=False,
                             reason="lstsq"),
    "firwin": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                              jax_jit=False, allow_dask_compute=True),
    "firwin2": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               jax_jit=False, allow_dask_compute=True),
    "fftconvolve": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"]),
    "freqs": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             jax_jit=False, allow_dask_compute=True),
    "freqs_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "freqz": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             jax_jit=False, allow_dask_compute=True),
    "freqz_sos": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "group_delay": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                   jax_jit=False, allow_dask_compute=True),
    "hilbert": xp_capabilities(
        cpu_only=True, exceptions=["cupy", "torch"],
        skip_backends=[("jax.numpy", "item assignment")],
    ),
    "hilbert2": xp_capabilities(
        cpu_only=True, exceptions=["cupy", "torch"],
        skip_backends=[("jax.numpy", "item assignment")],
    ),
    "invres": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "invresz": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "iircomb": xp_capabilities(xfail_backends=[("jax.numpy", "inaccurate")]),
    "iirfilter": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "lfiltic": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True),
    "lfilter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True, jax_jit=False),
    "lfilter_zi": xp_capabilities(cpu_only=True, allow_dask_compute=True,
                                  jax_jit=False),
    "lp2bp": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True,
                             skip_backends=[("jax.numpy", "in-place item assignment")]),
    "lp2bp_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lp2bs": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True,
                             skip_backends=[("jax.numpy", "in-place item assignment")]),
    "lp2bs_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lp2lp": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True, jax_jit=False),
    "lp2lp_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "lp2hp": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                             allow_dask_compute=True,
                             skip_backends=[("jax.numpy", "in-place item assignment")]),
    "lp2hp_zpk": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 allow_dask_compute=True, jax_jit=False),
    "minimum_phase": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                     allow_dask_compute=True, jax_jit=False),
    "medfilt": xp_capabilities(allow_dask_compute=True, jax_jit=False),
    "medfilt2d": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                 allow_dask_compute=True, jax_jit=False),
    "normalize": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "oaconvolve": xp_capabilities(
        cpu_only=True, exceptions=["cupy", "torch"],
        skip_backends=[("jax.numpy", "fails all around")],
        xfail_backends=[("dask.array", "wrong answer")],
    ),
    "order_filter": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                                    allow_dask_compute=True, jax_jit=False),
    "qspline1d": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                 jax_jit=False, allow_dask_compute=True),
    "qspline1d_eval": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                      jax_jit=False, allow_dask_compute=True),
    "remez": xp_capabilities(cpu_only=True, allow_dask_compute=True, jax_jit=False),
    "resample": xp_capabilities(
        cpu_only=True, exceptions=["cupy"],
        jax_jit=False, skip_backends=[("dask.array", "XXX something in dask")]
    ),
    "resample_poly": xp_capabilities(
        cpu_only=True, exceptions=["cupy"],
        jax_jit=False, skip_backends=[("dask.array", "XXX something in dask")]
    ),
    "residue": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "residuez": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "savgol_filter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                     reason="convolve1d is cpu-only"),
    "sos2zpk": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                               allow_dask_compute=True),
    "sos2tf": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "sosfilt": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                               allow_dask_compute=True),
    "sosfiltfilt": xp_capabilities(
        cpu_only=True, exceptions=["cupy"],
        skip_backends=[
            (
                "dask.array",
                "sosfiltfilt directly sets shape attributes on arrays"
                " which dask doesn't like"
            ),
            ("torch", "negative strides"),
            ("jax.numpy", "sosfilt works in-place"),
        ],
    ),
    "spline_filter": xp_capabilities(cpu_only=True, exceptions=["cupy"],
                                     jax_jit=False, allow_dask_compute=True),
    "sosfreqz": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                 jax_jit=False, allow_dask_compute=True),
    "tf2sos": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "tf2zpk": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "unique_roots": xp_capabilities(np_only=True, exceptions=["cupy"]),
    "upfirdn": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                               allow_dask_compute=True,
                               reason="Cython implementation"),
    "vectorstrength": xp_capabilities(cpu_only=True, exceptions=["cupy", "torch"],
                                      allow_dask_compute=True, jax_jit=False),
    "wiener": xp_capabilities(cpu_only=True, exceptions=["cupy", "jax.numpy"],
                              allow_dask_compute=True, jax_jit=False),
    "zpk2sos": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
    "zpk2tf": xp_capabilities(cpu_only=True, exceptions=["cupy"], jax_jit=False,
                              allow_dask_compute=True),
}
        
        
# ### decorate ###
for obj_name in _signal_api.__all__:
    bare_obj = getattr(_signal_api, obj_name)
    delegator = getattr(_delegators, obj_name + "_signature", None)

    if SCIPY_ARRAY_API and delegator is not None:
        f = delegate_xp(delegator, MODULE_NAME)(bare_obj)
    else:
        f = bare_obj

    if callable(f) and hasattr(f, "__name__"):
        capabilities = capabilities_overrides.get(
            obj_name, get_default_capabilities(obj_name, delegator)
        )
        f = capabilities(f)

    # add the decorated function to the namespace, to be imported in __init__.py
    vars()[obj_name] = f
