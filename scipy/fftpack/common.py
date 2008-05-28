"""This module takes care of exposing the fft implementation from available
backends to scipy.fftpack.

For each backend, the exposed implementation consists solely of the items in
__all__ .

The functions are set up in load_backend function, which initialize the
dictionary _FUNCS with working functions. If a backend does not implement one
of the function, a fallback from _DEF_BACKEND is automatically used."""

# Default backend: fftpack.
import fftpack as _DEF_BACKEND
_DEF_BACKEND_NAME = "fftpack"

__all__ = ["zfft", "drfft", "zfftnd", "zrfft", "init_convolution_kernel",
           "convolve", "convolve_z", "destroy_convolve_cache", "backend",
           "detailed_backend"]

_FFT_FUNCNAME = ["zfft", "drfft", "zfftnd", "zrfft"]

# Convolve needs to be treated as a whole, not per function, we keep the
# convolve functions separately.
_CONVOLVE_FUNCNAME = ["init_convolution_kernel",
                      "convolve", "convolve_z", "destroy_convolve_cache"]

# Dictionary of (name, callable)
_FUNCS = dict([(name, None) for name in _FFT_FUNCNAME + _CONVOLVE_FUNCNAME])

_FALLBACK = dict([(name, _DEF_BACKEND.__dict__[f]) for f in _FUNCS.keys()])

# Dictionary to keep func -> used backend association
_BACK_PER_FUNC = {}

def myimport(name):
    """Load a fft backend from its name.

    Name should be fftw3, etc..."""
    mod = __import__("scipy.fftpack.backends", fromlist = [name])

    # Because of nose trying to import backends, we do not generate ImporError
    # in backends, but we only set IS_INIT to true if the backend is available
    try:
        ret = mod.__dict__[name]
        if ret.IS_INIT:
            return ret
        else:
            raise ImportError("Could not import %s (was not initialized)"\
                              % name)
    except KeyError, e:
        raise ImportError(e)

def find_backend():
    """Try to import one of the backend from the list of possible backends, and
    return the first found.

    Returns the backend (module class), raise ImportError if nothing is
    found."""
    for backend in ["mkl", "fftw3", "fftw"]:
        try:
            mod = myimport(backend)
            return mod, backend
        except ImportError, e:
            pass
    raise ImportError("No backend found")

def load_backend(name = None):
    """Load backend and set functions accordingly.

    If name is None, all backend are tried one after the other, as defined in
    function find_backend.

    If the backend does not implement one function, the implementation in
    _DEF_BACKEND is used instead.

    If no optimized backend is found,  fallback is used for all the
    functions."""
    try:
        if name:
            mod = myimport(name)
        else:
            mod, name = find_backend()

        # Loading fft functions: each of them can be loaded independently
        for f in _FFT_FUNCNAME:
            try:
                _FUNCS[f] = mod.__dict__[f]
                _BACK_PER_FUNC[f] = name
            except KeyError:
                _FUNCS[f] = _DEF_BACKEND.__dict__[f]
                _BACK_PER_FUNC[f] = _DEF_BACKEND_NAME

        # Loading convolve: we try to load all of them: if any failure, we use
        # fallback for all of them
        try:
            for f in _CONVOLVE_FUNCNAME:
                _FUNCS[f] = mod.__dict__[f]
                _BACK_PER_FUNC[f] = name
        except KeyError:
            for f in _CONVOLVE_FUNCNAME:
                _FUNCS[f] = _DEF_BACKEND.__dict__[f]
                _BACK_PER_FUNC[f] = _DEF_BACKEND_NAME
    except ImportError, e:
        name = _DEF_BACKEND_NAME
        # If cannot load backend, just use default backend (fftpack)
        for f in _FUNCS.keys():
            _FUNCS[f] = _DEF_BACKEND.__dict__[f]
            _BACK_PER_FUNC[f] = _DEF_BACKEND_NAME

    return name

def backend():
    """Return backend name."""
    return _BACKEND_NAME

def detailed_backend():
    """Return backend name for each function.
    
    The return value is a dictionary of (func name, backend)"""
    return _BACK_PER_FUNC

_BACKEND_NAME = load_backend()

zfft = _FUNCS["zfft"]
drfft = _FUNCS["drfft"]
zfftnd = _FUNCS["zfftnd"]
zrfft = _FUNCS["zrfft"]
init_convolution_kernel = _FUNCS["init_convolution_kernel"]
convolve = _FUNCS["convolve"]
convolve_z = _FUNCS["convolve_z"]
destroy_convolve_cache = _FUNCS["destroy_convolve_cache"]
