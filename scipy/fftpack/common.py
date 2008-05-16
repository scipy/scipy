"""This module takes care of exposing the fft implementation from available
backends.

The exposed implementation consists solely in the __all__ items.

The functions are set up in load_backend function, which initialize the
dictionary _FUNCS with working functions. If a backend does not implement one
of the function, a fallback is automatically used."""
# Default backend: fftpack.
import fftpack as _DEF_BACKEND

__all__ = ["zfft", "drfft", "zfftnd", "zrfft", "init_convolution_kernel",
           "convolve", "convolve_z", "destroy_convolve_cache"]

_FFT_FUNCNAME = ["zfft", "drfft", "zfftnd", "zrfft"]
_CONVOLVE_FUNCNAME = ["init_convolution_kernel", 
                      "convolve", "convolve_z", "destroy_convolve_cache"]

# Dictionary of (name, callable)
_FUNCS = dict([(name, None) for name in __all__])

_FALLBACK = dict([(name, _DEF_BACKEND.__dict__[f]) for f in _FUNCS.keys()])

def myimport(name):
    """Load a fft backend from its name.
    
    Name should be fftw3, etc..."""
    mod = __import__("scipy.fftpack.backends", fromlist = [name])
    return mod.__dict__[name]

def load_backend(name):
    try:
        mod = myimport(name)
        # Loading fft functions: each of them can be loaded independently
        for f in _FFT_FUNCNAME:
            try:
                _FUNCS[f] = mod.__dict__[f]
                print "loading %s from %s" % (f, name)
            except KeyError:
                _FUNCS[f] = _DEF_BACKEND.__dict__[f]
                print "loading %s from %s" % (f, "def backend")

        # Loading convolve: we try to load all of them: if any failure, we use
        # fallback for all of them
        try:
            for f in _CONVOLVE_FUNCNAME:
                _FUNCS[f] = mod.__dict__[f]
                print "loading %s from %s" % (f, name)
        except KeyError:
            for f in _CONVOLVE_FUNCNAME:
                _FUNCS[f] = _DEF_BACKEND.__dict__[f]
                print "loading %s from %s" % (f, "def backend")
    except ImportError, e:
        # If cannot load backend, just use default backend (fftpack)
        print "%s: failed loading backend %s" % (e, name)
        for f in _FUNCS.keys():
            _FUNCS[f] = _DEF_BACKEND.__dict__[f]

load_backend("fftw3")

zfft = _FUNCS["zfft"]
drfft = _FUNCS["drfft"]
zfftnd = _FUNCS["zfftnd"]
zrfft = _FUNCS["zrfft"]
init_convolution_kernel = _FUNCS["init_convolution_kernel"]
convolve = _FUNCS["convolve"]
convolve_z = _FUNCS["convolve_z"]
destroy_convolve_cache = _FUNCS["destroy_convolve_cache"]
