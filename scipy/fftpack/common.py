import fftpack as _DEF_BACKEND

__all__ = ["zfft", "drfft", "zfftnd", "zrfft", "init_convolution_kernel",
           "convolve", "convolve_z", "destroy_convolve_cache"]

_FUNCS = dict([(name, None) for name in __all__])
_FALLBACK = dict([(name, _DEF_BACKEND.__dict__[f]) for f in _FUNCS.keys()])

def myimport(name):
    mod = __import__(name)
    comps = name.split('.')[1:]
    for c in comps:
        mod = getattr(mod, c)
    return mod

def load_backend(name):
    try:
        mod = myimport(name)
        for f in _FUNCS.keys():
            try:
                _FUNCS[f] = mod.__dict__[f]
                print "loading %s from %s" % (f, name)
            except KeyError:
                _FUNCS[f] = _DEF_BACKEND.__dict__[f]
                print "loading %s from %s" % (f, "def backend")
    except ImportError, e:
        print "%s: failed loading backend %s" % (e, name)
        for f in _FUNCS.keys():
            _FUNCS[f] = _DEF_BACKEND.__dict__[f]

load_backend("fftpack.backends.fftw3")

zfft = _FUNCS["zfft"]
drfft = _FUNCS["drfft"]
zfftnd = _FUNCS["zfftnd"]
zrfft = _FUNCS["zrfft"]
init_convolution_kernel = _FUNCS["init_convolution_kernel"]
convolve = _FUNCS["convolve"]
convolve_z = _FUNCS["convolve_z"]
destroy_convolve_cache = _FUNCS["destroy_convolve_cache"]
