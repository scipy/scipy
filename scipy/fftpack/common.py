import fftpack as _DEF_BACKEND

_FUNCS_NAMES = ["zfft", "drfft", "zfftnd", "zrfft", "init_convolution_kernel",
                "convolve", "convolve_z", "destroy_convolve_cache"]
_FUNCS = dict([(name, None) for name in _FUNCS_NAMES])
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

