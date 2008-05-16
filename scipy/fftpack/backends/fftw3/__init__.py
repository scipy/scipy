"""This module should never be imported directly."""

# We wrap around ImportError because nosetest tries to import this
# module and generates an error when the backend is not built

try:
    from _fftw3 import zfft_fftw3 as zfft, \
                       drfft_fftw3 as drfft, \
                       zfftnd_fftw3 as zfftnd
    IS_INIT = True
except ImportError:
    IS_INIT = False
