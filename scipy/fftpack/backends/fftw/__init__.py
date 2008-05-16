"""This module should never be imported directly."""

# We wrap around ImportError because nosetest tries to import this
# module and generates an error when the backend is not built
try:
    from _fftw import zfft_fftw as zfft, \
                      drfft_fftw as drfft, \
                      zfftnd_fftw as zfftnd
except ImportError:
    pass
