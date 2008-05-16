"""This module should never be imported directly."""

# We wrap around ImportError because nosetest tries to import this
# module and generates an error when the backend is not built
try:
    from _mkl import zfft_mkl as zfft, \
                     zfftnd_mkl as zfftnd
    IS_INIT = True
except ImportError:
    IS_INIT = False
