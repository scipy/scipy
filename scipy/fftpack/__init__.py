#
# fftpack - Discrete Fourier Transform algorithms.
#
# Created: Pearu Peterson, August,September 2002

from info import __all__,__doc__

from fftpack_version import fftpack_version as __version__

from basic import *
from pseudo_diffs import *
from helper import *

from numpy.dual import register_func
for k in ['fft', 'ifft', 'fftn', 'ifftn', 'fft2', 'ifft2']:
    register_func(k, eval(k))
del k, register_func

from realtransforms import *
__all__.extend(['dct', 'idct'])

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
