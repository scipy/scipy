# Last Change: Tue Dec 12 05:00 PM 2006 J
from info import __doc__

from lpc import lpc2 as lpc
from autocorr import autocorr_oneside_nofft, autocorr_fft
from segmentaxis import segment_axis

from numpy.testing import Tester
test = Tester().test
