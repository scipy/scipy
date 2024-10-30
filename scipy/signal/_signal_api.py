"""This is the 'bare' scipy.signal API.

This --- private! --- module only collects implementations of public  API
for _support_alternative_backends.
The latter --- also private! --- module adds delegation to CuPy etc and
re-exports decorated names to __init__.py
"""

from . import _sigtools, windows
from ._waveforms import *
from ._max_len_seq import max_len_seq
from ._upfirdn import upfirdn

from ._spline import sepfir2d

from ._spline_filters import *
from ._filter_design import *
from ._fir_filter_design import *
from ._ltisys import *
from ._lti_conversion import *
from ._signaltools import *
from ._savitzky_golay import savgol_coeffs, savgol_filter
from ._spectral_py import *
from ._short_time_fft import *
from ._peak_finding import *
from ._czt import *
from .windows import get_window  # keep this one in signal namespace


__all__ = [s for s in dir() if not s.startswith('_')]
