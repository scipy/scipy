"""This is the 'bare' scipy.signal API.

This --- private! --- module only collects implementations of public  API
for _support_alternative_backends.
The latter --- also private! --- module adds delegation to CuPy etc and
re-exports decorated names to __init__.py
"""

from . import _sigtools, windows         # noqa: F401
from . import _waveforms
from ._waveforms import *        # noqa: F403
from . import _max_len_seq
from ._max_len_seq import max_len_seq       # noqa: F401
from ._upfirdn import upfirdn         # noqa: F401

from ._spline import sepfir2d          # noqa: F401

from . import _spline_filters
from ._spline_filters import *         # noqa: F403
from . import _filter_design
from ._filter_design import *         # noqa: F403
from . import _fir_filter_design
from ._fir_filter_design import *         # noqa: F403
from . import _ltisys
from ._ltisys import *         # noqa: F403
from . import _lti_conversion
from ._lti_conversion import *         # noqa: F403
from . import _signaltools
from ._signaltools import *         # noqa: F403
from ._savitzky_golay import savgol_coeffs, savgol_filter  # noqa: F401
from . import _spectral_py
from ._spectral_py import *         # noqa: F403
from . import _short_time_fft
from ._short_time_fft import *         # noqa: F403
from . import _peak_finding
from ._peak_finding import *         # noqa: F403
from . import _czt
from ._czt import *         # noqa: F403
from ._whittaker import whittaker_henderson  # noqa: F401
from .windows import get_window  # keep this one in signal namespace  # noqa: F401


__all__ = []
__all__ += _waveforms.__all__
__all__ += _max_len_seq.__all__
__all__ += ['upfirdn']  # not _upfirdn.__all__ (includes private `_output_len`)
__all__ += ['sepfir2d']
__all__ += _spline_filters.__all__
__all__ += _filter_design.__all__
__all__ += _fir_filter_design.__all__
__all__ += _ltisys.__all__
__all__ += _lti_conversion.__all__
__all__ += _signaltools.__all__
__all__ += ['savgol_coeffs', 'savgol_filter']
__all__ += _spectral_py.__all__
__all__ += _short_time_fft.__all__
__all__ += _peak_finding.__all__
__all__ += _czt.__all__
__all__ += ['whittaker_henderson', 'get_window', 'windows']
