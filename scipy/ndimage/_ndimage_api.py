"""This is the 'bare' ndimage API.

This --- private! --- module only collects implementations of public ndimage API
for _support_alternative_backends.
The latter --- also private! --- module adds delegation to CuPy etc and
re-exports decorated names to __init__.py
"""

from . import _filters
from ._filters import *    # noqa: F403
from . import _fourier
from ._fourier import *   # noqa: F403
from . import _interpolation
from ._interpolation import *   # noqa: F403
from . import _measurements
from ._measurements import *   # noqa: F403
from . import _morphology
from ._morphology import *   # noqa: F403

__all__ = []
__all__ += _filters.__all__
__all__ += _fourier.__all__
__all__ += _interpolation.__all__
__all__ += _measurements.__all__
__all__ += _morphology.__all__
