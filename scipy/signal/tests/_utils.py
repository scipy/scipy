import pytest
import numpy
from scipy._lib._array_api_util import _NUMPY_INCLUDES_ARRAY_API


if _NUMPY_INCLUDES_ARRAY_API:
    import numpy.array_api

pytest_enable_array_api = pytest.mark.parametrize(
    "xp", [numpy, *((numpy.array_api,) if _NUMPY_INCLUDES_ARRAY_API else ())]
)
