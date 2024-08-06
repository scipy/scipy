import os
from functools import wraps

from hypothesis import strategies as st
from hypothesis.extra import array_api

from . import _version

__all__ = ["xp", "api_version", "xps"]


from scipy import sparse as sp
xp = sp.csr_array(sp.eye(2)).__array_namespace__()

# We monkey patch floats() to always disable subnormals as they are out-of-scope

_floats = st.floats


@wraps(_floats)
def floats(*a, **kw):
    kw["allow_subnormal"] = False
    return _floats(*a, **kw)


st.floats = floats


# We do the same with xps.from_dtype() - this is not strictly necessary, as
# the underlying floats() will never generate subnormals. We only do this
# because internal logic in xps.from_dtype() assumes xp.finfo() has its
# attributes as scalar floats, which is expected behaviour but disrupts many
# unrelated tests.
try:
    __from_dtype = array_api._from_dtype

    @wraps(__from_dtype)
    def _from_dtype(*a, **kw):
        kw["allow_subnormal"] = False
        return __from_dtype(*a, **kw)

    array_api._from_dtype = _from_dtype
except AttributeError:
    # Ignore monkey patching if Hypothesis changes the private API
    pass


api_version = os.getenv(
    "ARRAY_API_TESTS_VERSION", getattr(xp, "__array_api_version__", "2021.12")
)
xps = array_api.make_strategies_namespace(xp, api_version=api_version)

__version__ = _version.get_versions()["version"]
