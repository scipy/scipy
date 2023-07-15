import functools
import warnings
from importlib import import_module


__all__ = ["_deprecated"]


# Object to use as default value for arguments to be deprecated. This should
# be used over 'None' as the user could parse 'None' as a positional argument
_NoValue = object()

def _sub_module_deprecation(*, sub_package, module, private_module, all,
                            attribute):
    """Helper function for deprecating modules that are public but were
    intended to be private.

    Parameters
    ----------
    sub_package : str
        Subpackage the module belongs to eg. stats
    module : str
        Public but intended private module to deprecate
    private_module : str
        Private replacement for `module`
    all : list
        ``__all__`` belonging to `module`
    attribute : str
        The attribute in `module` being accessed
    """
    if attribute not in all:
        raise AttributeError(
            f"`scipy.{sub_package}.{module}` has no attribute `{attribute}`; furthermore, "
            f"`scipy.{sub_package}.{module}` is deprecated and will be removed in "
            "SciPy 2.0.0.")

    attr = getattr(import_module(f"scipy.{sub_package}"), attribute, None)

    if attr is not None:
        message = (f"Please import `{attribute}` from the `scipy.{sub_package}` namespace; "
                   f"the `scipy.{sub_package}.{module}` namespace is deprecated and "
                   "will be removed in SciPy 2.0.0.")
    else:
        message = (f"`scipy.{sub_package}.{module}.{attribute}` is deprecated along with "
                   f"the `scipy.{sub_package}.{module}` namespace. "
                   f"`scipy.{sub_package}.{module}.{attribute}` will be removed in SciPy 1.13.0, and "
                   f"the `scipy.{sub_package}.{module}` namespace will be removed in SciPy 2.0.0.")

    warnings.warn(message, category=DeprecationWarning, stacklevel=3)

    return getattr(import_module(f"scipy.{sub_package}.{private_module}"), attribute)
    

def _deprecated(msg, stacklevel=2):
    """Deprecate a function by emitting a warning on use."""
    def wrap(fun):
        if isinstance(fun, type):
            warnings.warn(
                f"Trying to deprecate class {fun!r}",
                category=RuntimeWarning, stacklevel=2)
            return fun

        @functools.wraps(fun)
        def call(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning,
                          stacklevel=stacklevel)
            return fun(*args, **kwargs)
        call.__doc__ = fun.__doc__
        return call

    return wrap


class _DeprecationHelperStr:
    """
    Helper class used by deprecate_cython_api
    """
    def __init__(self, content, message):
        self._content = content
        self._message = message

    def __hash__(self):
        return hash(self._content)

    def __eq__(self, other):
        res = (self._content == other)
        if res:
            warnings.warn(self._message, category=DeprecationWarning,
                          stacklevel=2)
        return res


def deprecate_cython_api(module, routine_name, new_name=None, message=None):
    """
    Deprecate an exported cdef function in a public Cython API module.

    Only functions can be deprecated; typedefs etc. cannot.

    Parameters
    ----------
    module : module
        Public Cython API module (e.g. scipy.linalg.cython_blas).
    routine_name : str
        Name of the routine to deprecate. May also be a fused-type
        routine (in which case its all specializations are deprecated).
    new_name : str
        New name to include in the deprecation warning message
    message : str
        Additional text in the deprecation warning message

    Examples
    --------
    Usually, this function would be used in the top-level of the
    module ``.pyx`` file:

    >>> from scipy._lib.deprecation import deprecate_cython_api
    >>> import scipy.linalg.cython_blas as mod
    >>> deprecate_cython_api(mod, "dgemm", "dgemm_new",
    ...                      message="Deprecated in Scipy 1.5.0")
    >>> del deprecate_cython_api, mod

    After this, Cython modules that use the deprecated function emit a
    deprecation warning when they are imported.

    """
    old_name = f"{module.__name__}.{routine_name}"

    if new_name is None:
        depdoc = "`%s` is deprecated!" % old_name
    else:
        depdoc = "`%s` is deprecated, use `%s` instead!" % \
                 (old_name, new_name)

    if message is not None:
        depdoc += "\n" + message

    d = module.__pyx_capi__

    # Check if the function is a fused-type function with a mangled name
    j = 0
    has_fused = False
    while True:
        fused_name = f"__pyx_fuse_{j}{routine_name}"
        if fused_name in d:
            has_fused = True
            d[_DeprecationHelperStr(fused_name, depdoc)] = d.pop(fused_name)
            j += 1
        else:
            break

    # If not, apply deprecation to the named routine
    if not has_fused:
        d[_DeprecationHelperStr(routine_name, depdoc)] = d.pop(routine_name)

