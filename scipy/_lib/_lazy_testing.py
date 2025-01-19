from collections.abc import Callable, Iterable, Sequence
from types import ModuleType
import pytest
from scipy._lib._array_api import is_jax


def lazy_xp_function(
    func: Callable,
    *,
    jax_jit: bool = True,
    static_argnums: int | Sequence[int] | None=None,
    static_argnames: str | Iterable[str] | None=None,
) -> None:
    """Tag a function, which must be imported in the test module globals,
    so that when any tests defined in the same module are executed with
    xp=jax.numpy the function is replaced with a jitted version of itself.

    This will be later expanded to provide test coverage for other lazy backends,
    e.g. Dask.

    Example::

      # test_mymodule.py:
      from scipy._lib._lazy_testing import lazy_xp_function
      from scipy.mymodule import myfunc

      lazy_xp_function(myfunc)

      def test_myfunc(xp):
          a = xp.asarray([1, 2])
          # When xp=jax.numpy, this is the same as
          # b = jax.jit(myfunc)(a)
          b = myfunc(a)

    Parameters
    ----------
    func : callable
        Function to be tested
    jax_jit : bool, optional
        Set to True to replace `func` with `jax.jit(func)` when calling the
        `patch_lazy_xp_functions` test helper with `xp=jax.numpy`.
        Set to False if `func` is only compatible with eager (non-jitted) JAX.
        Default: True.
    static_argnums : int | Sequence[int], optional
        Passed to jax.jit.
        Positional arguments to treat as static (trace- and compile-time constant).
        Default: infer from static_argnames using `inspect.signature(func)`.
    static_argnames : str | Iterable[str], optional
        Passed to jax.jit.
        Named arguments to treat as static (compile-time constant).
        Default: infer from static_argnums using `inspect.signature(func)`.

    Notes
    -----
    A test function can circumvent this monkey-patching system by calling `func` an
    attribute of the original module. You need to sanitize your code to
    make sure this does not happen.

    Example::

      import mymodule
      from mymodule import myfunc

      lazy_xp_function(myfunc)

      def test_myfunc(xp):
          a = xp.asarray([1, 2])
          b = myfunc(a)  # This is jitted when xp=jax.numpy
          c = mymodule.myfunc(a)  # This is not

    See Also
    --------
    patch_lazy_xp_functions
    jax.jit: https://jax.readthedocs.io/en/latest/_autosummary/jax.jit.html
    """
    if jax_jit:
        func._lazy_jax_jit_kwargs = {  # type: ignore[attr-defined]
            "static_argnums": static_argnums,
            "static_argnames": static_argnames,
        }


def patch_lazy_xp_functions(
    xp: ModuleType, request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch
) -> None:
    """If xp==jax.numpy, search for all functions which have been tagged by
    `lazy_xp_function` in the globals of the module that defines the current test
    and wrap them with `jax.jit`. Unwrap them at the end of the test.

    Parameters
    ----------
    xp: module
        Array namespace to be tested
    request: pytest.FixtureRequest
        Pytest fixture, as acquired by the test itself or by one of its fixtures.
    monkeypatch: pytest.MonkeyPatch
        Pytest fixture, as acquired by the test itself or by one of its fixtures.

    See Also
    --------
    lazy_xp_function
    https://docs.pytest.org/en/stable/reference/reference.html#std-fixture-request
    """
    if is_jax(xp):
        import jax

        globals_ = request.module.__dict__
        for name, func in globals_.items():
            kwargs = getattr(func, "_lazy_jax_jit_kwargs", None)
            if kwargs is not None:
                monkeypatch.setitem(globals_, name, jax.jit(func, **kwargs))
