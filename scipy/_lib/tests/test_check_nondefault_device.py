"""Unit tests for ``tools/check_nondefault_device.py``.

The checker lives in the ``tools`` directory of a source checkout.  Tests are
run from an installed copy of SciPy (where ``tools`` is absent), so the checker
is located relative to the pytest root directory (the repo root) and the tests
are skipped when it cannot be found.
"""
import importlib.util

import pytest


@pytest.fixture(scope="module")
def check_source(pytestconfig):
    tool = pytestconfig.rootpath / "tools" / "check_nondefault_device.py"
    if not tool.is_file():
        pytest.skip("check_nondefault_device.py requires a source checkout")
    spec = importlib.util.spec_from_file_location("check_nondefault_device", tool)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.check_source


@pytest.mark.parametrize("src", [
    "xp.zeros(n)",
    "xp.ones(n)",
    "xp.empty(n)",
    "xp.full(s, 0)",
    "xp.arange(n)",
    "xp.eye(n)",
    "xp.linspace(0, 1, n)",
    "xp.asarray(0)",
    "xp.asarray(-1000)",
    "xp.asarray([1, 2, 3])",
    "xp.asarray((1, 2))",
    "xp.fft.rfftfreq(n)",
    "xp.fft.fftfreq(n)",
    # nested creation is still flagged
    "xp.reshape(xp.arange(k), shape)",
])
def test_flagged(check_source, src):
    assert len(check_source(src)) == 1


@pytest.mark.parametrize("src", [
    # device is propagated explicitly
    "xp.zeros(n, device=xp_device(x))",
    "xp.full(s, 0, dtype=x.dtype, device=xp_device(x))",
    "xp.arange(n, device=dev)",
    # asarray of an existing array/expression infers the device
    "xp.asarray(x)",
    "xp.asarray(a.dtype)",
    "xp.asarray(some_func())",
    # *_like infers the device from its reference array
    "xp.zeros_like(x)",
    "xp.full_like(x, 0)",
    "xp.empty_like(x)",
    # **kwargs splat might supply device -> lenient, not flagged
    "xp.zeros(n, **kwargs)",
    # not the `xp` namespace
    "np.zeros(n)",
    "self._xp.zeros(n)",
    "mxp.zeros(n)",
    # not a creation function
    "xp.sum(x)",
])
def test_not_flagged(check_source, src):
    assert len(check_source(src)) == 0


@pytest.mark.parametrize("src", [
    "xp.zeros(n)  # skip device check",
    "xp.asarray([1, 2])  # skip device check",
])
def test_pragma_suppresses(check_source, src):
    assert len(check_source(src)) == 0


def test_multiple_violations_reported(check_source):
    src = "a = xp.zeros(n)\nb = xp.ones(n)\nc = xp.asarray(x)\n"
    violations = check_source(src)
    assert len(violations) == 2
    # line numbers are 1-based and point at the offending lines
    assert sorted(lineno for lineno, _, _ in violations) == [1, 2]


def test_syntax_error_is_ignored(check_source):
    assert check_source("def f(:\n") == []
