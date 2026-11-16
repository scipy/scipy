"""Shims that let the SciPy test suite import, collect and run when the
*optional* ``hypothesis`` dependency is not installed.

``hypothesis`` supplements the deterministic pytest suite with property-based
tests. Since it grew a Rust extension it is hard to install on platforms
without prebuilt wheels, so it is no longer a required test dependency (see
gh-25837). Test modules import the names they need from here instead of from
``hypothesis`` directly.

When ``hypothesis`` is available the names are simply re-exported. When it is
missing, ``given``/``settings`` become identity decorators (so the decorated
test is still collected as a plain function) and the strategy objects become
no-op stand-ins (so evaluating decorator arguments does not raise at import
time). Each property-based test must also carry
``@pytest.mark.skipif(not HAS_HYPOTHESIS, reason="hypothesis is not installed")``
so that it is reported as skipped rather than run.

This mirrors ``numpy.testing._private.hypothesis_helpers``.
"""

try:
    import hypothesis
    from hypothesis import given, settings, strategies
    from hypothesis.extra import numpy as npst
    HAS_HYPOTHESIS = True
except ImportError:
    HAS_HYPOTHESIS = False

    def given(*args, **kwargs):
        # Identity decorator: the test stays a plain function so pytest
        # collects it and the accompanying skipif marker takes effect.
        return lambda func: func

    def settings(*args, **kwargs):
        return lambda func: func

    class _Swallow:
        """Callable / attribute / ``|``-chainable no-op standing in for the
        ``hypothesis`` strategy objects used in decorator arguments."""

        def __call__(self, *args, **kwargs):
            return self

        def __getattr__(self, name):
            return self

        def __or__(self, other):
            return self

        __ror__ = __or__

    _swallow = _Swallow()
    strategies = npst = _swallow

    class _HypothesisStub:
        """Stand-in for the top-level ``hypothesis`` module, for tests that
        use ``@hypothesis.given(...)`` / ``hypothesis.strategies`` directly."""

        given = staticmethod(given)
        settings = staticmethod(settings)
        strategies = _swallow

        def __getattr__(self, name):
            return _swallow

    hypothesis = _HypothesisStub()

# Alias used by some test modules (``import hypothesis.strategies as st``).
st = strategies

__all__ = ["HAS_HYPOTHESIS", "hypothesis", "given", "settings", "strategies",
           "st", "npst"]
