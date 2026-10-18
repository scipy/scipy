"""Pytest plugin that shows the SciPy test environment variables in the header.

It is registered in the repository's ``pytest.ini`` via
``-p scipy._lib._pytest_header`` so that it is loaded before the report header
is printed, both for path-style invocations (``spin test scipy/...``) and for
``--pyargs`` ones (``spin test -t scipy.xxx``). A ``pytest_report_header`` hook
in ``scipy/conftest.py`` would only be seen for the former, because with
``--pyargs`` that conftest is loaded during collection, after the header.
"""
import os

# Environment variables that change how the test suite runs.
TEST_ENV_VARS = (
    "SCIPY_ARRAY_API",
    "SCIPY_DEVICE",
    "SCIPY_DEFAULT_DTYPE",
    "SCIPY_XSLOW",
    "SCIPY_HYPOTHESIS_PROFILE",
    "SCIPY_AVAILABLE_MEM",
)


def environment_report_header():
    """Return a one-line summary of the SciPy test environment variables."""
    values = ", ".join(f"{name}={os.environ.get(name, '<unset>')}"
                       for name in TEST_ENV_VARS)
    return f"SciPy test environment: {values}"


def pytest_report_header(config):
    return environment_report_header()
