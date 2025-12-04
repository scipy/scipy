"""Check for functions advertising alt backend support without tests.

This checks for functions/classes/methods ``f`` which have been decorated with
`xp_capabilities` advertising something other than `np_only=True`, but which
have no tests using the `xp` fixture that have been registered as testing ``f``
with one of `make_xp_test_case`, `make_xp_pytest_param`, or
`make_xp_pytest_marks`. The intention is that for a function to be documented
as supported on alternative backends, this capability needs to be tested. Functions
``f`` which are tested but whose tests have not yet been converted to relying
on `xp_capabilities` plus `make_xp_test_case` or one of its relatives can be
accounted for by explicitly adding the marker
`pytest.mark.uses_xp_capabilities(False, funcs=[f])` to the associated tests.

`make_xp_test_case(f)` and its relatives add the marker
`pytest.mark.uses_xp_capabilities(True, funcs=[f])`.
"""

import contextlib
import io
import importlib
import os
import pytest
import sys

from pathlib import Path

import scipy

from scipy._lib._array_api import xp_capabilities_table, SCIPY_ARRAY_API
from scipy._lib._public_api import PUBLIC_MODULES


if __name__ == "__main__":

    if not SCIPY_ARRAY_API:
        sys.exit(0)

    # `xp_capabilities_table` is populated lazily, so need to import
    # everything to make sure all entries are there.
    for module_name in PUBLIC_MODULES:
        if module_name == "scipy.odr":
            continue
        importlib.import_module(module_name)

    tested_functions = set()

    class MarkerCollector:
        def pytest_collection_modifyitems(self, session, config, items):
            # `make_xp_pytest_marks(funcs)` adds a mark
            # `pytest.mark.uses_xp_capabilities` which stores funcs in a kwarg.
            # We can use this to identify which functions have associated tests
            # using the `xp` fixture.
            for item in items:
                marker = item.get_closest_marker("uses_xp_capabilities")
                if marker is None:
                    continue
                funcs = marker.kwargs.get("funcs")
                if funcs is None:
                    continue
                tested_functions.update(funcs)

    with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
        # Suppress output because it is extremely long and not needed.
        pytest.main(
            ["--collect-only", str(Path(scipy.__file__).parent)],
            plugins=[MarkerCollector()]
        )

    incorrectly_decorated = [
        func for func, data in xp_capabilities_table.items()
        if not data["np_only"]  and func not in tested_functions
    ]
    if incorrectly_decorated:
        print(
            "There are functions documented as supported on one or more"
            " alternative backends with no known tests using the `xp`"
            " fixture."
        )
        for func in incorrectly_decorated:
            print(f"- {func.__module__}.{func.__qualname__}")
        sys.exit(1)

    print(
        "All tests documented as supported on one or more alternative backends"
        " have tests using the `xp` fixture."
    )
    sys.exit(0)
