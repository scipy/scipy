import sys
import uarray
import pytest  # type: ignore

from .tests import example_helpers


@pytest.fixture(autouse=True)
def add_namespaces(doctest_namespace):
    doctest_namespace["ua"] = uarray
    doctest_namespace["ex"] = example_helpers
