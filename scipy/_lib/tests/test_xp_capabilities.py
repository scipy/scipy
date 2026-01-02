import pytest

from scipy._lib.array_api_compat import numpy as np
from scipy._lib._array_api import (
    array_namespace, make_xp_pytest_param, xp_capabilities, _xp_copy_to_numpy,
    xp_assert_close, make_xp_test_case
)

local_capabilities_table = {}

# B is a child of A which inherits the method g which is array-agnostic
# so long as the method f is supported. A.f does not support the JAX jit but
# B.f does support the JAX jit. Test that this inheritence does not
# cause problems when testing with JAX jit. For verisimilitude, there is also
# a method h such that A.h and B.h both support that jit, making it so it
# makes sense to have both A and B in lazy_xp_modules.

@xp_capabilities(
    capabilities_table=local_capabilities_table,
    cpu_only=True,
    skip_backends=[("dask.array", ""), ("torch", "")],
    jax_jit=False,
    method_capabilities={
        "h": dict(jax_jit=True, cpu_only=True,
                  skip_backends=[("dask.array", ""), ("torch", "")]),
    },
)
class A:
    def __init__(self, x):
        xp = array_namespace(x)
        self._xp = xp
        self.x = np.asarray(x)

    def f(self, y):
        y = np.asarray(y)
        return self._xp.asarray(np.matmul(self.x, y))

    def g(self, y, z):
        return self.f(y) + self.f(z)

    def h(self, y):
        xp = array_namespace(y)
        y = xp.asarray(y)
        return xp.matmul(y, y)


@xp_capabilities(
    capabilities_table=local_capabilities_table,
    skip_backends=[("cupy", ""), ("dask.array", ""), ("torch", "")],
)
class B(A):
    def __init__(self, x):
        xp = array_namespace(x)
        self._xp = xp
        self.x = xp.asarray(x)

    def f(self, y):
        return self._xp.matmul(self.x, y)


lazy_xp_modules = [A, B]


@pytest.mark.parametrize(
    "cls",
    [
        make_xp_pytest_param(
            (A, "g"),
            capabilities_table=local_capabilities_table,
            additional_marks=pytest.mark.xfail(reason="spooky action at a distance"),
        ),
        make_xp_pytest_param((B, "g"), capabilities_table=local_capabilities_table),
    ],
)
def test_inheritance1(cls, xp):
    x = xp.asarray([1.1, 2.2, 3.3])
    y = xp.asarray([1.0, 2.0, 3.0])
    z = xp.asarray([3.0, 4.0, 5.0])
    foo = cls(x)
    observed = foo.g(y, z)
    expected = xp.asarray(44.0)[()]
    xp_assert_close(observed, expected)


@pytest.mark.parametrize(
    "cls",
    [
        make_xp_pytest_param((A, "h"), capabilities_table=local_capabilities_table),
        make_xp_pytest_param((B, "h"), capabilities_table=local_capabilities_table),
    ],
)
def test_inheritance2(cls, xp):
    x = xp.asarray([1.1, 2.2, 3.3])
    y = xp.asarray([1.0, 2.0, 3.0])
    foo = cls(x)
    observed = foo.h(y)
    expected = xp.asarray(14.0)[()]
    xp_assert_close(observed, expected)
