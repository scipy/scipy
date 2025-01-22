import pytest
from scipy._lib._array_api import array_namespace, is_jax, xp_assert_equal
from scipy._lib._lazy_testing import lazy_xp_function


def jittable(x):
    """A jittable function"""
    return x * 2.0


def non_jittable(x):
    """This function materializes the input array, so it will fail
    when wrapped in jax.jit
    """
    xp = array_namespace(x)
    if xp.any(x < 0.0):
        raise ValueError("Negative values not allowed")
    return x


def non_jittable2(x):
    return non_jittable(x)


def static_params(x, n, flag=False):
    """Function with static parameters that must not be jitted"""
    if flag and n > 0:  # This fails if n or flag are jitted arrays
        return x * 2.0
    else:
        return x * 3.0


def static_params1(x, n, flag=False):
    return static_params(x, n, flag)


def static_params2(x, n, flag=False):
    return static_params(x, n, flag)


def static_params3(x, n, flag=False):
    return static_params(x, n, flag)


lazy_xp_function(jittable)
lazy_xp_function(non_jittable2)
lazy_xp_function(static_params1, static_argnums=(1, 2))
lazy_xp_function(static_params2, static_argnames=("n", "flag"))
lazy_xp_function(static_params3, static_argnums=1, static_argnames="flag")


def test_lazy_xp_function(xp):
    x = xp.asarray([1.0, 2.0])

    xp_assert_equal(jittable(x), xp.asarray([2.0, 4.0]))

    xp_assert_equal(non_jittable(x), xp.asarray([1.0, 2.0]))  # Not jitted
    if is_jax(xp):
        with pytest.raises(
            TypeError, match="Attempted boolean conversion of traced array"
        ):
            non_jittable2(x)  # Jitted
    else:
        xp_assert_equal(non_jittable2(x), xp.asarray([1.0, 2.0]))


@pytest.mark.parametrize("func", [static_params1, static_params2, static_params3])
def test_lazy_xp_function_static_params(xp, func):
    x = xp.asarray([1.0, 2.0])
    xp_assert_equal(func(x, 1), xp.asarray([3.0, 6.0]))
    xp_assert_equal(func(x, 1, True), xp.asarray([2.0, 4.0]))
    xp_assert_equal(func(x, 1, False), xp.asarray([3.0, 6.0]))
    xp_assert_equal(func(x, 0, False), xp.asarray([3.0, 6.0]))
    xp_assert_equal(func(x, 1, flag=True), xp.asarray([2.0, 4.0]))
    xp_assert_equal(func(x, n=1, flag=True), xp.asarray([2.0, 4.0]))
