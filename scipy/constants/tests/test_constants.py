import pytest

import scipy.constants as sc
from scipy.conftest import array_api_compatible
from scipy._lib._array_api import xp_assert_equal, xp_assert_close


@array_api_compatible
def test_convert_temperature(xp):
    xp_assert_equal(sc.convert_temperature(xp.asarray(32.), 'f', 'Celsius'),
                    xp.asarray(0.0))
    xp_assert_equal(sc.convert_temperature(xp.asarray([0., 0.]), 'celsius', 'Kelvin'),
                    xp.asarray([273.15, 273.15]))
    xp_assert_equal(sc.convert_temperature(xp.asarray([0., 0.]), 'kelvin', 'c'),
                    xp.asarray([-273.15, -273.15]))
    xp_assert_equal(sc.convert_temperature(xp.asarray([32., 32.]), 'f', 'k'),
                    xp.asarray([273.15, 273.15]))
    xp_assert_equal(sc.convert_temperature(xp.asarray([273.15, 273.15]), 'kelvin', 'F'),
                    xp.asarray([32., 32.]))
    xp_assert_equal(sc.convert_temperature(xp.asarray([0., 0.]), 'C', 'fahrenheit'),
                    xp.asarray([32., 32.]))
    xp_assert_close(sc.convert_temperature(xp.asarray([0., 0.]), 'c', 'r'),
                    xp.asarray([491.67, 491.67]),
                    rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature(xp.asarray([491.67, 491.67]),
                                           'Rankine', 'C'),
                    xp.asarray([0., 0.]), rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature(xp.asarray([491.67, 491.67]), 'r', 'F'),
                    xp.asarray([32., 32.]), rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature(xp.asarray([32., 32.]), 'fahrenheit', 'R'),
                    xp.asarray([491.67, 491.67]), rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature(xp.asarray([273.15, 273.15]), 'K', 'R'),
                    xp.asarray([491.67, 491.67]), rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature(xp.asarray([491.67, 0.]),
                                           'rankine', 'kelvin'),
                    xp.asarray([273.15, 0.]), rtol=0., atol=1e-13)


@array_api_compatible
def test_convert_temperature_errors(xp):
    with pytest.raises(NotImplementedError, match="old_scale="):
        sc.convert_temperature(1, old_scale="cheddar", new_scale="kelvin")
    with pytest.raises(NotImplementedError, match="new_scale="):
        sc.convert_temperature(1, old_scale="kelvin", new_scale="brie")


@array_api_compatible
def test_lambda_to_nu(xp):
    xp_assert_equal(sc.lambda2nu(xp.asarray([sc.speed_of_light, 1])),
                    xp.asarray([1, sc.speed_of_light]))


@array_api_compatible
def test_nu_to_lambda(xp):
    xp_assert_equal(sc.nu2lambda(xp.asarray([sc.speed_of_light, 1])),
                    xp.asarray([1, sc.speed_of_light]))

