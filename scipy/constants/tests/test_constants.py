import scipy.constants as sc
from scipy.conftest import array_api_compatible
from scipy._lib._array_api import xp_assert_equal, xp_assert_close


@array_api_compatible
def test_convert_temperature(xp):
    xp_assert_equal(sc.convert_temperature(32, 'f', 'Celsius'), 0.0)
    xp_assert_equal(sc.convert_temperature([0, 0], 'celsius', 'Kelvin'),
                    [273.15, 273.15])
    xp_assert_equal(sc.convert_temperature([0, 0], 'kelvin', 'c'),
                    [-273.15, -273.15])
    xp_assert_equal(sc.convert_temperature([32, 32], 'f', 'k'), [273.15, 273.15])
    xp_assert_equal(sc.convert_temperature([273.15, 273.15], 'kelvin', 'F'),
                    [32.0, 32.0])
    xp_assert_equal(sc.convert_temperature([0, 0], 'C', 'fahrenheit'), [32.0, 32.0])
    xp_assert_close(sc.convert_temperature([0, 0], 'c', 'r'), [491.67, 491.67],
                    rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature([491.67, 491.67], 'Rankine', 'C'),
                    [0., 0.], rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature([491.67, 491.67], 'r', 'F'),
                    [32., 32.], rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature([32, 32], 'fahrenheit', 'R'),
                    [491.67, 491.67], rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature([273.15, 273.15], 'K', 'R'),
                    [491.67, 491.67], rtol=0., atol=1e-13)
    xp_assert_close(sc.convert_temperature([491.67, 0.], 'rankine', 'kelvin'),
                    [273.15, 0.], rtol=0., atol=1e-13)


@array_api_compatible
def test_lambda_to_nu(xp):
    xp_assert_equal(sc.lambda2nu([sc.speed_of_light, 1]), [1, sc.speed_of_light])


@array_api_compatible
def test_nu_to_lambda(xp):
    xp_assert_equal(sc.nu2lambda([sc.speed_of_light, 1]), [1, sc.speed_of_light])

