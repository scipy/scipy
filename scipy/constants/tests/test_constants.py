from __future__ import division, print_function, absolute_import

import warnings

from numpy.testing import run_module_suite, assert_equal, assert_allclose
import scipy.constants as sc


def test_convert_temperature():
    assert_equal(sc.convert_temperature(32, 'f', 'Celsius'), 0)
    assert_equal(sc.convert_temperature([0, 0], 'celsius', 'Kelvin'),
                 [273.15, 273.15])
    assert_equal(sc.convert_temperature([0, 0], 'kelvin', 'c'),
                 [-273.15, -273.15])
    assert_equal(sc.convert_temperature([32, 32], 'f', 'k'), [273.15, 273.15])
    assert_equal(sc.convert_temperature([273.15, 273.15], 'kelvin', 'F'),
                 [32, 32])
    assert_equal(sc.convert_temperature([0, 0], 'C', 'fahrenheit'), [32, 32])
    assert_allclose(sc.convert_temperature([0, 0], 'c', 'r'), [491.67, 491.67],
                    rtol=0., atol=1e-13)
    assert_allclose(sc.convert_temperature([491.67, 491.67], 'Rankine', 'C'),
                    [0., 0.], rtol=0., atol=1e-13)
    assert_allclose(sc.convert_temperature([491.67, 491.67], 'r', 'F'),
                    [32., 32.], rtol=0., atol=1e-13)
    assert_allclose(sc.convert_temperature([32, 32], 'fahrenheit', 'R'),
                    [491.67, 491.67], rtol=0., atol=1e-13)
    assert_allclose(sc.convert_temperature([273.15, 273.15], 'K', 'R'),
                    [491.67, 491.67], rtol=0., atol=1e-13)
    assert_allclose(sc.convert_temperature([491.67, 0.], 'rankine', 'kelvin'),
                    [273.15, 0.], rtol=0., atol=1e-13)


def test_fahrenheit_to_celcius():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert_equal(sc.F2C(32), 0)
        assert_equal(sc.F2C([32, 32]), [0, 0])


def test_celcius_to_kelvin():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert_equal(sc.C2K([0, 0]), [273.15, 273.15])


def test_kelvin_to_celcius():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert_equal(sc.K2C([0, 0]), [-273.15, -273.15])


def test_fahrenheit_to_kelvin():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert_equal(sc.F2K([32, 32]), [273.15, 273.15])


def test_kelvin_to_fahrenheit():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert_equal(sc.K2F([273.15, 273.15]), [32, 32])


def test_celcius_to_fahrenheit():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert_equal(sc.C2F([0, 0]), [32, 32])


def test_lambda_to_nu():
    assert_equal(sc.lambda2nu(sc.speed_of_light), 1)


def test_nu_to_lambda():
    assert_equal(sc.nu2lambda(1), sc.speed_of_light)

def test_ft_to_meters():
    assert_equal(sc.ft_to_meters([1,1]),[0.3048,0.3048])

def test_meters_to_ft():
    assert_equal(sc.meters_to_ft([1,1]), [1/0.3048, 1/0.30048])

def test_inches_to_ft():
    assert_equal(sc.inches_to_ft([1,1]),[1/12,1/12])

def test_ft_to_inches():
    assert_equal(sc.ft_to_inches([1,1]),[12,12])

def test_in_to_meters():
    assert_equal(sc.in_to_meters([1,1]),[0.0254,0.0254])

def test_meters_to_in():
    assert_equal(sc.meters_to_in([1,1]),[1/0.0254,1/0.0254])

def test_lb_to_kg():
    assert_equal(sc.lb_to_kg([1/0.453592, 0.453592]),[1,1])

def test_kg_to_lb():
    assert_equal(sc.kg_to_lb([1,1]),[0.453592,0.453592])

def test_oz_to_g():
    assert_equal(sc.oz_to_g([28.3495,28.3495]),[1,1])

def test_g_to_oz():
    assert_equal(sc.g_to_oz([1,1]),[28.3495,28.3495])

def test_sec_to_min():
    assert_equal(sc.sec_to_min([1,1]),[1/60,1/60])

def test_min_to_sec():
    assert_equal(sc.min_to_sec([1,1]),[60,60])

def test_sec_to_h():
    assert_equal(sc.sec_to_h([1,1]), [1/3600,1/3600])

def test_h_to_sec():
    assert_equal(sc.h_to_sec([1,1]),[3600,3600])

def test_sec_to_day():
    assert_equal(sc.sec_to_day([1,1]), [1/24*3600, 1/24*3600])

def test_day_to_sec():
    assert_equal(sc.day_to_sec([1,1]),[24*3600, 24*3600])

def test_sec_to_week():
    assert_equal(sc.sec_to_week([1,1]), [1/7*24*3600,1/7*24*3600])

def test_week_to_sec():
    assert_equal(sc.week_to_sec([1,1]),[7*24*3600,7*24*3600])

def test_sec_to_year():
    assert_equal(sc.sec_to_year([1,1]),[1/365*24*3600,1/365*24*3600])

def test_year_to_sec():
    assert_equal(sc.year_to_sec([1,1]),[365*24*3600,365*24*3600])

def test_min_to_h():
    assert_equal(sc.min_to_h([1,1]),[1/60,1/60])

def test_h_to_min():
    assert_equal(sc.h_to_min([1,1]).[60,60])

def test_min_to_day():
    assert_equal(sc.min_to_day([1,1]),[1/24*60,24*60])

def test_day_to_min():
    assert_equal(sc.day_to_min([1,1]),[24*60,24*60])

def test_min_to_week():
    assert_equal(sc.min_to_week([1,1]).[1/7*24*60,1/7*24*60])

def test_week_to_min():
    assert_equal(sc.week_to_min([1,1]),[7*24*60,7*24*60])

def test_min_to_year():
    assert_equal(sc.min_to_year([1,1]),[1/365*24*60, 1/365*24*60])

def test_year_to_min():
    assert_equal(sc.year_to_min([1,1]),[365*24*60,365*24*60])

def test_h_to_day():
    assert_equal(sc.h_to_day([1,1]),[1/24,1/24])

def test_day_to_h():
    assert_equal(sc.day_to_h([1,1]),[24,24])

def test_h_to_week():
    assert_equal(sc.h_to_week([1,1]),[1/24*7,1/24*7])

def test_week_to_h():
    assert_equal(sc.week_to_h([1,1]),[24*7,24*7])

def test_h_to_year():
    assert_equal(sc.h_to_year([1,1]),[1/365*24,1/365*24])

def test_year_to_h():
    assert_equal(sc.year_to_h([1,1]),[365*24, 365*24])

def test_day_to_week():
    assert_equal(sc.day_to_week([1,1]), [1/7,1/7])

def test_week_to_day():
    assert_equal(sc.week_to_day([1,1]),[7,7])


if __name__ == "__main__":
    run_module_suite()
