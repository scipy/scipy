from __future__ import division, print_function, absolute_import

from numpy.testing import run_module_suite, assert_equal
import scipy.constants as sc


def test_fahrenheit_to_celcius():
    assert_equal(sc.F2C(32), 0)
    assert_equal(sc.F2C([32, 32]), [0, 0])


def test_celcius_to_kelvin():
    assert_equal(sc.C2K([0, 0]), [273.15, 273.15])


def test_kelvin_to_celcius():
    assert_equal(sc.K2C([0, 0]), [-273.15, -273.15])


def test_fahrenheit_to_kelvin():
    assert_equal(sc.F2K([32, 32]), [273.15, 273.15])


def test_kelvin_to_fahrenheit():
    assert_equal(sc.K2F([273.15, 273.15]), [32, 32])


def test_celcius_to_fahrenheit():
    assert_equal(sc.C2F([0, 0]), [32, 32])


def test_lambda_to_nu():
    assert_equal(sc.lambda2nu(sc.speed_of_light), 1)


def test_nu_to_lambda():
    assert_equal(sc.nu2lambda(1), sc.speed_of_light)


if __name__ == "__main__":
    run_module_suite()
