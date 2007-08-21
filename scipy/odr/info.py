"""Orthogonal Distance Regression

Introduction
============

Why Orthogonal Distance Regression (ODR)? Sometimes one has measurement errors
in the explanatory variable, not just the response variable. Ordinary Least
Squares (OLS) fitting procedures treat the data for explanatory variables as
fixed. Furthermore, OLS procedures require that the response variable be an
explicit function of the explanatory variables; sometimes making the equation
explicit is unwieldy and introduces errors. ODR can handle both of these cases
with ease and can even reduce to the OLS case if necessary.

ODRPACK is a FORTRAN-77 library for performing ODR with possibly non-linear
fitting functions. It uses a modified trust-region Levenberg-Marquardt-type
algorithm to estimate the function parameters. The fitting functions are
provided by Python functions operating on NumPy arrays. The required derivatives
may be provided by Python functions as well or may be numerically estimated.
ODRPACK can do explicit or implicit ODR fits or can do OLS. Input and output
variables may be multi-dimensional. Weights can be provided to account for
different variances of the observations (even covariances between dimensions of
the variables).

odr provides two interfaces: a single function and a set of high-level classes
that wrap that function. Please refer to their docstrings for more information.
While the docstring of the function, odr, does not have a full explanation of
its arguments, the classes do, and the arguments with the same name usually have
the same requirements. Furthermore, it is highly suggested that one at least
skim the ODRPACK User's Guide.  Know Thy Algorithm.


Use
===

See the docstrings of odr.odrpack and the functions and classes for
usage instructions. The ODRPACK User's Guide is also quite helpful. It can be
found on one of the ODRPACK's original author's website:

    http://www.boulder.nist.gov/mcsd/Staff/JRogers/odrpack.html

Robert Kern
robert.kern@gmail.com
"""
postpone_import = 1
