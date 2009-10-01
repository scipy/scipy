Special functions (:mod:`scipy.special`)
========================================

.. sectionauthor:: Travis E. Oliphant

.. currentmodule:: scipy.special

The main feature of the :mod:`scipy.special` package is the definition of
numerous special functions of mathematical physics. Available
functions include airy, elliptic, bessel, gamma, beta, hypergeometric,
parabolic cylinder, mathieu, spheroidal wave, struve, and
kelvin. There are also some low-level stats functions that are not
intended for general use as an easier interface to these functions is
provided by the ``stats`` module. Most of these functions can take
array arguments and return array results following the same
broadcasting rules as other math functions in Numerical Python. Many
of these functions also accept complex numbers as input. For a
complete list of the available functions with a one-line description
type ``>>> help(special).`` Each function also has its own
documentation accessible using help.  If you don't see a function you
need, consider writing it and contributing it to the library. You can
write the function in either C, Fortran, or Python. Look in the source
code of the library for examples of each of these kinds of functions.
