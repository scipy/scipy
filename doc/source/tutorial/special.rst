Special functions (:mod:`scipy.special`)
========================================

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

Bessel functions of real order(:func:`jn`, :func:`jn_zeros`)
--------------------------------------------------------------------------
Bessel functions are a family of solutions to Bessel's differential equation
with real or complex order alpha:

.. math::
   x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} + (x^2 - \alpha^2)y = 0

Among other uses, these functions arise in wave propagation problems such as
the vibrational modes of a thin drum head.  Here is an example of a circular
drum head anchored at the edge:

.. plot::

   >>> from scipy import special
   >>> def drumhead_height(n, k, distance, angle, t):
   ...    kth_zero = special.jn_zeros(n, k)[-1]
   ...    return np.cos(t) * np.cos(n*angle) * special.jn(n, distance*kth_zero)
   >>> theta = np.r_[0:2*np.pi:50j]
   >>> radius = np.r_[0:1:50j]
   >>> x = np.array([r * np.cos(theta) for r in radius])
   >>> y = np.array([r * np.sin(theta) for r in radius])
   >>> z = np.array([drumhead_height(1, 1, r, theta, 0.5) for r in radius])

   >>> import matplotlib.pyplot as plt
   >>> from mpl_toolkits.mplot3d import Axes3D
   >>> from matplotlib import cm
   >>> fig = plt.figure()
   >>> ax = Axes3D(fig)
   >>> ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
   >>> ax.set_xlabel('X')
   >>> ax.set_ylabel('Y')
   >>> ax.set_zlabel('Z')
   >>> plt.show()

..   :caption: Vibrating drum head using
..             :obj:`scipy.special.jn`
