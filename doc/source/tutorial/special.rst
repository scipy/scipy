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
------------------------------------------------------------

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

.. _special-orthogonal-polynomials:
   
Orthogonal Polynomials
----------------------

A family of polynomials :math:`p_n` defined on an interval :math:`(a,
b)` is orthogonal if

.. math::

   \int_a^b p_m(x)p_n(x)w(x)dx = 0

for :math:`m \ne n`. The function :math:`w(x) \geq 0` is called the
weight function. Orthogonal polynomials are primarily used for
quadrature and approximating functions. Users interested in
approximating functions should investigate `numpy.polynomial`. Scipy
provides functions for evaluating orthogonal polynomials as well as
functions related to using them in quadrature.

Orthogonal polynomials are useful in quadrature because of the
following result: let :math:`x_k` be the roots of an orthogonal
polynomial :math:`p_n`, and define weights

.. math::

   w_k = \int_a^b \frac{p_n(x)}{(x - x_k)p_n'(x_k)}w(x)dx.

Then the quadrature rule

.. math::

   \int_a^b f(x)w(x)dx \approx \sum_{k = 1}^n w_kf(x_k)

integrates polynomials up to degree :math:`2n - 1` exactly. Because of
this result it is of interest to calculate the roots of orthogonal
polynomials and the weights associated with them. The interested
reader can consult section 3.5(v) in [dlmf]_ for more information.

The following table summarizes the properties of the orthogonal
polynomials supported by Scipy. It is largely reproduced from [as]_.

==================================== ============================= ========================= =================================== ==========================
Name                                 Notation                      Interval                  Weight                              Remarks
==================================== ============================= ========================= =================================== ==========================
Jacobi                               :math:`P_n^{(\alpha, \beta)}` :math:`(-1, 1)`           :math:`(1 - x)^\alpha(1 + x)^\beta` :math:`\alpha, \beta > -1`
Shifted Jacobi                       :math:`G_n(p, q, x)`          :math:`(0, 1)`            :math:`(1 - x)^{p - q}x^{q - 1}`    :math:`p - q > -1, q > 0`
Ultraspherical (Gegenbauer)          :math:`C_n^{(\alpha)}(x)`     :math:`(-1, 1)`           :math:`(1 - x^2)^{\alpha - 1/2}`    :math:`\alpha > -1/2`
Chebyshev of the first kind          :math:`T_n(x)`                :math:`(-1, 1)`           :math:`(1 - x^2)^{-1/2}`
Chebyshev of the second kind         :math:`U_n(x)`                :math:`(-1, 1)`           :math:`(1 - x^2)^{1/2}`
Chebyshev of the first kind          :math:`C_n(x)`                :math:`(-2, 2)`           :math:`(1 - x^2/4)^{-1/2}`
Chebyshev of the second kind         :math:`S_n(x)`                :math:`(-2, 2)`           :math:`(1 - x^2/4)^{1/2}`
Shifted Chebyshev of the first kind  :math:`T_n^*(x)`              :math:`(0, 1)`            :math:`(x - x^2)^{-1/2}`
Shifted Chebyshev of the second kind :math:`U_n^*(x)`              :math:`(0, 1)`            :math:`(x - x^2)^{1/2}`
Legendre (Spherical)                 :math:`P_n(x)`                :math:`(-1, 1)`           :math:`1`
Shifted Legendre                     :math:`P_n^*(x)`              :math:`(0, 1)`            :math:`1`
Generalized Laguerre                 :math:`L_n^{(\alpha)}(x)`     :math:`(0, \infty)`       :math:`e^{-x}x^\alpha`              :math:`\alpha > -1`
Laguerre                             :math:`L_n(x)`                :math:`(0, \infty)`       :math:`e^{-x}`
Hermite (Physicist's)                :math:`H_n(x)`                :math:`(-\infty, \infty)` :math:`e^{-x^2}`
Hermite (Probabilist's)              :math:`He_n(x)`               :math:`(-\infty, \infty)` :math:`e^{-x^2/2}`
==================================== ============================= ========================= =================================== ==========================

The following functions are available for evaluating orthogonal
polynomials, getting their roots and associated weights, and
evaluating their weight functions.

==================================== ================== =================== ====================
Name                                 Evaluation         Roots and weights   Weight function
==================================== ================== =================== ====================
Jacobi                               `eval_jacobi`      `roots_jacobi`      `weight_jacobi`
Shifted Jacobi                       `eval_sh_jacobi`   `roots_sh_jacobi`   `weight_sh_jacobi`
Ultraspherical (Gegenbauer)          `eval_gegenbauer`  `roots_gegenbauer`  `weight_gegenbauer`
Chebyshev of the first kind          `eval_chebyt`      `roots_chebyt`      `weight_chebyt`
Chebyshev of the second kind         `eval_chebyu`      `roots_chebyu`      `weight_chebyu`
Chebyshev of the first kind          `eval_chebyc`      `roots_chebyc`      `weight_chebyc`
Chebyshev of the second kind         `eval_chebys`      `roots_chebys`      `weight_chebys`
Shifted Chebyshev of the first kind  `eval_sh_chebyt`   `roots_sh_chebyt`   `weight_sh_chebyt`
Shifted Chebyshev of the second kind `eval_sh_chebyu`   `roots_sh_chebyu`   `weight_sh_chebyu`
Legendre (Spherical)                 `eval_legendre`    `roots_legendre`    `weight_legendre`
Shifted Legendre                     `eval_sh_legendre` `roots_sh_legendre` `weight_sh_legendre`
Generalized Laguerre                 `eval_genlaguerre` `roots_genlaguerre` `weight_genlaguerre`
Laguerre                             `eval_laguerre`    `roots_laguerre`    `weight_laguerre`
Hermite (Physicist's)                `eval_hermite`     `roots_hermite`     `weight_hermite`
Hermite (Probabilist's)              `eval_hermitenorm` `roots_hermitenorm` `weight_hermitenorm`
==================================== ================== =================== ====================

As an example of using these functions, suppose that we wished to
compute the Laplace transform of :math:`\sin(x)`:

.. math::

   \int_0^\infty e^{-sx}\sin(x)dx

at :math:`s = 1`. This looks like integrating :math:`\sin(x)` against
the Laguerre weight function :math:`w(x) = e^{-x}`, so we do:
   
.. code::

   >>> import scipy.special as sc
   >>> x, w = sc.roots_laguerre(20)
   >>> np.sum(np.sin(x)*w)
   >>> np.sum(np.sin(x)*w)
   0.49999999999998251

The exact answer is :math:`1/2`.

References
----------

.. [dlmf] NIST, "Digital Library of Mathematical Functions",
	  http://dlmf.nist.gov/
.. [as] Abramowitz, Stegun, "Handbook of Mathematical Functions",
	National Bureau of Standards, 1964.
