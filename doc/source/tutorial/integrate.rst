Integration (:mod:`scipy.integrate`)
====================================

.. sectionauthor:: Travis E. Oliphant

.. currentmodule:: scipy.integrate

The :mod:`scipy.integrate` sub-package provides several integration
techniques including an ordinary differential equation integrator. An
overview of the module is provided by the help command:

.. literalinclude:: examples/4-1


General integration (:func:`quad`)
----------------------------------

The function :obj:`quad` is provided to integrate a function of one
variable between two points. The points can be :math:`\pm\infty`
(:math:`\pm` ``inf``) to indicate infinite limits. For example,
suppose you wish to integrate a bessel function ``jv(2.5, x)`` along
the interval :math:`[0, 4.5].`

.. math::

    I=\int_{0}^{4.5}J_{2.5}\left(x\right)\, dx.


This could be computed using :obj:`quad`:

    >>> import scipy.integrate as integrate
    >>> import scipy.special as special
    >>> result = integrate.quad(lambda x: special.jv(2.5,x), 0, 4.5)
    >>> result
    (1.1178179380783249, 7.8663172481899801e-09)

    >>> from numpy import sqrt, sin, cos, pi
    >>> I = sqrt(2/pi)*(18.0/27*sqrt(2)*cos(4.5) - 4.0/27*sqrt(2)*sin(4.5) +
    ...                 sqrt(2*pi) * special.fresnel(3/sqrt(pi))[0])
    >>> I
    1.117817938088701

    >>> print(abs(result[0]-I))
    1.03761443881e-11

The first argument to quad is a "callable" Python object (i.e., a
function, method, or class instance). Notice the use of a lambda-
function in this case as the argument. The next two arguments are the
limits of integration. The return value is a tuple, with the first
element holding the estimated value of the integral and the second
element holding an upper bound on the error. Notice, that in this
case, the true value of this integral is

.. math::

    I=\sqrt{\frac{2}{\pi}}\left(\frac{18}{27}\sqrt{2}\cos\left(4.5\right)-\frac{4}{27}\sqrt{2}\sin\left(4.5\right)+\sqrt{2\pi}\textrm{Si}\left(\frac{3}{\sqrt{\pi}}\right)\right),

where

.. math::

    \textrm{Si}\left(x\right)=\int_{0}^{x}\sin\left(\frac{\pi}{2}t^{2}\right)\, dt.

is the Fresnel sine integral. Note that the numerically-computed integral is
within :math:`1.04\times10^{-11}` of the exact result --- well below the
reported error bound.


If the function to integrate takes additional parameters, they can be provided
in the `args` argument. Suppose that the following integral shall be calculated:

.. math::

    I(a,b)=\int_{0}^{1} ax^2+b \, dx.


This integral can be evaluated by using the following code:

>>> from scipy.integrate import quad
>>> def integrand(x, a, b):
...     return a*x**2 + b
...
>>> a = 2
>>> b = 1
>>> I = quad(integrand, 0, 1, args=(a,b))
>>> I
(1.6666666666666667, 1.8503717077085944e-14)


Infinite inputs are also allowed in :obj:`quad` by using :math:`\pm`
``inf`` as one of the arguments. For example, suppose that a numerical
value for the exponential integral:

.. math::

    E_{n}\left(x\right)=\int_{1}^{\infty}\frac{e^{-xt}}{t^{n}}\, dt.

is desired (and the fact that this integral can be computed as
``special.expn(n,x)`` is forgotten). The functionality of the function
:obj:`special.expn <scipy.special.expn>` can be replicated by defining a new function
``vec_expint`` based on the routine :obj:`quad`:

    >>> from scipy.integrate import quad
    >>> def integrand(t, n, x):
    ...     return np.exp(-x*t) / t**n
    ...

    >>> def expint(n, x):
    ...     return quad(integrand, 1, np.inf, args=(n, x))[0]
    ...

    >>> vec_expint = np.vectorize(expint)

    >>> vec_expint(3, np.arange(1.0, 4.0, 0.5))
    array([ 0.1097,  0.0567,  0.0301,  0.0163,  0.0089,  0.0049])
    >>> import scipy.special as special
    >>> special.expn(3, np.arange(1.0,4.0,0.5))
    array([ 0.1097,  0.0567,  0.0301,  0.0163,  0.0089,  0.0049])

The function which is integrated can even use the quad argument (though the
error bound may underestimate the error due to possible numerical error in the
integrand from the use of :obj:`quad` ). The integral in this case is

.. math::

    I_{n}=\int_{0}^{\infty}\int_{1}^{\infty}\frac{e^{-xt}}{t^{n}}\, dt\, dx=\frac{1}{n}.

>>> result = quad(lambda x: expint(3, x), 0, np.inf)
>>> print(result)
(0.33333333324560266, 2.8548934485373678e-09)

>>> I3 = 1.0/3.0
>>> print(I3)
0.333333333333

>>> print(I3 - result[0])
8.77306560731e-11

This last example shows that multiple integration can be handled using
repeated calls to :func:`quad`.


General multiple integration (:func:`dblquad`, :func:`tplquad`, :func:`nquad`)
------------------------------------------------------------------------------

The mechanics for double and triple integration have been wrapped up into the
functions :obj:`dblquad` and :obj:`tplquad`. These functions take the function
to  integrate and four, or six arguments, respectively. The limits of all
inner integrals need to be defined as functions.

An example of using double integration to compute several values of
:math:`I_{n}` is shown below:

    >>> from scipy.integrate import quad, dblquad
    >>> def I(n):
    ...     return dblquad(lambda t, x: np.exp(-x*t)/t**n, 0, np.inf, lambda x: 1, lambda x: np.inf)
    ...

    >>> print(I(4))
    (0.2500000000043577, 1.29830334693681e-08)
    >>> print(I(3))
    (0.33333333325010883, 1.3888461883425516e-08)
    >>> print(I(2))
    (0.4999999999985751, 1.3894083651858995e-08)


As example for non-constant limits consider the integral

.. math::

    I=\int_{y=0}^{1/2}\int_{x=0}^{1-2y} x y \, dx\, dy=\frac{1}{96}.


This integral can be evaluated using the expression below (Note the use of the
non-constant lambda functions for the upper limit of the inner integral):

>>> from scipy.integrate import dblquad
>>> area = dblquad(lambda x, y: x*y, 0, 0.5, lambda x: 0, lambda x: 1-2*x)
>>> area
(0.010416666666666668, 1.1564823173178715e-16)


For n-fold integration, scipy provides the function :obj:`nquad`. The
integration bounds are an iterable object: either a list of constant bounds,
or a list of functions for the non-constant integration bounds. The order of
integration (and therefore the bounds) is from the innermost integral to the
outermost one.

The integral from above

.. math::

    I_{n}=\int_{0}^{\infty}\int_{1}^{\infty}\frac{e^{-xt}}{t^{n}}\, dt\, dx=\frac{1}{n}

can be calculated as

>>> from scipy import integrate
>>> N = 5
>>> def f(t, x):
...    return np.exp(-x*t) / t**N
...
>>> integrate.nquad(f, [[1, np.inf],[0, np.inf]])
(0.20000000000002294, 1.2239614263187945e-08)

Note that the order of arguments for `f` must match the order of the
integration bounds; i.e., the inner integral with respect to :math:`t` is on
the interval :math:`[1, \infty]` and the outer integral with respect to
:math:`x` is on the interval :math:`[0, \infty]`.

Non-constant integration bounds can be treated in a similar manner; the
example from above

.. math::

    I=\int_{y=0}^{1/2}\int_{x=0}^{1-2y} x y \, dx\, dy=\frac{1}{96}.

can be evaluated by means of

>>> from scipy import integrate
>>> def f(x, y):
...     return x*y
...
>>> def bounds_y():
...     return [0, 0.5]
...
>>> def bounds_x(y):
...     return [0, 1-2*y]
...
>>> integrate.nquad(f, [bounds_x, bounds_y])
(0.010416666666666668, 4.101620128472366e-16)

which is the same result as before.

Gaussian quadrature
-------------------

A few functions are also provided in order to perform simple Gaussian
quadrature over a fixed interval. The first is :obj:`fixed_quad`, which
performs fixed-order Gaussian quadrature. The second function is
:obj:`quadrature`, which performs Gaussian quadrature of multiple
orders until the difference in the integral estimate is beneath some
tolerance supplied by the user. These functions both use the module
``scipy.special.orthogonal``, which can calculate the roots and quadrature
weights of a large variety of orthogonal polynomials (the polynomials
themselves are available as special functions returning instances of
the polynomial class --- e.g., :obj:`special.legendre <scipy.special.legendre>`).


Romberg Integration
-------------------

Romberg's method [WPR]_ is another method for numerically evaluating an
integral. See the help function for :func:`romberg` for further details.


Integrating using Samples
-------------------------

If the samples are equally-spaced and the number of samples available
is :math:`2^{k}+1` for some integer :math:`k`, then Romberg :obj:`romb`
integration can be used to obtain high-precision estimates of the
integral using the available samples. Romberg integration uses the
trapezoid rule at step-sizes related by a power of two and then
performs Richardson extrapolation on these estimates to approximate
the integral with a higher degree of accuracy.

In case of arbitrary spaced samples, the two functions :obj:`trapezoid`
and :obj:`simpson` are available. They are using Newton-Coates formulas
of order 1 and 2 respectively to perform integration. The trapezoidal rule
approximates the function as a straight line between adjacent points, while
Simpson's rule approximates the function between three adjacent points as a
parabola.

For an odd number of samples that are equally spaced Simpson's rule is exact
if the function is a polynomial of order 3 or less. If the samples are not
equally spaced, then the result is exact only if the function is a polynomial
of order 2 or less.

>>> import numpy as np
>>> def f1(x):
...    return x**2
...
>>> def f2(x):
...    return x**3
...
>>> x = np.array([1,3,4])
>>> y1 = f1(x)
>>> from scipy import integrate
>>> I1 = integrate.simpson(y1, x)
>>> print(I1)
21.0


This corresponds exactly to

.. math::

    \int_{1}^{4} x^2 \, dx = 21,

whereas integrating the second function

>>> y2 = f2(x)
>>> I2 = integrate.simpson(y2, x)
>>> print(I2)
61.5

does not correspond to

.. math::

    \int_{1}^{4} x^3 \, dx = 63.75

because the order of the polynomial in f2 is larger than two.

.. _quad-callbacks:

Faster integration using low-level callback functions
-----------------------------------------------------

A user desiring reduced integration times may pass a C function
pointer through `scipy.LowLevelCallable` to `quad`, `dblquad`,
`tplquad` or `nquad` and it will be integrated and return a result in
Python.  The performance increase here arises from two factors.  The
primary improvement is faster function evaluation, which is provided
by compilation of the function itself.  Additionally we have a speedup
provided by the removal of function calls between C and Python in
:obj:`quad`.  This method may provide a speed improvements of ~2x for
trivial functions such as sine but can produce a much more noticeable
improvements (10x+) for more complex functions.  This feature then, is
geared towards a user with numerically intensive integrations willing
to write a little C to reduce computation time significantly.

The approach can be used, for example, via `ctypes` in a few simple steps:

1.) Write an integrand function in C with the function signature
``double f(int n, double *x, void *user_data)``, where ``x`` is an
array containing the point the function f is evaluated at, and ``user_data``
to arbitrary additional data you want to provide.

.. code-block:: c

   /* testlib.c */
   double f(int n, double *x, void *user_data) {
       double c = *(double *)user_data;
       return c + x[0] - x[1] * x[2]; /* corresponds to c + x - y * z */
   }

2.) Now compile this file to a shared/dynamic library (a quick search will help
with this as it is OS-dependent). The user must link any math libraries,
etc., used.  On linux this looks like::

    $ gcc -shared -fPIC -o testlib.so testlib.c

The output library will be referred to as ``testlib.so``, but it may have a
different file extension. A library has now been created that can be loaded
into Python with `ctypes`.

3.) Load shared library into Python using `ctypes` and set ``restypes`` and
``argtypes`` - this allows SciPy to interpret the function correctly:

.. code:: python

   import os, ctypes
   from scipy import integrate, LowLevelCallable

   lib = ctypes.CDLL(os.path.abspath('testlib.so'))
   lib.f.restype = ctypes.c_double
   lib.f.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)

   c = ctypes.c_double(1.0)
   user_data = ctypes.cast(ctypes.pointer(c), ctypes.c_void_p)

   func = LowLevelCallable(lib.f, user_data)

The last ``void *user_data`` in the function is optional and can be omitted
(both in the C function and ctypes argtypes) if not needed. Note that the
coordinates are passed in as an array of doubles rather than a separate argument.

4.) Now integrate the library function as normally, here using `nquad`:

>>> integrate.nquad(func, [[0, 10], [-10, 0], [-1, 1]])
(1200.0, 1.1102230246251565e-11)

The Python tuple is returned as expected in a reduced amount of time.  All
optional parameters can be used with this method including specifying
singularities, infinite bounds, etc.

Ordinary differential equations (:func:`solve_ivp`)
---------------------------------------------------

Integrating a set of ordinary differential equations (ODEs) given
initial conditions is another useful example. The function
:obj:`solve_ivp` is available in SciPy for integrating a first-order
vector differential equation:

.. math::

    \frac{d\mathbf{y}}{dt}=\mathbf{f}\left(\mathbf{y},t\right),

given initial conditions :math:`\mathbf{y}\left(0\right)=y_{0}`, where
:math:`\mathbf{y}` is a length :math:`N` vector and :math:`\mathbf{f}`
is a mapping from :math:`\mathcal{R}^{N}` to :math:`\mathcal{R}^{N}.`
A higher-order ordinary differential equation can always be reduced to
a differential equation of this type by introducing intermediate
derivatives into the :math:`\mathbf{y}` vector.

For example, suppose it is desired to find the solution to the
following second-order differential equation:

.. math::

    \frac{d^{2}w}{dz^{2}}-zw(z)=0

with initial conditions :math:`w\left(0\right)=\frac{1}{\sqrt[3]{3^{2}}\Gamma\left(\frac{2}{3}\right)}` and :math:`\left.\frac{dw}{dz}\right|_{z=0}=-\frac{1}{\sqrt[3]{3}\Gamma\left(\frac{1}{3}\right)}.` It is known that the solution to this differential equation with these
boundary conditions is the Airy function

.. math::

    w=\textrm{Ai}\left(z\right),

which gives a means to check the integrator using :func:`special.airy <scipy.special.airy>`.

First, convert this ODE into standard form by setting
:math:`\mathbf{y}=\left[\frac{dw}{dz},w\right]` and :math:`t=z`. Thus,
the differential equation becomes

.. math::

    \frac{d\mathbf{y}}{dt}=\left[\begin{array}{c} ty_{1}\\ y_{0}\end{array}\right]=\left[\begin{array}{cc} 0 & t\\ 1 & 0\end{array}\right]\left[\begin{array}{c} y_{0}\\ y_{1}\end{array}\right]=\left[\begin{array}{cc} 0 & t\\ 1 & 0\end{array}\right]\mathbf{y}.

In other words,

.. math::

    \mathbf{f}\left(\mathbf{y},t\right)=\mathbf{A}\left(t\right)\mathbf{y}.

As an interesting reminder, if :math:`\mathbf{A}\left(t\right)`
commutes with :math:`\int_{0}^{t}\mathbf{A}\left(\tau\right)\, d\tau`
under matrix multiplication, then this linear differential equation
has an exact solution using the matrix exponential:

.. math::

    \mathbf{y}\left(t\right)=\exp\left(\int_{0}^{t}\mathbf{A}\left(\tau\right)d\tau\right)\mathbf{y}\left(0\right),

However, in this case, :math:`\mathbf{A}\left(t\right)` and its integral do not commute.

This differential equation can be solved using the function :obj:`solve_ivp`.
It requires the derivative, *fprime*, the time span `[t_start, t_end]`
and the initial conditions vector, *y0*, as input arguments and returns
an object whose *y* field is an array with consecutive solution values as
columns. The initial conditions are therefore given in the first output column.

>>> from scipy.integrate import solve_ivp
>>> from scipy.special import gamma, airy
>>> y1_0 = +1 / 3**(2/3) / gamma(2/3)
>>> y0_0 = -1 / 3**(1/3) / gamma(1/3)
>>> y0 = [y0_0, y1_0]
>>> def func(t, y):
...     return [t*y[1],y[0]]
...
>>> t_span = [0, 4]
>>> sol1 = solve_ivp(func, t_span, y0)
>>> print("sol1.t: {}".format(sol1.t))
sol1.t:    [0.         0.10097672 1.04643602 1.91060117 2.49872472 3.08684827
 3.62692846 4.        ]

As it can be seen `solve_ivp` determines its time steps automatically if not
specified otherwise. To compare the solution of `solve_ivp` with the `airy`
function the time vector created by `solve_ivp` is passed to the `airy` function.

>>> print("sol1.y[1]: {}".format(sol1.y[1]))
sol1.y[1]: [0.35502805 0.328952   0.12801343 0.04008508 0.01601291 0.00623879
 0.00356316 0.00405982]
>>> print("airy(sol.t)[0]:  {}".format(airy(sol1.t)[0]))
airy(sol.t)[0]: [0.35502805 0.328952   0.12804768 0.03995804 0.01575943 0.00562799
 0.00201689 0.00095156]

The solution of `solve_ivp` with its standard parameters shows a big deviation
to the airy function. To minimize this deviation, relative and absolute
tolerances can be used.

>>> rtol, atol = (1e-8, 1e-8)
>>> sol2 = solve_ivp(func, t_span, y0, rtol=rtol, atol=atol)
>>> print("sol2.y[1][::6]: {}".format(sol2.y[1][0::6]))
sol2.y[1][::6]: [0.35502805 0.19145234 0.06368989 0.0205917  0.00554734 0.00106409]
>>> print("airy(sol2.t)[0][::6]: {}".format(airy(sol2.t)[0][::6]))
airy(sol2.t)[0][::6]: [0.35502805 0.19145234 0.06368989 0.0205917  0.00554733 0.00106406]

To specify user defined time points for the solution of `solve_ivp`, `solve_ivp`
offers two possibilities that can also be used complementarily. By passing the `t_eval`
option to the function call `solve_ivp` returns the solutions of these time points
of `t_eval` in its output.

>>> import numpy as np
>>> t = np.linspace(0, 4, 100)
>>> sol3 = solve_ivp(func, t_span, y0, t_eval=t)

If the jacobian matrix of function is known, it can be passed to the `solve_ivp`
to achieve better results. Please be aware however that the default integration method
`RK45` does not support jacobian matrices and thereby another integration method has
to be chosen. One of the integration methods that support a jacobian matrix is the for
example the `Radau` method of following example.

>>> def gradient(t, y):
...     return [[0,t], [1,0]]
>>> sol4 = solve_ivp(func, t_span, y0, method='Radau', jac=gradient)

Solving a system with a banded Jacobian matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`odeint` can be told that the Jacobian is *banded*.  For a large
system of differential equations that are known to be stiff, this
can improve performance significantly.

As an example, we'll solve the 1-D Gray-Scott partial
differential equations using the method of lines [MOL]_.  The Gray-Scott equations
for the functions :math:`u(x, t)` and :math:`v(x, t)` on the interval
:math:`x \in [0, L]` are

.. math::

    \begin{split}
    \frac{\partial u}{\partial t} = D_u \frac{\partial^2 u}{\partial x^2} - uv^2 + f(1-u) \\
    \frac{\partial v}{\partial t} = D_v \frac{\partial^2 v}{\partial x^2} + uv^2 - (f + k)v \\
    \end{split}

where :math:`D_u` and :math:`D_v` are the diffusion coefficients of the
components :math:`u` and :math:`v`, respectively, and :math:`f` and :math:`k`
are constants.  (For more information about the system, see
http://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/)

We'll assume Neumann (i.e., "no flux") boundary conditions:

.. math::

    \frac{\partial u}{\partial x}(0,t) = 0, \quad
    \frac{\partial v}{\partial x}(0,t) = 0, \quad
    \frac{\partial u}{\partial x}(L,t) = 0, \quad
    \frac{\partial v}{\partial x}(L,t) = 0

To apply the method of lines, we discretize the :math:`x` variable by defining
the uniformly spaced grid of :math:`N` points :math:`\left\{x_0, x_1, \ldots, x_{N-1}\right\}`, with
:math:`x_0 = 0` and :math:`x_{N-1} = L`.
We define :math:`u_j(t) \equiv u(x_k, t)` and :math:`v_j(t) \equiv v(x_k, t)`, and
replace the :math:`x` derivatives with finite differences.  That is,

.. math::

    \frac{\partial^2 u}{\partial x^2}(x_j, t) \rightarrow
        \frac{u_{j-1}(t) - 2 u_{j}(t) + u_{j+1}(t)}{(\Delta x)^2}

We then have a system of :math:`2N` ordinary differential equations:

.. math::
   :label: interior

    \begin{split}
    \frac{du_j}{dt} = \frac{D_u}{(\Delta x)^2} \left(u_{j-1} - 2 u_{j} + u_{j+1}\right)
          -u_jv_j^2 + f(1 - u_j) \\
    \frac{dv_j}{dt} = \frac{D_v}{(\Delta x)^2} \left(v_{j-1} - 2 v_{j} + v_{j+1}\right)
          + u_jv_j^2 - (f + k)v_j
    \end{split}

For convenience, the :math:`(t)` arguments have been dropped.

To enforce the boundary conditions, we introduce "ghost" points
:math:`x_{-1}` and :math:`x_N`, and define :math:`u_{-1}(t) \equiv u_1(t)`,
:math:`u_N(t) \equiv u_{N-2}(t)`; :math:`v_{-1}(t)` and :math:`v_N(t)`
are defined analogously.

Then

.. math::
   :label: boundary0

    \begin{split}
    \frac{du_0}{dt} = \frac{D_u}{(\Delta x)^2} \left(2u_{1} - 2 u_{0}\right)
          -u_0v_0^2 + f(1 - u_0) \\
    \frac{dv_0}{dt} = \frac{D_v}{(\Delta x)^2} \left(2v_{1} - 2 v_{0}\right)
          + u_0v_0^2 - (f + k)v_0
    \end{split}

and

.. math::
   :label: boundaryL

    \begin{split}
    \frac{du_{N-1}}{dt} = \frac{D_u}{(\Delta x)^2} \left(2u_{N-2} - 2 u_{N-1}\right)
          -u_{N-1}v_{N-1}^2 + f(1 - u_{N-1}) \\
    \frac{dv_{N-1}}{dt} = \frac{D_v}{(\Delta x)^2} \left(2v_{N-2} - 2 v_{N-1}\right)
          + u_{N-1}v_{N-1}^2 - (f + k)v_{N-1}
    \end{split}

Our complete system of :math:`2N` ordinary differential equations is :eq:`interior`
for :math:`k = 1, 2, \ldots, N-2`, along with :eq:`boundary0` and :eq:`boundaryL`.

We can now starting implementing this system in code.  We must combine
:math:`\{u_k\}` and :math:`\{v_k\}` into a single vector of length :math:`2N`.
The two obvious choices are
:math:`\{u_0, u_1, \ldots, u_{N-1}, v_0, v_1, \ldots, v_{N-1}\}`
and
:math:`\{u_0, v_0, u_1, v_1, \ldots, u_{N-1}, v_{N-1}\}`.
Mathematically, it does not matter, but the choice affects how
efficiently `odeint` can solve the system.  The reason is in how
the order affects the pattern of the nonzero elements of the Jacobian matrix.


When the variables are ordered
as :math:`\{u_0, u_1, \ldots, u_{N-1}, v_0, v_1, \ldots, v_{N-1}\}`,
the pattern of nonzero elements of the Jacobian matrix is

.. math::

    \begin{smallmatrix}
       * & * & 0 & 0 & 0 & 0 & 0  &  * & 0 & 0 & 0 & 0 & 0 & 0 \\
       * & * & * & 0 & 0 & 0 & 0  &  0 & * & 0 & 0 & 0 & 0 & 0 \\
       0 & * & * & * & 0 & 0 & 0  &  0 & 0 & * & 0 & 0 & 0 & 0 \\
       0 & 0 & * & * & * & 0 & 0  &  0 & 0 & 0 & * & 0 & 0 & 0 \\
       0 & 0 & 0 & * & * & * & 0  &  0 & 0 & 0 & 0 & * & 0 & 0 \\
       0 & 0 & 0 & 0 & * & * & *  &  0 & 0 & 0 & 0 & 0 & * & 0 \\
       0 & 0 & 0 & 0 & 0 & * & *  &  0 & 0 & 0 & 0 & 0 & 0 & * \\
       * & 0 & 0 & 0 & 0 & 0 & 0  &  * & * & 0 & 0 & 0 & 0 & 0 \\
       0 & * & 0 & 0 & 0 & 0 & 0  &  * & * & * & 0 & 0 & 0 & 0 \\
       0 & 0 & * & 0 & 0 & 0 & 0  &  0 & * & * & * & 0 & 0 & 0 \\
       0 & 0 & 0 & * & 0 & 0 & 0  &  0 & 0 & * & * & * & 0 & 0 \\
       0 & 0 & 0 & 0 & * & 0 & 0  &  0 & 0 & 0 & * & * & * & 0 \\
       0 & 0 & 0 & 0 & 0 & * & 0  &  0 & 0 & 0 & 0 & * & * & * \\
       0 & 0 & 0 & 0 & 0 & 0 & *  &  0 & 0 & 0 & 0 & ) & * & * \\
    \end{smallmatrix}

The Jacobian pattern with variables interleaved
as :math:`\{u_0, v_0, u_1, v_1, \ldots, u_{N-1}, v_{N-1}\}` is

.. math::
    \begin{smallmatrix}
       * & * & * & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
       * & * & 0 & * & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
       * & 0 & * & * & * & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
       0 & * & * & * & 0 & * & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
       0 & 0 & * & 0 & * & * & * & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
       0 & 0 & 0 & * & * & * & 0 & * & 0 & 0 & 0 & 0 & 0 & 0 \\
       0 & 0 & 0 & 0 & * & 0 & * & * & * & 0 & 0 & 0 & 0 & 0 \\
       0 & 0 & 0 & 0 & 0 & * & * & * & 0 & * & 0 & 0 & 0 & 0 \\
       0 & 0 & 0 & 0 & 0 & 0 & * & 0 & * & * & * & 0 & 0 & 0 \\
       0 & 0 & 0 & 0 & 0 & 0 & 0 & * & * & * & 0 & * & 0 & 0 \\
       0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & * & 0 & * & * & * & 0 \\
       0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & * & * & * & 0 & * \\
       0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & * & 0 & * & * \\
       0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & * & * & * \\
    \end{smallmatrix}

In both cases, there are just five nontrivial diagonals, but
when the variables are interleaved, the bandwidth is much
smaller.
That is, the main diagonal and the two diagonals immediately
above and the two immediately below the main diagonal
are the nonzero diagonals.
This is important, because the inputs ``mu`` and ``ml``
of `odeint` are the upper and lower bandwidths of the
Jacobian matrix.  When the variables are interleaved,
``mu`` and ``ml`` are 2.  When the variables are stacked
with :math:`\{v_k\}` following :math:`\{u_k\}`, the upper
and lower bandwidths are :math:`N`.

With that decision made, we can write the function that
implements the system of differential equations.

First, we define the functions for the source and reaction
terms of the system::

    def G(u, v, f, k):
        return f * (1 - u) - u*v**2

    def H(u, v, f, k):
        return -(f + k) * v + u*v**2

Next, we define the function that computes the right-hand side
of the system of differential equations::

    def grayscott1d(y, t, f, k, Du, Dv, dx):
        """
        Differential equations for the 1-D Gray-Scott equations.

        The ODEs are derived using the method of lines.
        """
        # The vectors u and v are interleaved in y.  We define
        # views of u and v by slicing y.
        u = y[::2]
        v = y[1::2]

        # dydt is the return value of this function.
        dydt = np.empty_like(y)

        # Just like u and v are views of the interleaved vectors
        # in y, dudt and dvdt are views of the interleaved output
        # vectors in dydt.
        dudt = dydt[::2]
        dvdt = dydt[1::2]

        # Compute du/dt and dv/dt.  The end points and the interior points
        # are handled separately.
        dudt[0]    = G(u[0],    v[0],    f, k) + Du * (-2.0*u[0] + 2.0*u[1]) / dx**2
        dudt[1:-1] = G(u[1:-1], v[1:-1], f, k) + Du * np.diff(u,2) / dx**2
        dudt[-1]   = G(u[-1],   v[-1],   f, k) + Du * (- 2.0*u[-1] + 2.0*u[-2]) / dx**2
        dvdt[0]    = H(u[0],    v[0],    f, k) + Dv * (-2.0*v[0] + 2.0*v[1]) / dx**2
        dvdt[1:-1] = H(u[1:-1], v[1:-1], f, k) + Dv * np.diff(v,2) / dx**2
        dvdt[-1]   = H(u[-1],   v[-1],   f, k) + Dv * (-2.0*v[-1] + 2.0*v[-2]) / dx**2

        return dydt

We won't implement a function to compute the Jacobian, but we will tell
`odeint` that the Jacobian matrix is banded.  This allows the underlying
solver (LSODA) to avoid computing values that it knows are zero.  For a large
system, this improves the performance significantly, as demonstrated in the
following ipython session.

First, we define the required inputs::

    In [31]: y0 = np.random.randn(5000)

    In [32]: t = np.linspace(0, 50, 11)

    In [33]: f = 0.024

    In [34]: k = 0.055

    In [35]: Du = 0.01

    In [36]: Dv = 0.005

    In [37]: dx = 0.025

Time the computation without taking advantage of the banded structure
of the Jacobian matrix::

    In [38]: %timeit sola = odeint(grayscott1d, y0, t, args=(f, k, Du, Dv, dx))
    1 loop, best of 3: 25.2 s per loop

Now set ``ml=2`` and ``mu=2``, so `odeint` knows that the Jacobian matrix
is banded::

    In [39]: %timeit solb = odeint(grayscott1d, y0, t, args=(f, k, Du, Dv, dx), ml=2, mu=2)
    10 loops, best of 3: 191 ms per loop

That is quite a bit faster!

Let's ensure that they have computed the same result::

    In [41]: np.allclose(sola, solb)
    Out[41]: True

References
~~~~~~~~~~

.. [WPR] https://en.wikipedia.org/wiki/Romberg's_method

.. [MOL] https://en.wikipedia.org/wiki/Method_of_lines
