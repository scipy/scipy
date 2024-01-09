Optimization (:mod:`scipy.optimize`)
====================================

.. sectionauthor:: Travis E. Oliphant

.. sectionauthor:: Pauli Virtanen

.. sectionauthor:: Denis Laxalde

.. currentmodule:: scipy.optimize

.. contents::

The :mod:`scipy.optimize` package provides several commonly used
optimization algorithms. A detailed listing is available:
:mod:`scipy.optimize` (can also be found by ``help(scipy.optimize)``).


Unconstrained minimization of multivariate scalar functions (:func:`minimize`)
------------------------------------------------------------------------------

The :func:`minimize` function provides a common interface to unconstrained
and constrained minimization algorithms for multivariate scalar functions
in `scipy.optimize`. To demonstrate the minimization function, consider the
problem of minimizing the Rosenbrock function of :math:`N` variables:

.. math::

    f\left(\mathbf{x}\right)=\sum_{i=1}^{N-1}100\left(x_{i+1}-x_{i}^{2}\right)^{2}+\left(1-x_{i}\right)^{2}.

The minimum value of this function is 0 which is achieved when
:math:`x_{i}=1.`

Note that the Rosenbrock function and its derivatives are included in
`scipy.optimize`. The implementations shown in the following sections
provide examples of how to define an objective function as well as its
jacobian and hessian functions. Objective functions in `scipy.optimize`
expect a numpy array as their first parameter which is to be optimized
and must return a float value. The exact calling signature must be
``f(x, *args)`` where ``x`` represents a numpy array and ``args``
a tuple of additional arguments supplied to the objective function.

Nelder-Mead Simplex algorithm (``method='Nelder-Mead'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the example below, the :func:`minimize` routine is used
with the *Nelder-Mead* simplex algorithm (selected through the ``method``
parameter):

    >>> import numpy as np
    >>> from scipy.optimize import minimize

    >>> def rosen(x):
    ...     """The Rosenbrock function"""
    ...     return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

    >>> x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    >>> res = minimize(rosen, x0, method='nelder-mead',
    ...                options={'xatol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 339
             Function evaluations: 571

    >>> print(res.x)
    [1. 1. 1. 1. 1.]

The simplex algorithm is probably the simplest way to minimize a fairly
well-behaved function. It requires only function evaluations and is a good
choice for simple minimization problems. However, because it does not use
any gradient evaluations, it may take longer to find the minimum.

Another optimization algorithm that needs only function calls to find
the minimum is *Powell*'s method available by setting ``method='powell'`` in
:func:`minimize`.

To demonstrate how to supply additional arguments to an objective function,
let us minimize the Rosenbrock function with an additional scaling factor `a`
and an offset `b`:

.. math::

    f\left(\mathbf{x}, a, b\right)=\sum_{i=1}^{N-1}a\left(x_{i+1}-x_{i}^{2}\right)^{2}+\left(1-x_{i}\right)^{2} + b.

Again using the :func:`minimize` routine this can be solved by the following
code block for the example parameters ``a=0.5`` and ``b=1``.

    >>> def rosen_with_args(x, a, b):
    ...     """The Rosenbrock function with additional arguments"""
    ...     return sum(a*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0) + b

    >>> x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    >>> res = minimize(rosen_with_args, x0, method='nelder-mead',
    ...	               args=(0.5, 1.), options={'xatol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 1.000000
             Iterations: 319 # may vary
             Function evaluations: 525 # may vary

    >>> print(res.x)
    [1.         1.         1.         1.         0.99999999]

As an alternative to using the ``args`` parameter of :func:`minimize`, simply
wrap the objective function in a new function that accepts only ``x``. This
approach is also useful when it is necessary to pass additional parameters to
the objective function as keyword arguments.

    >>> def rosen_with_args(x, a, *, b):  # b is a keyword-only argument
    ...     return sum(a*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0) + b
    >>> def wrapped_rosen_without_args(x):
    ...     return rosen_with_args(x, 0.5, b=1.)  # pass in `a` and `b`
    >>> x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    >>> res = minimize(wrapped_rosen_without_args, x0, method='nelder-mead',
    ...                options={'xatol': 1e-8,})
    >>> print(res.x)
    [1.         1.         1.         1.         0.99999999]

Another alternative is to use :py:func:`functools.partial`.

    >>> from functools import partial
    >>> partial_rosen = partial(rosen_with_args, a=0.5, b=1.)
    >>> res = minimize(partial_rosen, x0, method='nelder-mead',
    ...                options={'xatol': 1e-8,})
    >>> print(res.x)
    [1.         1.         1.         1.         0.99999999]

Broyden-Fletcher-Goldfarb-Shanno algorithm (``method='BFGS'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to converge more quickly to the solution, this routine uses
the gradient of the objective function. If the gradient is not given
by the user, then it is estimated using first-differences. The
Broyden-Fletcher-Goldfarb-Shanno (BFGS) method typically requires
fewer function calls than the simplex algorithm even when the gradient
must be estimated.

To demonstrate this algorithm, the Rosenbrock function is again used.
The gradient of the Rosenbrock function is the vector:

.. math::
   :nowrap:

    \begin{eqnarray*} \frac{\partial f}{\partial x_{j}} & = & \sum_{i=1}^{N}200\left(x_{i}-x_{i-1}^{2}\right)\left(\delta_{i,j}-2x_{i-1}\delta_{i-1,j}\right)-2\left(1-x_{i-1}\right)\delta_{i-1,j}.\\  & = & 200\left(x_{j}-x_{j-1}^{2}\right)-400x_{j}\left(x_{j+1}-x_{j}^{2}\right)-2\left(1-x_{j}\right).\end{eqnarray*}

This expression is valid for the interior derivatives. Special cases
are

.. math::
   :nowrap:

    \begin{eqnarray*} \frac{\partial f}{\partial x_{0}} & = & -400x_{0}\left(x_{1}-x_{0}^{2}\right)-2\left(1-x_{0}\right),\\ \frac{\partial f}{\partial x_{N-1}} & = & 200\left(x_{N-1}-x_{N-2}^{2}\right).\end{eqnarray*}

A Python function which computes this gradient is constructed by the
code-segment:

    >>> def rosen_der(x):
    ...     xm = x[1:-1]
    ...     xm_m1 = x[:-2]
    ...     xm_p1 = x[2:]
    ...     der = np.zeros_like(x)
    ...     der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
    ...     der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    ...     der[-1] = 200*(x[-1]-x[-2]**2)
    ...     return der

This gradient information is specified in the :func:`minimize` function
through the ``jac`` parameter as illustrated below.


    >>> res = minimize(rosen, x0, method='BFGS', jac=rosen_der,
    ...                options={'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 25                     # may vary
             Function evaluations: 30
             Gradient evaluations: 30
    >>> res.x
    array([1., 1., 1., 1., 1.])

Avoiding Redundant Calculation
""""""""""""""""""""""""""""""

It is common for the objective function and its gradient to share parts of the
calculation. For instance, consider the following problem.

    >>> def f(x):
    ...    return -expensive(x[0])**2
    >>>
    >>> def df(x):
    ...     return -2 * expensive(x[0]) * dexpensive(x[0])
    >>>
    >>> def expensive(x):
    ...     # this function is computationally expensive!
    ...     expensive.count += 1  # let's keep track of how many times it runs
    ...     return np.sin(x)
    >>> expensive.count = 0
    >>>
    >>> def dexpensive(x):
    ...     return np.cos(x)
    >>>
    >>> res = minimize(f, 0.5, jac=df)
    >>> res.fun
    -0.9999999999999174
    >>> res.nfev, res.njev
    6, 6
    >>> expensive.count
    12

Here, ``expensive`` is called 12 times: six times in the objective function and
six times from the gradient. One way of reducing redundant calculations is to
create a single function that returns both the objective function and the
gradient.

    >>> def f_and_df(x):
    ...     expensive_value = expensive(x[0])
    ...     return (-expensive_value**2,  # objective function
    ...             -2*expensive_value*dexpensive(x[0]))  # gradient
    >>>
    >>> expensive.count = 0  # reset the counter
    >>> res = minimize(f_and_df, 0.5, jac=True)
    >>> res.fun
    -0.9999999999999174
    >>> expensive.count
    6

When we call minimize, we specify ``jac==True`` to indicate that the provided
function returns both the objective function and its gradient. While
convenient, not all :mod:`scipy.optimize` functions support this feature,
and moreover, it is only for sharing calculations between the function and its
gradient, whereas in some problems we will want to share calculations with the
Hessian (second derivative of the objective function) and constraints. A more
general approach is to memoize the expensive parts of the calculation. In
simple situations, this can be accomplished with the
:func:`functools.lru_cache` wrapper.

    >>> from functools import lru_cache
    >>> expensive.count = 0  # reset the counter
    >>> expensive = lru_cache(expensive)
    >>> res = minimize(f, 0.5, jac=df)
    >>> res.fun
    -0.9999999999999174
    >>> expensive.count
    6


Newton-Conjugate-Gradient algorithm (``method='Newton-CG'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Newton-Conjugate Gradient algorithm is a modified Newton's
method and uses a conjugate gradient algorithm to (approximately) invert
the local Hessian [NW]_.  Newton's method is based on fitting the function
locally to a quadratic form:

.. math::

    f\left(\mathbf{x}\right)\approx f\left(\mathbf{x}_{0}\right)+\nabla f\left(\mathbf{x}_{0}\right)\cdot\left(\mathbf{x}-\mathbf{x}_{0}\right)+\frac{1}{2}\left(\mathbf{x}-\mathbf{x}_{0}\right)^{T}\mathbf{H}\left(\mathbf{x}_{0}\right)\left(\mathbf{x}-\mathbf{x}_{0}\right).

where :math:`\mathbf{H}\left(\mathbf{x}_{0}\right)` is a matrix of second-derivatives (the Hessian). If the Hessian is
positive definite then the local minimum of this function can be found
by setting the gradient of the quadratic form to zero, resulting in

.. math::

    \mathbf{x}_{\textrm{opt}}=\mathbf{x}_{0}-\mathbf{H}^{-1}\nabla f.

The inverse of the Hessian is evaluated using the conjugate-gradient
method. An example of employing this method to minimizing the
Rosenbrock function is given below. To take full advantage of the
Newton-CG method, a function which computes the Hessian must be
provided. The Hessian matrix itself does not need to be constructed,
only a vector which is the product of the Hessian with an arbitrary
vector needs to be available to the minimization routine. As a result,
the user can provide either a function to compute the Hessian matrix,
or a function to compute the product of the Hessian with an arbitrary
vector.


Full Hessian example:
"""""""""""""""""""""

The Hessian of the Rosenbrock function is

.. math::
   :nowrap:

    \begin{eqnarray*} H_{ij}=\frac{\partial^{2}f}{\partial x_{i}\partial x_{j}} & = & 200\left(\delta_{i,j}-2x_{i-1}\delta_{i-1,j}\right)-400x_{i}\left(\delta_{i+1,j}-2x_{i}\delta_{i,j}\right)-400\delta_{i,j}\left(x_{i+1}-x_{i}^{2}\right)+2\delta_{i,j},\\  & = & \left(202+1200x_{i}^{2}-400x_{i+1}\right)\delta_{i,j}-400x_{i}\delta_{i+1,j}-400x_{i-1}\delta_{i-1,j},\end{eqnarray*}

if :math:`i,j\in\left[1,N-2\right]` with :math:`i,j\in\left[0,N-1\right]` defining the :math:`N\times N` matrix. Other non-zero entries of the matrix are

.. math::
   :nowrap:

    \begin{eqnarray*} \frac{\partial^{2}f}{\partial x_{0}^{2}} & = & 1200x_{0}^{2}-400x_{1}+2,\\ \frac{\partial^{2}f}{\partial x_{0}\partial x_{1}}=\frac{\partial^{2}f}{\partial x_{1}\partial x_{0}} & = & -400x_{0},\\ \frac{\partial^{2}f}{\partial x_{N-1}\partial x_{N-2}}=\frac{\partial^{2}f}{\partial x_{N-2}\partial x_{N-1}} & = & -400x_{N-2},\\ \frac{\partial^{2}f}{\partial x_{N-1}^{2}} & = & 200.\end{eqnarray*}

For example, the Hessian when :math:`N=5` is

.. math::

    \mathbf{H}=\begin{bmatrix} 1200x_{0}^{2}+2\mkern-2em\\&1200x_{1}^{2}+202\mkern-2em\\&&1200x_{1}^{2}+202\mkern-2em\\&&&1200x_{3}^{2}+202\mkern-1em\\&&&&200\end{bmatrix}-400\begin{bmatrix} x_1 & x_0 \\ x_0 & x_2 & x_1 \\ & x_1 & x_3 & x_2\\ & & x_2 & x_4 & x_3 \\ & & & x_3 & 0\end{bmatrix}.

The code which computes this Hessian along with the code to minimize
the function using Newton-CG method is shown in the following example:

    >>> def rosen_hess(x):
    ...     x = np.asarray(x)
    ...     H = np.diag(-400*x[:-1],1) - np.diag(400*x[:-1],-1)
    ...     diagonal = np.zeros_like(x)
    ...     diagonal[0] = 1200*x[0]**2-400*x[1]+2
    ...     diagonal[-1] = 200
    ...     diagonal[1:-1] = 202 + 1200*x[1:-1]**2 - 400*x[2:]
    ...     H = H + np.diag(diagonal)
    ...     return H

    >>> res = minimize(rosen, x0, method='Newton-CG',
    ...                jac=rosen_der, hess=rosen_hess,
    ...                options={'xtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 19                       # may vary
             Function evaluations: 22
             Gradient evaluations: 19
             Hessian evaluations: 19
    >>> res.x
    array([1.,  1.,  1.,  1.,  1.])


Hessian product example:
""""""""""""""""""""""""

For larger minimization problems, storing the entire Hessian matrix can
consume considerable time and memory. The Newton-CG algorithm only needs
the product of the Hessian times an arbitrary vector. As a result, the user
can supply code to compute this product rather than the full Hessian by
giving a ``hess`` function which take the minimization vector as the first
argument and the arbitrary vector as the second argument (along with extra
arguments passed to the function to be minimized). If possible, using
Newton-CG with the Hessian product option is probably the fastest way to
minimize the function.

In this case, the product of the Rosenbrock Hessian with an arbitrary
vector is not difficult to compute. If :math:`\mathbf{p}` is the arbitrary
vector, then :math:`\mathbf{H}\left(\mathbf{x}\right)\mathbf{p}` has
elements:

.. math::

    \mathbf{H}\left(\mathbf{x}\right)\mathbf{p}=\begin{bmatrix} \left(1200x_{0}^{2}-400x_{1}+2\right)p_{0}-400x_{0}p_{1}\\ \vdots\\ -400x_{i-1}p_{i-1}+\left(202+1200x_{i}^{2}-400x_{i+1}\right)p_{i}-400x_{i}p_{i+1}\\ \vdots\\ -400x_{N-2}p_{N-2}+200p_{N-1}\end{bmatrix}.

Code which makes use of this Hessian product to minimize the
Rosenbrock function using :func:`minimize` follows:

    >>> def rosen_hess_p(x, p):
    ...     x = np.asarray(x)
    ...     Hp = np.zeros_like(x)
    ...     Hp[0] = (1200*x[0]**2 - 400*x[1] + 2)*p[0] - 400*x[0]*p[1]
    ...     Hp[1:-1] = -400*x[:-2]*p[:-2]+(202+1200*x[1:-1]**2-400*x[2:])*p[1:-1] \
    ...                -400*x[1:-1]*p[2:]
    ...     Hp[-1] = -400*x[-2]*p[-2] + 200*p[-1]
    ...     return Hp

    >>> res = minimize(rosen, x0, method='Newton-CG',
    ...                jac=rosen_der, hessp=rosen_hess_p,
    ...                options={'xtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 20                    # may vary
             Function evaluations: 23
             Gradient evaluations: 20
             Hessian evaluations: 44
    >>> res.x
    array([1., 1., 1., 1., 1.])


According to [NW]_ p. 170 the ``Newton-CG`` algorithm can be inefficient
when the Hessian is ill-conditioned because of the poor quality search directions
provided by the method in those situations. The method ``trust-ncg``,
according to the authors, deals more effectively with this problematic situation
and will be described next.

Trust-Region Newton-Conjugate-Gradient Algorithm (``method='trust-ncg'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Newton-CG`` method is a line search method: it finds a direction
of search minimizing a quadratic approximation of the function and then uses
a line search algorithm to find the (nearly) optimal step size in that direction.
An alternative approach is to, first, fix the step size limit :math:`\Delta` and then find the
optimal step :math:`\mathbf{p}` inside the given trust-radius by solving
the following quadratic subproblem:

.. math::
   :nowrap:

   \begin{eqnarray*}
      \min_{\mathbf{p}} f\left(\mathbf{x}_{k}\right)+\nabla f\left(\mathbf{x}_{k}\right)\cdot\mathbf{p}+\frac{1}{2}\mathbf{p}^{T}\mathbf{H}\left(\mathbf{x}_{k}\right)\mathbf{p};&\\
      \text{subject to: } \|\mathbf{p}\|\le \Delta.&
    \end{eqnarray*}

The solution is then updated :math:`\mathbf{x}_{k+1} = \mathbf{x}_{k} + \mathbf{p}` and
the trust-radius :math:`\Delta` is adjusted according to the degree of agreement of the quadratic
model with the real function. This family of methods is known as trust-region methods.
The ``trust-ncg`` algorithm is a trust-region method that uses a conjugate gradient algorithm
to solve the trust-region subproblem [NW]_.


Full Hessian example:
"""""""""""""""""""""

    >>> res = minimize(rosen, x0, method='trust-ncg',
    ...                jac=rosen_der, hess=rosen_hess,
    ...                options={'gtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 20                    # may vary
             Function evaluations: 21
             Gradient evaluations: 20
             Hessian evaluations: 19
    >>> res.x
    array([1., 1., 1., 1., 1.])

Hessian product example:
""""""""""""""""""""""""

    >>> res = minimize(rosen, x0, method='trust-ncg',
    ...                jac=rosen_der, hessp=rosen_hess_p,
    ...                options={'gtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 20                    # may vary
             Function evaluations: 21
             Gradient evaluations: 20
             Hessian evaluations: 0
    >>> res.x
    array([1., 1., 1., 1., 1.])

Trust-Region Truncated Generalized Lanczos / Conjugate Gradient Algorithm (``method='trust-krylov'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to the ``trust-ncg`` method, the ``trust-krylov`` method is a method
suitable for large-scale problems as it uses the hessian only as linear
operator by means of matrix-vector products.
It solves the quadratic subproblem more accurately than the ``trust-ncg``
method.

.. math::
   :nowrap:

   \begin{eqnarray*}
      \min_{\mathbf{p}} f\left(\mathbf{x}_{k}\right)+\nabla f\left(\mathbf{x}_{k}\right)\cdot\mathbf{p}+\frac{1}{2}\mathbf{p}^{T}\mathbf{H}\left(\mathbf{x}_{k}\right)\mathbf{p};&\\
      \text{subject to: } \|\mathbf{p}\|\le \Delta.&
    \end{eqnarray*}

This method wraps the [TRLIB]_ implementation of the [GLTR]_ method solving
exactly a trust-region subproblem restricted to a truncated Krylov subspace.
For indefinite problems it is usually better to use this method as it reduces
the number of nonlinear iterations at the expense of few more matrix-vector
products per subproblem solve in comparison to the ``trust-ncg`` method.

Full Hessian example:
"""""""""""""""""""""

    >>> res = minimize(rosen, x0, method='trust-krylov',
    ...                jac=rosen_der, hess=rosen_hess,
    ...                options={'gtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 19                    # may vary
             Function evaluations: 20
             Gradient evaluations: 20
             Hessian evaluations: 18
    >>> res.x
    array([1., 1., 1., 1., 1.])

Hessian product example:
""""""""""""""""""""""""

    >>> res = minimize(rosen, x0, method='trust-krylov',
    ...                jac=rosen_der, hessp=rosen_hess_p,
    ...                options={'gtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 19                    # may vary
             Function evaluations: 20
             Gradient evaluations: 20
             Hessian evaluations: 0
    >>> res.x
    array([1., 1., 1., 1., 1.])

.. [TRLIB] F. Lenders, C. Kirches, A. Potschka: "trlib: A vector-free
           implementation of the GLTR method for iterative solution of
           the trust region problem", :arxiv:`1611.04718`

.. [GLTR]  N. Gould, S. Lucidi, M. Roma, P. Toint: "Solving the
           Trust-Region Subproblem using the Lanczos Method",
           SIAM J. Optim., 9(2), 504--525, (1999).
           :doi:`10.1137/S1052623497322735`


Trust-Region Nearly Exact Algorithm (``method='trust-exact'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All methods ``Newton-CG``, ``trust-ncg`` and ``trust-krylov`` are suitable for dealing with
large-scale problems (problems with thousands of variables). That is because the conjugate
gradient algorithm approximately solve the trust-region subproblem (or invert the Hessian)
by iterations without the explicit Hessian factorization. Since only the product of the Hessian
with an arbitrary vector is needed, the algorithm is specially suited for dealing
with sparse Hessians, allowing low storage requirements and significant time savings for
those sparse problems.

For medium-size problems, for which the storage and factorization cost of the Hessian are not critical,
it is possible to obtain a solution within fewer iteration by solving the trust-region subproblems
almost exactly. To achieve that, a certain nonlinear equations is solved iteratively for each quadratic
subproblem [CGT]_. This solution requires usually 3 or 4 Cholesky factorizations of the
Hessian matrix. As the result, the method converges in fewer number of iterations
and takes fewer evaluations of the objective function than the other implemented
trust-region methods. The Hessian product option is not supported by this algorithm. An
example using the Rosenbrock function follows:


    >>> res = minimize(rosen, x0, method='trust-exact',
    ...                jac=rosen_der, hess=rosen_hess,
    ...                options={'gtol': 1e-8, 'disp': True})
    Optimization terminated successfully.
             Current function value: 0.000000
             Iterations: 13                    # may vary
             Function evaluations: 14
             Gradient evaluations: 13
             Hessian evaluations: 14
    >>> res.x
    array([1., 1., 1., 1., 1.])


.. [NW] J. Nocedal, S.J. Wright "Numerical optimization."
	2nd edition. Springer Science (2006).
.. [CGT] Conn, A. R., Gould, N. I., & Toint, P. L.
        "Trust region methods". Siam. (2000). pp. 169-200.


.. _tutorial-sqlsp:

Constrained minimization of multivariate scalar functions (:func:`minimize`)
----------------------------------------------------------------------------

The :func:`minimize` function provides algorithms for constrained minimization,
namely ``'trust-constr'`` ,  ``'SLSQP'`` and ``'COBYLA'``. They require the constraints
to be defined using slightly different structures. The method ``'trust-constr'`` requires
the  constraints to be defined as a sequence of objects :func:`LinearConstraint` and
:func:`NonlinearConstraint`. Methods ``'SLSQP'`` and ``'COBYLA'``, on the other hand,
require constraints to be defined  as a sequence of dictionaries, with keys
``type``, ``fun`` and ``jac``.

As an example let us consider the constrained minimization of the Rosenbrock function:

.. math::
   :nowrap:

     \begin{eqnarray*} \min_{x_0, x_1} & ~~100\left(x_{1}-x_{0}^{2}\right)^{2}+\left(1-x_{0}\right)^{2} &\\
                     \text{subject to: } & x_0 + 2 x_1 \leq 1 & \\
		                         & x_0^2 + x_1 \leq 1  & \\
		                         & x_0^2 - x_1 \leq 1  & \\
					 & 2 x_0 + x_1 = 1 & \\
					 & 0 \leq  x_0  \leq 1 & \\
					 & -0.5 \leq  x_1  \leq 2.0. & \end{eqnarray*}

This optimization problem has the unique solution :math:`[x_0, x_1] = [0.4149,~ 0.1701]`,
for which only the first and fourth constraints are active.


Trust-Region Constrained Algorithm (``method='trust-constr'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The trust-region constrained method deals with constrained minimization problems of the form:

.. math::
   :nowrap:

     \begin{eqnarray*} \min_x & f(x) & \\
          \text{subject to: } & ~~~ c^l  \leq c(x) \leq c^u, &\\
           &  x^l  \leq x \leq x^u. & \end{eqnarray*}

When :math:`c^l_j = c^u_j` the method reads the :math:`j`-th constraint as an
equality constraint and deals with it accordingly. Besides that, one-sided constraint
can be specified by setting the upper or lower bound to ``np.inf`` with the appropriate sign.

The implementation is based on [EQSQP]_ for equality-constraint problems and on [TRIP]_
for problems with inequality constraints. Both are trust-region type algorithms suitable
for large-scale problems.


Defining Bounds Constraints:
""""""""""""""""""""""""""""

The bound constraints  :math:`0 \leq  x_0  \leq 1` and :math:`-0.5 \leq  x_1  \leq 2.0`
are defined using a :func:`Bounds` object.

    >>> from scipy.optimize import Bounds
    >>> bounds = Bounds([0, -0.5], [1.0, 2.0])

Defining Linear Constraints:
""""""""""""""""""""""""""""

The constraints :math:`x_0 + 2 x_1 \leq 1`
and :math:`2 x_0 + x_1 = 1` can be written in the linear constraint standard format:

.. math::
   :nowrap:

     \begin{equation*} \begin{bmatrix}-\infty \\1\end{bmatrix} \leq
      \begin{bmatrix} 1& 2 \\ 2& 1\end{bmatrix}
       \begin{bmatrix} x_0 \\x_1\end{bmatrix} \leq
        \begin{bmatrix} 1 \\ 1\end{bmatrix},\end{equation*}

and defined using a :func:`LinearConstraint` object.

    >>> from scipy.optimize import LinearConstraint
    >>> linear_constraint = LinearConstraint([[1, 2], [2, 1]], [-np.inf, 1], [1, 1])

Defining Nonlinear Constraints:
"""""""""""""""""""""""""""""""
The nonlinear constraint:

.. math::
   :nowrap:

     \begin{equation*} c(x) =
     \begin{bmatrix} x_0^2 + x_1 \\ x_0^2 - x_1\end{bmatrix}
      \leq
      \begin{bmatrix} 1 \\ 1\end{bmatrix}, \end{equation*}

with Jacobian matrix:

.. math::
   :nowrap:

     \begin{equation*} J(x) =
     \begin{bmatrix} 2x_0 & 1 \\ 2x_0 & -1\end{bmatrix},\end{equation*}

and linear combination of the Hessians:

.. math::
   :nowrap:

     \begin{equation*} H(x, v) = \sum_{i=0}^1 v_i \nabla^2 c_i(x) =
     v_0\begin{bmatrix} 2 & 0 \\ 0 & 0\end{bmatrix} +
     v_1\begin{bmatrix} 2 & 0 \\ 0 & 0\end{bmatrix},
     \end{equation*}

is defined using a :func:`NonlinearConstraint` object.

    >>> def cons_f(x):
    ...     return [x[0]**2 + x[1], x[0]**2 - x[1]]
    >>> def cons_J(x):
    ...     return [[2*x[0], 1], [2*x[0], -1]]
    >>> def cons_H(x, v):
    ...     return v[0]*np.array([[2, 0], [0, 0]]) + v[1]*np.array([[2, 0], [0, 0]])
    >>> from scipy.optimize import NonlinearConstraint
    >>> nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1, jac=cons_J, hess=cons_H)

Alternatively, it is also possible to define the Hessian :math:`H(x, v)`
as a sparse matrix,

    >>> from scipy.sparse import csc_matrix
    >>> def cons_H_sparse(x, v):
    ...     return v[0]*csc_matrix([[2, 0], [0, 0]]) + v[1]*csc_matrix([[2, 0], [0, 0]])
    >>> nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1,
    ...                                            jac=cons_J, hess=cons_H_sparse)

or as a :obj:`~scipy.sparse.linalg.LinearOperator` object.

    >>> from scipy.sparse.linalg import LinearOperator
    >>> def cons_H_linear_operator(x, v):
    ...     def matvec(p):
    ...         return np.array([p[0]*2*(v[0]+v[1]), 0])
    ...     return LinearOperator((2, 2), matvec=matvec)
    >>> nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1,
    ...                                           jac=cons_J, hess=cons_H_linear_operator)

When the evaluation of the Hessian :math:`H(x, v)`
is difficult to implement or computationally infeasible, one may use :class:`HessianUpdateStrategy`.
Currently available strategies are :class:`BFGS` and :class:`SR1`.

    >>> from scipy.optimize import BFGS
    >>> nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1, jac=cons_J, hess=BFGS())

Alternatively, the Hessian may be approximated using finite differences.

    >>> nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1, jac=cons_J, hess='2-point')

The Jacobian of the constraints can be approximated by finite differences as well. In this case,
however, the Hessian cannot be computed with finite differences and needs to
be provided by the user or defined using :class:`HessianUpdateStrategy`.

    >>> nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, 1, jac='2-point', hess=BFGS())


Solving the Optimization Problem:
"""""""""""""""""""""""""""""""""
The optimization problem is solved using:

    >>> x0 = np.array([0.5, 0])
    >>> res = minimize(rosen, x0, method='trust-constr', jac=rosen_der, hess=rosen_hess,
    ...                constraints=[linear_constraint, nonlinear_constraint],
    ...                options={'verbose': 1}, bounds=bounds)
    # may vary
    `gtol` termination condition is satisfied.
    Number of iterations: 12, function evaluations: 8, CG iterations: 7, optimality: 2.99e-09, constraint violation: 1.11e-16, execution time: 0.016 s.
    >>> print(res.x)
    [0.41494531 0.17010937]

When needed, the objective function Hessian can be defined using a :obj:`~scipy.sparse.linalg.LinearOperator` object,

    >>> def rosen_hess_linop(x):
    ...     def matvec(p):
    ...         return rosen_hess_p(x, p)
    ...     return LinearOperator((2, 2), matvec=matvec)
    >>> res = minimize(rosen, x0, method='trust-constr', jac=rosen_der, hess=rosen_hess_linop,
    ...                constraints=[linear_constraint, nonlinear_constraint],
    ...                options={'verbose': 1}, bounds=bounds)
    # may vary
    `gtol` termination condition is satisfied.
    Number of iterations: 12, function evaluations: 8, CG iterations: 7, optimality: 2.99e-09, constraint violation: 1.11e-16, execution time: 0.018 s.
    >>> print(res.x)
    [0.41494531 0.17010937]

or a Hessian-vector product through the parameter ``hessp``.

    >>> res = minimize(rosen, x0, method='trust-constr', jac=rosen_der, hessp=rosen_hess_p,
    ...                constraints=[linear_constraint, nonlinear_constraint],
    ...                options={'verbose': 1}, bounds=bounds)
    # may vary
    `gtol` termination condition is satisfied.
    Number of iterations: 12, function evaluations: 8, CG iterations: 7, optimality: 2.99e-09, constraint violation: 1.11e-16, execution time: 0.018 s.
    >>> print(res.x)
    [0.41494531 0.17010937]

Alternatively, the first and second derivatives of the objective function can be approximated.
For instance,  the Hessian can be approximated with :func:`SR1` quasi-Newton approximation
and the gradient with finite differences.

    >>> from scipy.optimize import SR1
    >>> res = minimize(rosen, x0, method='trust-constr',  jac="2-point", hess=SR1(),
    ...                constraints=[linear_constraint, nonlinear_constraint],
    ...                options={'verbose': 1}, bounds=bounds)
    # may vary
    `gtol` termination condition is satisfied.
    Number of iterations: 12, function evaluations: 24, CG iterations: 7, optimality: 4.48e-09, constraint violation: 0.00e+00, execution time: 0.016 s.
    >>> print(res.x)
    [0.41494531 0.17010937]


.. [TRIP] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal. 1999.
    An interior point algorithm for large-scale nonlinear  programming.
    SIAM Journal on Optimization 9.4: 877-900.
.. [EQSQP] Lalee, Marucha, Jorge Nocedal, and Todd Plantega. 1998. On the
    implementation of an algorithm for large-scale equality constrained
    optimization. SIAM Journal on Optimization 8.3: 682-706.

Sequential Least SQuares Programming (SLSQP) Algorithm (``method='SLSQP'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The SLSQP method deals with constrained minimization problems of the form:

.. math::
   :nowrap:

     \begin{eqnarray*} \min_x & f(x) \\
          \text{subject to: } & c_j(x) =  0  ,  &j \in \mathcal{E}\\
            & c_j(x) \geq 0  ,  &j \in \mathcal{I}\\
           &  \text{lb}_i  \leq x_i \leq \text{ub}_i , &i = 1,...,N. \end{eqnarray*}

Where :math:`\mathcal{E}` or :math:`\mathcal{I}` are sets of indices
containing equality and inequality constraints.

Both linear and nonlinear constraints are defined as dictionaries with keys ``type``, ``fun`` and ``jac``.

    >>> ineq_cons = {'type': 'ineq',
    ...              'fun' : lambda x: np.array([1 - x[0] - 2*x[1],
    ...                                          1 - x[0]**2 - x[1],
    ...                                          1 - x[0]**2 + x[1]]),
    ...              'jac' : lambda x: np.array([[-1.0, -2.0],
    ...                                          [-2*x[0], -1.0],
    ...                                          [-2*x[0], 1.0]])}
    >>> eq_cons = {'type': 'eq',
    ...            'fun' : lambda x: np.array([2*x[0] + x[1] - 1]),
    ...            'jac' : lambda x: np.array([2.0, 1.0])}


And the optimization problem is solved with:

    >>> x0 = np.array([0.5, 0])
    >>> res = minimize(rosen, x0, method='SLSQP', jac=rosen_der,
    ...                constraints=[eq_cons, ineq_cons], options={'ftol': 1e-9, 'disp': True},
    ...                bounds=bounds)
    # may vary
    Optimization terminated successfully.    (Exit mode 0)
                Current function value: 0.342717574857755
                Iterations: 5
                Function evaluations: 6
                Gradient evaluations: 5
    >>> print(res.x)
    [0.41494475 0.1701105 ]

Most of the options available for the method ``'trust-constr'`` are not available
for ``'SLSQP'``.

Global optimization
-------------------

Global optimization aims to find the global minimum of a function within given
bounds, in the presence of potentially many local minima. Typically, global
minimizers efficiently search the parameter space, while using a local
minimizer (e.g., :func:`minimize`) under the hood.  SciPy contains a
number of good global optimizers.  Here, we'll use those on the same objective
function, namely the (aptly named) ``eggholder`` function::

   >>> def eggholder(x):
   ...     return (-(x[1] + 47) * np.sin(np.sqrt(abs(x[0]/2 + (x[1]  + 47))))
   ...             -x[0] * np.sin(np.sqrt(abs(x[0] - (x[1]  + 47)))))

   >>> bounds = [(-512, 512), (-512, 512)]

This function looks like an egg carton::

   >>> import matplotlib.pyplot as plt
   >>> from mpl_toolkits.mplot3d import Axes3D

   >>> x = np.arange(-512, 513)
   >>> y = np.arange(-512, 513)
   >>> xgrid, ygrid = np.meshgrid(x, y)
   >>> xy = np.stack([xgrid, ygrid])

   >>> fig = plt.figure()
   >>> ax = fig.add_subplot(111, projection='3d')
   >>> ax.view_init(45, -45)
   >>> ax.plot_surface(xgrid, ygrid, eggholder(xy), cmap='terrain')
   >>> ax.set_xlabel('x')
   >>> ax.set_ylabel('y')
   >>> ax.set_zlabel('eggholder(x, y)')
   >>> plt.show()

.. plot:: tutorial/examples/optimize_global_2.py
   :alt: "A 3-D plot shown from a three-quarter view. The function is very noisy with dozens of valleys and peaks. There is no clear min or max discernible from this view and it's not possible to see all the local peaks and valleys from this view."
   :align: center
   :include-source: 0

We now use the global optimizers to obtain the minimum and the function value
at the minimum. We'll store the results in a dictionary so we can compare
different optimization results later.

   >>> from scipy import optimize
   >>> results = dict()
   >>> results['shgo'] = optimize.shgo(eggholder, bounds)
   >>> results['shgo']
        fun: -935.3379515604197  # may vary
       funl: array([-935.33795156])
    message: 'Optimization terminated successfully.'
       nfev: 42
        nit: 2
      nlfev: 37
      nlhev: 0
      nljev: 9
    success: True
          x: array([439.48096952, 453.97740589])
         xl: array([[439.48096952, 453.97740589]])

   >>> results['DA'] = optimize.dual_annealing(eggholder, bounds)
   >>> results['DA']
        fun: -956.9182316237413  # may vary
    message: ['Maximum number of iteration reached']
       nfev: 4091
       nhev: 0
        nit: 1000
       njev: 0
          x: array([482.35324114, 432.87892901])

All optimizers return an ``OptimizeResult``, which in addition to the solution
contains information on the number of function evaluations, whether the
optimization was successful, and more.  For brevity, we won't show the full
output of the other optimizers::

   >>> results['DE'] = optimize.differential_evolution(eggholder, bounds)

:func:`shgo` has a second method, which returns all local minima rather than
only what it thinks is the global minimum::

   >>> results['shgo_sobol'] = optimize.shgo(eggholder, bounds, n=200, iters=5,
   ...                                       sampling_method='sobol')

We'll now plot all found minima on a heatmap of the function::

   >>> fig = plt.figure()
   >>> ax = fig.add_subplot(111)
   >>> im = ax.imshow(eggholder(xy), interpolation='bilinear', origin='lower',
   ...                cmap='gray')
   >>> ax.set_xlabel('x')
   >>> ax.set_ylabel('y')
   >>>
   >>> def plot_point(res, marker='o', color=None):
   ...     ax.plot(512+res.x[0], 512+res.x[1], marker=marker, color=color, ms=10)

   >>> plot_point(results['DE'], color='c')  # differential_evolution - cyan
   >>> plot_point(results['DA'], color='w')  # dual_annealing.        - white

   >>> # SHGO produces multiple minima, plot them all (with a smaller marker size)
   >>> plot_point(results['shgo'], color='r', marker='+')
   >>> plot_point(results['shgo_sobol'], color='r', marker='x')
   >>> for i in range(results['shgo_sobol'].xl.shape[0]):
   ...     ax.plot(512 + results['shgo_sobol'].xl[i, 0],
   ...             512 + results['shgo_sobol'].xl[i, 1],
   ...             'ro', ms=2)

   >>> ax.set_xlim([-4, 514*2])
   >>> ax.set_ylim([-4, 514*2])
   >>> plt.show()

.. plot:: tutorial/examples/optimize_global_1.py
   :align: center
   :alt: "This X-Y plot is a heatmap with the Z value denoted with the lowest points as black and the highest values as white. The image resembles a chess board rotated 45 degrees but heavily smoothed. A red dot is located at many of the minima on the grid resulting from the SHGO optimizer. SHGO shows the global minima as a red X in the top right. A local minima found with dual annealing is a white circle marker in the top left. A different local minima found with basinhopping is a yellow marker in the top center. The code is plotting the differential evolution result as a cyan circle, but it is not visible on the plot. At a glance it's not clear which of these valleys is the true global minima."
   :include-source: 0

Least-squares minimization (:func:`least_squares`)
--------------------------------------------------

SciPy is capable of solving robustified bound-constrained nonlinear
least-squares problems:

.. math::
   :nowrap:

   \begin{align}
   &\min_\mathbf{x} \frac{1}{2} \sum_{i = 1}^m \rho\left(f_i(\mathbf{x})^2\right) \\
   &\text{subject to }\mathbf{lb} \leq \mathbf{x} \leq \mathbf{ub}
   \end{align}

Here :math:`f_i(\mathbf{x})` are smooth functions from
:math:`\mathbb{R}^n` to :math:`\mathbb{R}`, we refer to them as residuals.
The purpose of a scalar-valued function :math:`\rho(\cdot)` is to reduce the
influence of outlier residuals and contribute to robustness of the solution,
we refer to it as a loss function. A linear loss function gives a standard
least-squares problem. Additionally, constraints in a form of lower and upper
bounds on some of :math:`x_j` are allowed.

All methods specific to least-squares minimization utilize a :math:`m \times n`
matrix of partial derivatives called Jacobian and defined as
:math:`J_{ij} = \partial f_i / \partial x_j`. It is highly recommended to
compute this matrix analytically and pass it to :func:`least_squares`,
otherwise, it will be estimated by finite differences, which takes a lot of
additional time and can be very inaccurate in hard cases.

Function :func:`least_squares` can be used for fitting a function
:math:`\varphi(t; \mathbf{x})` to empirical data :math:`\{(t_i, y_i), i = 0, \ldots, m-1\}`.
To do this, one should simply precompute residuals as
:math:`f_i(\mathbf{x}) = w_i (\varphi(t_i; \mathbf{x}) - y_i)`, where :math:`w_i`
are weights assigned to each observation.

Example of solving a fitting problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we consider an enzymatic reaction [1]_. There are 11 residuals defined as

.. math::
    f_i(x) = \frac{x_0 (u_i^2 + u_i x_1)}{u_i^2 + u_i x_2 + x_3} - y_i, \quad i = 0, \ldots, 10,

where :math:`y_i` are measurement values and :math:`u_i` are values of
the independent variable. The unknown vector of parameters is
:math:`\mathbf{x} = (x_0, x_1, x_2, x_3)^T`. As was said previously, it is
recommended to compute Jacobian matrix in a closed form:

.. math::
   :nowrap:

    \begin{align}
    &J_{i0} = \frac{\partial f_i}{\partial x_0} = \frac{u_i^2 + u_i x_1}{u_i^2 + u_i x_2 + x_3} \\
    &J_{i1} = \frac{\partial f_i}{\partial x_1} = \frac{u_i x_0}{u_i^2 + u_i x_2 + x_3} \\
    &J_{i2} = \frac{\partial f_i}{\partial x_2} = -\frac{x_0 (u_i^2 + u_i x_1) u_i}{(u_i^2 + u_i x_2 + x_3)^2} \\
    &J_{i3} = \frac{\partial f_i}{\partial x_3} = -\frac{x_0 (u_i^2 + u_i x_1)}{(u_i^2 + u_i x_2 + x_3)^2}
    \end{align}

We are going to use the "hard" starting point defined in [2]_. To find a
physically meaningful solution, avoid potential division by zero and assure
convergence to the global minimum we impose constraints
:math:`0 \leq x_j \leq 100, j = 0, 1, 2, 3`.

The code below implements least-squares estimation of :math:`\mathbf{x}` and
finally plots the original data and the fitted model function:

.. plot::
    :alt: "This code plots an X-Y time-series. The series starts in the lower left at (0, 0) and rapidly trends up to the maximum of 0.2 then flattens out. The fitted model is shown as a smooth orange trace and is well fit to the data."

    >>> from scipy.optimize import least_squares

    >>> def model(x, u):
    ...     return x[0] * (u ** 2 + x[1] * u) / (u ** 2 + x[2] * u + x[3])

    >>> def fun(x, u, y):
    ...     return model(x, u) - y

    >>> def jac(x, u, y):
    ...     J = np.empty((u.size, x.size))
    ...     den = u ** 2 + x[2] * u + x[3]
    ...     num = u ** 2 + x[1] * u
    ...     J[:, 0] = num / den
    ...     J[:, 1] = x[0] * u / den
    ...     J[:, 2] = -x[0] * num * u / den ** 2
    ...     J[:, 3] = -x[0] * num / den ** 2
    ...     return J

    >>> u = np.array([4.0, 2.0, 1.0, 5.0e-1, 2.5e-1, 1.67e-1, 1.25e-1, 1.0e-1,
    ...               8.33e-2, 7.14e-2, 6.25e-2])
    >>> y = np.array([1.957e-1, 1.947e-1, 1.735e-1, 1.6e-1, 8.44e-2, 6.27e-2,
    ...               4.56e-2, 3.42e-2, 3.23e-2, 2.35e-2, 2.46e-2])
    >>> x0 = np.array([2.5, 3.9, 4.15, 3.9])
    >>> res = least_squares(fun, x0, jac=jac, bounds=(0, 100), args=(u, y), verbose=1)
    # may vary
    `ftol` termination condition is satisfied.
    Function evaluations 130, initial cost 4.4383e+00, final cost 1.5375e-04, first-order optimality 4.92e-08.
    >>> res.x
    array([ 0.19280596,  0.19130423,  0.12306063,  0.13607247])

    >>> import matplotlib.pyplot as plt
    >>> u_test = np.linspace(0, 5)
    >>> y_test = model(res.x, u_test)
    >>> plt.plot(u, y, 'o', markersize=4, label='data')
    >>> plt.plot(u_test, y_test, label='fitted model')
    >>> plt.xlabel("u")
    >>> plt.ylabel("y")
    >>> plt.legend(loc='lower right')
    >>> plt.show()

.. [1] J. Kowalik and J. F. Morrison, “Analysis of kinetic data for allosteric enzyme reactions as
   a nonlinear regression problem”, Math. Biosci., vol. 2, pp. 57-66, 1968.
.. [2] B. M. Averick et al., “The MINPACK-2 Test Problem Collection”.

Further examples
^^^^^^^^^^^^^^^^

Three interactive examples below illustrate usage of :func:`least_squares` in
greater detail.

1. `Large-scale bundle adjustment in scipy <https://scipy-cookbook.readthedocs.io/items/bundle_adjustment.html>`_
   demonstrates large-scale capabilities of :func:`least_squares` and how to
   efficiently compute finite difference approximation of sparse Jacobian.
2. `Robust nonlinear regression in scipy <https://scipy-cookbook.readthedocs.io/items/robust_regression.html>`_
   shows how to handle outliers with a robust loss function in a nonlinear
   regression.
3. `Solving a discrete boundary-value problem in scipy <https://scipy-cookbook.readthedocs.io/items/discrete_bvp.html>`_
   examines how to solve a large system of equations and use bounds to achieve
   desired properties of the solution.

For the details about mathematical algorithms behind the implementation refer
to documentation of :func:`least_squares`.


Univariate function minimizers (:func:`minimize_scalar`)
--------------------------------------------------------

Often only the minimum of an univariate function (i.e., a function that
takes a scalar as input) is needed. In these circumstances, other
optimization techniques have been developed that can work faster. These are
accessible from the :func:`minimize_scalar` function, which proposes several
algorithms.


Unconstrained minimization (``method='brent'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are, actually, two methods that can be used to minimize an univariate
function: `brent` and `golden`, but `golden` is included only for academic
purposes and should rarely be used. These can be respectively selected
through the `method` parameter in :func:`minimize_scalar`. The `brent`
method uses Brent's algorithm for locating a minimum. Optimally, a bracket
(the `bracket` parameter) should be given which contains the minimum desired. A
bracket is a triple :math:`\left( a, b, c \right)` such that :math:`f
\left( a \right) > f \left( b \right) < f \left( c \right)` and :math:`a <
b < c` . If this is not given, then alternatively two starting points can
be chosen and a bracket will be found from these points using a simple
marching algorithm. If these two starting points are not provided, `0` and
`1` will be used (this may not be the right choice for your function and
result in an unexpected minimum being returned).

Here is an example:

    >>> from scipy.optimize import minimize_scalar
    >>> f = lambda x: (x - 2) * (x + 1)**2
    >>> res = minimize_scalar(f, method='brent')
    >>> print(res.x)
    1.0


Bounded minimization (``method='bounded'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Very often, there are constraints that can be placed on the solution space
before minimization occurs. The `bounded` method in :func:`minimize_scalar`
is an example of a constrained minimization procedure that provides a
rudimentary interval constraint for scalar functions. The interval
constraint allows the minimization to occur only between two fixed
endpoints, specified using the mandatory `bounds` parameter.

For example, to find the minimum of :math:`J_{1}\left( x \right)` near
:math:`x=5` , :func:`minimize_scalar` can be called using the interval
:math:`\left[ 4, 7 \right]` as a constraint. The result is
:math:`x_{\textrm{min}}=5.3314` :

    >>> from scipy.special import j1
    >>> res = minimize_scalar(j1, bounds=(4, 7), method='bounded')
    >>> res.x
    5.33144184241


Custom minimizers
-----------------

Sometimes, it may be useful to use a custom method as a (multivariate
or univariate) minimizer, for example, when using some library wrappers
of :func:`minimize` (e.g., :func:`basinhopping`).

We can achieve that by, instead of passing a method name, passing
a callable (either a function or an object implementing a `__call__`
method) as the `method` parameter.

Let us consider an (admittedly rather virtual) need to use a trivial
custom multivariate minimization method that will just search the
neighborhood in each dimension independently with a fixed step size::

    >>> from scipy.optimize import OptimizeResult
    >>> def custmin(fun, x0, args=(), maxfev=None, stepsize=0.1,
    ...         maxiter=100, callback=None, **options):
    ...     bestx = x0
    ...     besty = fun(x0)
    ...     funcalls = 1
    ...     niter = 0
    ...     improved = True
    ...     stop = False
    ...
    ...     while improved and not stop and niter < maxiter:
    ...         improved = False
    ...         niter += 1
    ...         for dim in range(np.size(x0)):
    ...             for s in [bestx[dim] - stepsize, bestx[dim] + stepsize]:
    ...                 testx = np.copy(bestx)
    ...                 testx[dim] = s
    ...                 testy = fun(testx, *args)
    ...                 funcalls += 1
    ...                 if testy < besty:
    ...                     besty = testy
    ...                     bestx = testx
    ...                     improved = True
    ...             if callback is not None:
    ...                 callback(bestx)
    ...             if maxfev is not None and funcalls >= maxfev:
    ...                 stop = True
    ...                 break
    ...
    ...     return OptimizeResult(fun=besty, x=bestx, nit=niter,
    ...                           nfev=funcalls, success=(niter > 1))
    >>> x0 = [1.35, 0.9, 0.8, 1.1, 1.2]
    >>> res = minimize(rosen, x0, method=custmin, options=dict(stepsize=0.05))
    >>> res.x
    array([1., 1., 1., 1., 1.])

This will work just as well in case of univariate optimization::

    >>> def custmin(fun, bracket, args=(), maxfev=None, stepsize=0.1,
    ...         maxiter=100, callback=None, **options):
    ...     bestx = (bracket[1] + bracket[0]) / 2.0
    ...     besty = fun(bestx)
    ...     funcalls = 1
    ...     niter = 0
    ...     improved = True
    ...     stop = False
    ...
    ...     while improved and not stop and niter < maxiter:
    ...         improved = False
    ...         niter += 1
    ...         for testx in [bestx - stepsize, bestx + stepsize]:
    ...             testy = fun(testx, *args)
    ...             funcalls += 1
    ...             if testy < besty:
    ...                 besty = testy
    ...                 bestx = testx
    ...                 improved = True
    ...         if callback is not None:
    ...             callback(bestx)
    ...         if maxfev is not None and funcalls >= maxfev:
    ...             stop = True
    ...             break
    ...
    ...     return OptimizeResult(fun=besty, x=bestx, nit=niter,
    ...                           nfev=funcalls, success=(niter > 1))
    >>> def f(x):
    ...    return (x - 2)**2 * (x + 2)**2
    >>> res = minimize_scalar(f, bracket=(-3.5, 0), method=custmin,
    ...                       options=dict(stepsize = 0.05))
    >>> res.x
    -2.0


Root finding
------------

Scalar functions
^^^^^^^^^^^^^^^^

If one has a single-variable equation, there are multiple different root
finding algorithms that can be tried. Most of these algorithms require the
endpoints of an interval in which a root is expected (because the function
changes signs). In general, :obj:`brentq` is the best choice, but the other
methods may be useful in certain circumstances or for academic purposes.
When a bracket is not available, but one or more derivatives are available,
then :obj:`newton` (or ``halley``, ``secant``) may be applicable.
This is especially the case if the function is defined on a subset of the
complex plane, and the bracketing methods cannot be used.


Fixed-point solving
^^^^^^^^^^^^^^^^^^^

A problem closely related to finding the zeros of a function is the
problem of finding a fixed point of a function. A fixed point of a
function is the point at which evaluation of the function returns the
point: :math:`g\left(x\right)=x.` Clearly, the fixed point of :math:`g`
is the root of :math:`f\left(x\right)=g\left(x\right)-x.`
Equivalently, the root of :math:`f` is the fixed point of
:math:`g\left(x\right)=f\left(x\right)+x.` The routine
:obj:`fixed_point` provides a simple iterative method using Aitkens
sequence acceleration to estimate the fixed point of :math:`g` given a
starting point.

Sets of equations
^^^^^^^^^^^^^^^^^

Finding a root of a set of non-linear equations can be achieved using the
:func:`root` function. Several methods are available, amongst which ``hybr``
(the default) and ``lm``, which, respectively, use the hybrid method of Powell
and the Levenberg-Marquardt method from MINPACK.

The following example considers the single-variable transcendental
equation

.. math::

    x+2\cos\left(x\right)=0,

a root of which can be found as follows::

    >>> import numpy as np
    >>> from scipy.optimize import root
    >>> def func(x):
    ...     return x + 2 * np.cos(x)
    >>> sol = root(func, 0.3)
    >>> sol.x
    array([-1.02986653])
    >>> sol.fun
    array([ -6.66133815e-16])

Consider now a set of non-linear equations

.. math::
   :nowrap:

    \begin{eqnarray*}
    x_{0}\cos\left(x_{1}\right) & = & 4,\\
    x_{0}x_{1}-x_{1} & = & 5.
    \end{eqnarray*}

We define the objective function so that it also returns the Jacobian and
indicate this by setting the ``jac`` parameter to ``True``. Also, the
Levenberg-Marquardt solver is used here.

::

    >>> def func2(x):
    ...     f = [x[0] * np.cos(x[1]) - 4,
    ...          x[1]*x[0] - x[1] - 5]
    ...     df = np.array([[np.cos(x[1]), -x[0] * np.sin(x[1])],
    ...                    [x[1], x[0] - 1]])
    ...     return f, df
    >>> sol = root(func2, [1, 1], jac=True, method='lm')
    >>> sol.x
    array([ 6.50409711,  0.90841421])


Root finding for large problems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Methods ``hybr`` and ``lm`` in :func:`root` cannot deal with a very large
number of variables (*N*), as they need to calculate and invert a dense *N
x N* Jacobian matrix on every Newton step. This becomes rather inefficient
when *N* grows.

Consider, for instance, the following problem: we need to solve the
following integrodifferential equation on the square
:math:`[0,1]\times[0,1]`:

.. math::

   (\partial_x^2 + \partial_y^2) P + 5 \left(\int_0^1\int_0^1\cosh(P)\,dx\,dy\right)^2 = 0

with the boundary condition :math:`P(x,1) = 1` on the upper edge and
:math:`P=0` elsewhere on the boundary of the square. This can be done
by approximating the continuous function *P* by its values on a grid,
:math:`P_{n,m}\approx{}P(n h, m h)`, with a small grid spacing
*h*. The derivatives and integrals can then be approximated; for
instance :math:`\partial_x^2 P(x,y)\approx{}(P(x+h,y) - 2 P(x,y) +
P(x-h,y))/h^2`. The problem is then equivalent to finding the root of
some function ``residual(P)``, where ``P`` is a vector of length
:math:`N_x N_y`.

Now, because :math:`N_x N_y` can be large, methods ``hybr`` or ``lm`` in
:func:`root` will take a long time to solve this problem. The solution can,
however, be found using one of the large-scale solvers, for example
``krylov``, ``broyden2``, or ``anderson``. These use what is known as the
inexact Newton method, which instead of computing the Jacobian matrix
exactly, forms an approximation for it.

The problem we have can now be solved as follows:

.. plot::
    :alt: "This code generates a 2-D heatmap with Z values from 0 to 1. The graph resembles a smooth, dark blue-green, U shape, with an open yellow top. The right, bottom, and left edges have a value near zero and the top has a value close to 1. The center of the solution space has a value close to 0.8."

    import numpy as np
    from scipy.optimize import root
    from numpy import cosh, zeros_like, mgrid, zeros

    # parameters
    nx, ny = 75, 75
    hx, hy = 1./(nx-1), 1./(ny-1)

    P_left, P_right = 0, 0
    P_top, P_bottom = 1, 0

    def residual(P):
       d2x = zeros_like(P)
       d2y = zeros_like(P)

       d2x[1:-1] = (P[2:]   - 2*P[1:-1] + P[:-2]) / hx/hx
       d2x[0]    = (P[1]    - 2*P[0]    + P_left)/hx/hx
       d2x[-1]   = (P_right - 2*P[-1]   + P[-2])/hx/hx

       d2y[:,1:-1] = (P[:,2:] - 2*P[:,1:-1] + P[:,:-2])/hy/hy
       d2y[:,0]    = (P[:,1]  - 2*P[:,0]    + P_bottom)/hy/hy
       d2y[:,-1]   = (P_top   - 2*P[:,-1]   + P[:,-2])/hy/hy

       return d2x + d2y + 5*cosh(P).mean()**2

    # solve
    guess = zeros((nx, ny), float)
    sol = root(residual, guess, method='krylov', options={'disp': True})
    #sol = root(residual, guess, method='broyden2', options={'disp': True, 'max_rank': 50})
    #sol = root(residual, guess, method='anderson', options={'disp': True, 'M': 10})
    print('Residual: %g' % abs(residual(sol.x)).max())

    # visualize
    import matplotlib.pyplot as plt
    x, y = mgrid[0:1:(nx*1j), 0:1:(ny*1j)]
    plt.pcolormesh(x, y, sol.x, shading='gouraud')
    plt.colorbar()
    plt.show()


Still too slow? Preconditioning.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When looking for the zero of the functions :math:`f_i({\bf x}) = 0`,
*i = 1, 2, ..., N*, the ``krylov`` solver spends most of the
time inverting the Jacobian matrix,

.. math:: J_{ij} = \frac{\partial f_i}{\partial x_j} .

If you have an approximation for the inverse matrix
:math:`M\approx{}J^{-1}`, you can use it for *preconditioning* the
linear-inversion problem. The idea is that instead of solving
:math:`J{\bf s}={\bf y}` one solves :math:`MJ{\bf s}=M{\bf y}`: since
matrix :math:`MJ` is "closer" to the identity matrix than :math:`J`
is, the equation should be easier for the Krylov method to deal with.

The matrix *M* can be passed to :func:`root` with method ``krylov`` as an
option ``options['jac_options']['inner_M']``. It can be a (sparse) matrix
or a :obj:`scipy.sparse.linalg.LinearOperator` instance.

For the problem in the previous section, we note that the function to
solve consists of two parts: the first one is the application of the
Laplace operator, :math:`[\partial_x^2 + \partial_y^2] P`, and the second
is the integral. We can actually easily compute the Jacobian corresponding
to the Laplace operator part: we know that in 1-D

.. math::

   \partial_x^2 \approx \frac{1}{h_x^2} \begin{pmatrix}
   -2 & 1 & 0 & 0 \cdots \\
   1 & -2 & 1 & 0 \cdots \\
   0 & 1 & -2 & 1 \cdots \\
   \ldots
   \end{pmatrix}
   = h_x^{-2} L

so that the whole 2-D operator is represented by

.. math::

   J_1 = \partial_x^2 + \partial_y^2
   \simeq
   h_x^{-2} L \otimes I + h_y^{-2} I \otimes L

The matrix :math:`J_2` of the Jacobian corresponding to the integral
is more difficult to calculate, and since *all* of it entries are
nonzero, it will be difficult to invert. :math:`J_1` on the other hand
is a relatively simple matrix, and can be inverted by
:obj:`scipy.sparse.linalg.splu` (or the inverse can be approximated by
:obj:`scipy.sparse.linalg.spilu`). So we are content to take
:math:`M\approx{}J_1^{-1}` and hope for the best.

In the example below, we use the preconditioner :math:`M=J_1^{-1}`.

.. literalinclude:: examples/newton_krylov_preconditioning.py

Resulting run, first without preconditioning::

  0:  |F(x)| = 803.614; step 1; tol 0.000257947
  1:  |F(x)| = 345.912; step 1; tol 0.166755
  2:  |F(x)| = 139.159; step 1; tol 0.145657
  3:  |F(x)| = 27.3682; step 1; tol 0.0348109
  4:  |F(x)| = 1.03303; step 1; tol 0.00128227
  5:  |F(x)| = 0.0406634; step 1; tol 0.00139451
  6:  |F(x)| = 0.00344341; step 1; tol 0.00645373
  7:  |F(x)| = 0.000153671; step 1; tol 0.00179246
  8:  |F(x)| = 6.7424e-06; step 1; tol 0.00173256
  Residual 3.57078908664e-07
  Evaluations 317

and then with preconditioning::

  0:  |F(x)| = 136.993; step 1; tol 7.49599e-06
  1:  |F(x)| = 4.80983; step 1; tol 0.00110945
  2:  |F(x)| = 0.195942; step 1; tol 0.00149362
  3:  |F(x)| = 0.000563597; step 1; tol 7.44604e-06
  4:  |F(x)| = 1.00698e-09; step 1; tol 2.87308e-12
  Residual 9.29603061195e-11
  Evaluations 77

Using a preconditioner reduced the number of evaluations of the
``residual`` function by a factor of *4*. For problems where the
residual is expensive to compute, good preconditioning can be crucial
--- it can even decide whether the problem is solvable in practice or
not.

Preconditioning is an art, science, and industry. Here, we were lucky
in making a simple choice that worked reasonably well, but there is a
lot more depth to this topic than is shown here.

Linear programming (:func:`linprog`)
------------------------------------

The function :func:`linprog` can minimize a linear objective function
subject to linear equality and inequality constraints. This kind of
problem is well known as linear programming. Linear programming solves
problems of the following form:

.. math::

        \min_x \ & c^T x \\
        \mbox{such that} \ & A_{ub} x \leq b_{ub},\\
        & A_{eq} x = b_{eq},\\
        & l \leq x \leq u ,

where :math:`x` is a vector of decision variables; :math:`c`, :math:`b_{ub}`,
:math:`b_{eq}`, :math:`l`, and :math:`u` are vectors; and :math:`A_{ub}` and
:math:`A_{eq}` are matrices.

In this tutorial, we will try to solve a typical linear programming
problem using :func:`linprog`.

Linear programming example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider the following simple linear programming problem:

.. math::
        \max_{x_1, x_2, x_3, x_4} \ & 29x_1 + 45x_2 \\
        \mbox{such that} \
        & x_1 -x_2 -3x_3 \leq 5\\
        & 2x_1 -3x_2 -7x_3 + 3x_4 \geq 10\\
        & 2x_1 + 8x_2 + x_3 = 60\\
        & 4x_1 + 4x_2 + x_4 = 60\\
        & 0 \leq x_0\\
        & 0 \leq x_1 \leq 5\\
        & x_2 \leq 0.5\\
        & -3 \leq x_3\\

We need some mathematical manipulations to convert the target problem to the form accepted by :func:`linprog`.

First of all, let's consider the objective function.
We want to maximize the objective
function, but :func:`linprog` can only accept a minimization problem. This is easily remedied by converting the maximize
:math:`29x_1 + 45x_2` to minimizing :math:`-29x_1 -45x_2`. Also, :math:`x_3, x_4` are not shown in the objective
function. That means the weights corresponding with :math:`x_3, x_4` are zero. So, the objective function can be
converted to:

.. math::
        \min_{x_1, x_2, x_3, x_4} \ -29x_1 -45x_2 + 0x_3 + 0x_4

If we define the vector of decision variables :math:`x = [x_1, x_2, x_3, x_4]^T`, the objective weights vector :math:`c` of :func:`linprog` in this problem
should be

.. math::
        c = [-29, -45, 0, 0]^T

Next, let's consider the two inequality constraints. The first one is a "less than" inequality, so it is already in the form accepted by `linprog`.
The second one is a "greater than" inequality, so we need to multiply both sides by :math:`-1` to convert it to a "less than" inequality.
Explicitly showing zero coefficients, we have:

.. math::
        x_1 -x_2 -3x_3 + 0x_4  &\leq 5\\
        -2x_1 + 3x_2 + 7x_3 - 3x_4 &\leq -10\\

These equations can be converted to matrix form:

.. math::
    A_{ub} x \leq b_{ub}\\

where

.. math::
   :nowrap:

    \begin{equation*} A_{ub} =
    \begin{bmatrix} 1 & -1 & -3 & 0 \\
                    -2 & 3 & 7 & -3
    \end{bmatrix}
    \end{equation*}

.. math::
   :nowrap:

    \begin{equation*} b_{ub} =
    \begin{bmatrix} 5 \\
                    -10
    \end{bmatrix}
    \end{equation*}

Next, let's consider the two equality constraints. Showing zero weights explicitly, these are:

.. math::
        2x_1 + 8x_2 + 1x_3 + 0x_4 &= 60\\
        4x_1 + 4x_2 + 0x_3 + 1x_4 &= 60\\

These equations can be converted to matrix form:

.. math::
    A_{eq} x = b_{eq}\\

where

.. math::
   :nowrap:

    \begin{equation*} A_{eq} =
    \begin{bmatrix} 2 & 8 & 1 & 0 \\
                    4 & 4 & 0 & 1
    \end{bmatrix}
    \end{equation*}

.. math::
   :nowrap:

    \begin{equation*} b_{eq} =
    \begin{bmatrix} 60 \\
                    60
    \end{bmatrix}
    \end{equation*}

Lastly, let's consider the separate inequality constraints on individual decision variables, which are known as
"box constraints" or "simple bounds". These constraints can be applied using the bounds argument of :func:`linprog`.
As noted in the :func:`linprog` documentation, the default value of bounds is ``(0, None)``, meaning that the
lower bound on each decision variable is 0, and the upper bound on each decision variable is infinity:
all the decision variables are non-negative. Our bounds are different, so we will need to specify the lower and upper bound on each
decision variable as a tuple and group these tuples into a list.


Finally, we can solve the transformed problem using :func:`linprog`.

::

    >>> import numpy as np
    >>> from scipy.optimize import linprog
    >>> c = np.array([-29.0, -45.0, 0.0, 0.0])
    >>> A_ub = np.array([[1.0, -1.0, -3.0, 0.0],
    ...                 [-2.0, 3.0, 7.0, -3.0]])
    >>> b_ub = np.array([5.0, -10.0])
    >>> A_eq = np.array([[2.0, 8.0, 1.0, 0.0],
    ...                 [4.0, 4.0, 0.0, 1.0]])
    >>> b_eq = np.array([60.0, 60.0])
    >>> x0_bounds = (0, None)
    >>> x1_bounds = (0, 5.0)
    >>> x2_bounds = (-np.inf, 0.5)  # +/- np.inf can be used instead of None
    >>> x3_bounds = (-3.0, None)
    >>> bounds = [x0_bounds, x1_bounds, x2_bounds, x3_bounds]
    >>> result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds)
    >>> print(result.message)
    The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is At lower/fixed bound)

The result states that our problem is infeasible, meaning that there is no solution vector that satisfies all the
constraints. That doesn't necessarily mean we did anything wrong; some problems truly are infeasible.
Suppose, however, that we were to decide that our bound constraint on :math:`x_1` was too tight and that it could be loosened
to :math:`0 \leq x_1 \leq 6`. After adjusting our code ``x1_bounds = (0, 6)`` to reflect the change and executing it again:

::

    >>> x1_bounds = (0, 6)
    >>> bounds = [x0_bounds, x1_bounds, x2_bounds, x3_bounds]
    >>> result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds)
    >>> print(result.message)
    Optimization terminated successfully. (HiGHS Status 7: Optimal)

The result shows the optimization was successful.
We can check the objective value (``result.fun``) is same as :math:`c^Tx`:

::

    >>> x = np.array(result.x)
    >>> obj = result.fun
    >>> print(c @ x)
    -505.97435889013434  # may vary
    >>> print(obj)
    -505.97435889013434  # may vary

We can also check that all constraints are satisfied within reasonable tolerances:

::

    >>> print(b_ub - (A_ub @ x).flatten())  # this is equivalent to result.slack
    [ 6.52747190e-10, -2.26730279e-09]  # may vary
    >>> print(b_eq - (A_eq @ x).flatten())  # this is equivalent to result.con
    [ 9.78840831e-09, 1.04662945e-08]]  # may vary
    >>> print([0 <= result.x[0], 0 <= result.x[1] <= 6.0, result.x[2] <= 0.5, -3.0 <= result.x[3]])
    [True, True, True, True]


Assignment problems
-------------------

Linear sum assignment problem example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider the problem of selecting students for a swimming medley relay team.
We have a table showing times for each swimming style of five students:

==========  ===========  ============  ===========  ===============================
 Student    backstroke   breaststroke  butterfly    freestyle
==========  ===========  ============  ===========  ===============================
 A          43.5           47.1         48.4        38.2
 B          45.5           42.1         49.6        36.8
 C          43.4           39.1         42.1        43.2
 D          46.5           44.1         44.5        41.2
 E          46.3           47.8         50.4        37.2
==========  ===========  ============  ===========  ===============================

We need to choose a student for each of the four swimming styles such that
the total relay time is minimized.
This is a typical linear sum assignment problem. We can use :func:`linear_sum_assignment` to solve it.

The linear sum assignment problem is one of the most famous combinatorial optimization problems.
Given a "cost matrix" :math:`C`, the problem is to choose

- exactly one element from each row
- without choosing more than one element from any column
- such that the sum of the chosen elements is minimized

In other words, we need to assign each row to one column such that the sum of
the corresponding entries is minimized.

Formally, let :math:`X` be a boolean matrix where :math:`X[i,j] = 1` iff row  :math:`i` is assigned to column :math:`j`.
Then the optimal assignment has cost

.. math::

    \min \sum_i \sum_j C_{i,j} X_{i,j}

The first step is to define the cost matrix.
In this example, we want to assign each swimming style to a student.
:func:`linear_sum_assignment` is able to assign each row of a cost matrix to a column.
Therefore, to form the cost matrix, the table above needs to be transposed so that the rows
correspond with swimming styles and the columns correspond with students:

::

    >>> import numpy as np
    >>> cost = np.array([[43.5, 45.5, 43.4, 46.5, 46.3],
    ...                  [47.1, 42.1, 39.1, 44.1, 47.8],
    ...                  [48.4, 49.6, 42.1, 44.5, 50.4],
    ...                  [38.2, 36.8, 43.2, 41.2, 37.2]])

We can solve the assignment problem with :func:`linear_sum_assignment`:

::

    >>> from scipy.optimize import linear_sum_assignment
    >>> row_ind, col_ind = linear_sum_assignment(cost)

The ``row_ind`` and ``col_ind`` are optimal assigned matrix indexes of the cost matrix:

::

    >>> row_ind
    array([0, 1, 2, 3])
    >>> col_ind
    array([0, 2, 3, 1])

The optimal assignment is:

::

    >>> styles = np.array(["backstroke", "breaststroke", "butterfly", "freestyle"])[row_ind]
    >>> students = np.array(["A", "B", "C", "D", "E"])[col_ind]
    >>> dict(zip(styles, students))
    {'backstroke': 'A', 'breaststroke': 'C', 'butterfly': 'D', 'freestyle': 'B'}

The optimal total medley time is:

::

    >>> cost[row_ind, col_ind].sum()
    163.89999999999998

Note that this result is not the same as the sum of the minimum times for each swimming style:

::

    >>> np.min(cost, axis=1).sum()
    161.39999999999998

because student "C" is the best swimmer in both "breaststroke" and "butterfly" style.
We cannot assign student "C" to both styles, so we assigned student C to the "breaststroke" style
and D to the "butterfly" style to minimize the total time.

.. rubric:: References

Some further reading and related software, such as Newton-Krylov [KK]_,
PETSc [PP]_, and PyAMG [AMG]_:

.. [KK] D.A. Knoll and D.E. Keyes, "Jacobian-free Newton-Krylov methods",
        J. Comp. Phys. 193, 357 (2004). :doi:`10.1016/j.jcp.2003.08.010`

.. [PP] PETSc https://www.mcs.anl.gov/petsc/ and its Python bindings
        https://bitbucket.org/petsc/petsc4py/

.. [AMG] PyAMG (algebraic multigrid preconditioners/solvers)
         https://github.com/pyamg/pyamg/issues

.. _tutorial-optimize_milp:

Mixed integer linear programming
---------------------------------

Knapsack problem example
^^^^^^^^^^^^^^^^^^^^^^^^^

The knapsack problem is a well known combinatorial optimization problem.
Given a set of items, each with a size and a value, the problem is to choose
the items that maximize the total value under the condition that the total size
is below a certain threshold.

Formally, let

- :math:`x_i` be a boolean variable that indicates whether item :math:`i` is
  included in the knapsack,

- :math:`n` be the total number of items,

- :math:`v_i` be the value of item :math:`i`,

- :math:`s_i` be the size of item :math:`i`, and

- :math:`C` be the capacity of the knapsack.

Then the problem is:

.. math::

    \max \sum_i^n  v_{i} x_{i}

.. math::

    \text{subject to} \sum_i^n s_{i} x_{i} \leq C,  x_{i} \in {0, 1}

Although the objective function and inequality constraints are linear in the
*decision variables* :math:`x_i`, this differs from a typical linear
programming problem in that the decision variables can only assume integer
values.  Specifically, our decision variables can only be :math:`0` or
:math:`1`, so this is known as a *binary integer linear program* (BILP). Such
a problem falls within the larger class of *mixed integer linear programs*
(MILPs), which we we can solve with :func:`milp`.

In our example, there are 8 items to choose from, and the size and value of
each is specified as follows.

::

    >>> import numpy as np
    >>> from scipy import optimize
    >>> sizes = np.array([21, 11, 15, 9, 34, 25, 41, 52])
    >>> values = np.array([22, 12, 16, 10, 35, 26, 42, 53])

We need to constrain our eight decision variables to be binary. We do so
by adding a :class:`Bounds`: constraint to ensure that they lie between
:math:`0` and :math:`1`, and we apply "integrality" constraints to ensure that
they are *either* :math:`0` *or* :math:`1`.

::

    >>> bounds = optimize.Bounds(0, 1)  # 0 <= x_i <= 1
    >>> integrality = np.full_like(values, True)  # x_i are integers

The knapsack capacity constraint is specified using :class:`LinearConstraint`.

::

    >>> capacity = 100
    >>> constraints = optimize.LinearConstraint(A=sizes, lb=0, ub=capacity)

If we are following the usual rules of linear algebra, the input ``A`` should
be a  two-dimensional matrix, and the lower and upper bounds ``lb`` and ``ub``
should be one-dimensional vectors, but :class:`LinearConstraint` is forgiving
as long as the inputs can be broadcast to consistent shapes.

Using the variables defined above, we can solve the knapsack problem using
:func:`milp`. Note that :func:`milp` minimizes the objective function, but we
want to maximize the total value, so we set `c` to be negative of the values.

::

    >>> from scipy.optimize import milp
    >>> res = milp(c=-values, constraints=constraints,
    ...            integrality=integrality, bounds=bounds)

Let's check the result:

::

    >>> res.success
    True
    >>> res.x
    array([1., 1., 0., 1., 1., 1., 0., 0.])

This means that we should select the items 1, 2, 4, 5, 6 to optimize the total
value under the size constraint. Note that this is different from we would have
obtained had we solved the *linear programming relaxation* (without integrality
constraints) and attempted to round the decision variables.

::

    >>> from scipy.optimize import milp
    >>> res = milp(c=-values, constraints=constraints,
    ...            integrality=False, bounds=bounds)
    >>> res.x
    array([1.        , 1.        , 1.        , 1.        ,
           0.55882353, 1.        , 0.        , 0.        ])

If we were to round this solution up to
``array([1., 1., 1., 1., 1., 1., 0., 0.])``, our knapsack would be over the
capacity constraint, whereas if we were to round down to
``array([1., 1., 1., 1., 0., 1., 0., 0.])``, we would have a sub-optimal
solution.

For more MILP tutorials, see the Jupyter notebooks on SciPy Cookbooks:

- `Compressed Sensing l1 program <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/LinearAndMixedIntegerLinearProgramming/compressed_sensing_milp_tutorial_1.ipynb>`_
- `Compressed Sensing l0 program <https://nbviewer.org/github/scipy/scipy-cookbook/blob/main/ipython/LinearAndMixedIntegerLinearProgramming/compressed_sensing_milp_tutorial_2.ipynb>`_
