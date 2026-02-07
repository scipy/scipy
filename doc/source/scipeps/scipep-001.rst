
+----------+------------------------------------------------------------------------+
| PEP:     |  1                                                                     |
+----------+------------------------------------------------------------------------+
| Title:   | Introduction of Optimizer and Function classes for scalar minimisation |
+----------+------------------------------------------------------------------------+
| Authors: | Andrew Nelson (andyfaff@gmail.com),                                    |
|          | Scott Sievert (scott@stsievert.com)                                    |
+----------+------------------------------------------------------------------------+
| Status:  | Draft                                                                  |
+----------+------------------------------------------------------------------------+
| Created: | Feb 2018                                                               |
+----------+------------------------------------------------------------------------+


*Abstract*---We propose rewriting the minimization API in SciPy in terms of
classes. We believe this will provide benefits to users and enhance the
maintainability of SciPy, and any changes will be backwards compatible.
Preliminary speed tests do not indicate a significant slowdown.

Motivation
==========

We propose rewriting the optimizer framework in class form. The new framework
will:

- be easier to maintain, test, and develop.
- preserve backwards compatibility, and not significantly degrade
  performance
- expose a new API to interact with an optimization, and easily create new
  optimizers

Optimization is a critical part of SciPy's API, and we believe we can enhance
this API while retaining backwards compatibility and improving maintainability.

``Optimizer`` and ``Function`` classes will be introduced. Backwards
compatibility for existing scalar minimizer functions will remain, but will
have their core-functionality carried out by a corresponding ``Optimizer``
(e.g., the existing ``fmin`` function will use the new ``NelderMead`` class).
Usage of ``optimize.minimize`` will be unchanged.

The ``Function`` class is responsible for calculating the function, gradient
and Hessian (and will implement numerical differentiation if gradient/Hessian
implementation not provided). The ``Function`` class is general and can be used
to map between arbitrary dimensions, including scalar and vector functions.

The ``Optimizer`` class is used to optimize a ``Function``. Nearly all
optimizers have some fundamental iterative behavior. As such, the ``Optimizer``
will be iterable, which allows stepwise progression through the problem via
``Optimizer.__next__``. At each iteration the solution state is available to
the user, which can be used for many purposes including: user defined halting
criteria, modification of solver hyper-parameters, tracking solution
trajectories, etc. Running the optimizer to completion is achievable with the
``solve`` or ``__call__`` methods.

Different optimization algorithms can inherit from ``Optimizer``, with each of
the subclasses overriding the ``__next__`` method to represent the core of their
iterative technique. For some solvers, each iteration is implemented in
C/Fortran with the main optimization loop in Python (e.g., LBFGSB). We are not
proposing to replace those external calls at this time.

Other optimizers run the complete optimization in external C/Fortran code
(e.g., ``leastsq`` which calls ``minpack``). These methods can run the entire
optimization in external code. Future work may involve performing a single
optimization step in C/Fortran by extraction of iteration logic into Python or
Cython, but is beyond the scope of this proposal (and would also require
rigorous benchmarking and testing).

The proposed changes will be transparent to an end-user of ``minimize`` and
``fmin``. The intermediate or advanced user will appreciate the ``Optimizer``
and ``Function`` classes.

Timeline
--------

1. The Optimizer, Function classes are privately introduced, with subclasses
   for some minimizers (NelderMead, LBFGS, BFGS, DifferentialEvolution) added.
   They will be used as the core functionality of ``fmin``, etc.
2. Subsequent private classes for the remaining minimizers will be created.
   This will allow fine tuning of the Optimizer and Function classes.
3. In later releases the ``Optimizer`` and ``Function`` classes will be
   available in the SciPy public API.

Rationale
=========

Simplified maintenance
----------------------

The maintenance burden of the new classes will be significantly reduced compared
to the current state of scipy.optimize. It will be easier to develop new
features and provide more comprehensive testing.

The main reason for this is class inheritance. Improvements made to the base
``Optimizer`` class mean that all that all inheriting objects improve. Currently
such changes have to be made in each minimizer, which leads to code duplication,
and the attendant risk of bugs being introduced.

For example:

* Placing numerical differentiation in the Function class

  * would allow either absolute or relative delta change to be made easily,
    and in one place.

  * and, with the current API, require modifications and extra keywords for
    all minimizer functions.

* The user wishes to halt optimization early (issues `#4384
  <https://github.com/scipy/scipy/issues/4384>`_, `#7306
  <https://github.com/scipy/scipy/issues/7306>`_).

  * This would be simply achieved in the new framework via

    * using ``Optimizer.__next__`` to observe different states

    * the user raising ``StopIteration`` in a callback or function eval

  * Proposal: handled in ``Optimizer.solve``, easy to change

  * SciPy currently: each scalar minimizer would have to undergo significant
    changes to implement this, with a try/except around every
    function/callback, and a large amount of duplicate code.

* More comprehensive testing than currently achievable is enabled.

  * Instance methods are common to all classes, and the methods have less
    branching.

  * Deep testing of a single base class method means that all inheriting
    classes are then covered. With the current monolithic minimizer functions
    it is harder to write tests to cover every eventuality.

  * For example with the ``StopIteration`` example given above, the Exception
    could be raised in many places, each of which would have to be tested,
    with slightly different tests for each scalar minimizer.

User interaction
----------------

Introducing a class for optimizations allows intermediate and expert users
to interact more closely with their optimization, a very practical issue.
The main advantage of the ``Optimizer`` class is that it is iterable. This
means that a user is able to control the graining of how the optimization
proceeds. The user can perform step wise iteration; at each step the solution,
as well as solver state, is accessible. The user can therefore:

* access the current value of ``func`` without another function call

* define their own convergence critieria (e.g., accuracy for a
  machine learning model, or total time taken)

* change the optimizer hyper-parameters (e.g. mutation constant
  in differential evolution) during iteration

* encode the gradient (e.g., as in `QSGD`_ and `TernGrad`_)

* allow the user to compute other statistics partway through their
  optimization

  * Especially when the function is expensive to calculate and the
    statistics depend on the current value of the function and derivative

* deal with KeyboardInterrupts without losing state

* halt when they want to

Currently this is not possible, as the optimizers are at the extreme end of
coarse graining. In otherwords they must run through up to `maxiter` steps at
once, optimizer internals are not accessible at each step, and all optimizer
internals are lost when a function returns. The proposed ``Optimizer`` class
can use a coarse grained approach if necessary, by asking ``Optimizer.solve``
to advance by a set number of iterations. Even then the optimizer state is
retained.

.. _QSGD: https://arxiv.org/abs/1610.02132
.. _TernGrad: https://arxiv.org/pdf/1705.07878.pdf

Simplified halting
------------------
The user can halt whenever necessary if they advance in a stepwise fashion.
If they decide to move to completion (using ``Optimizer.solve``) it is still
possible to halt by raising ``StopIteration`` during function evaluation, or in
the callback. Crucially ``Optimizer`` state is retained, and it's still possible
to do further steps. This approach would be very difficult to implement across
the existing function minimizers, such as ``fmin`` (#4384, #7306).

Use as a context manager
------------------------

Each ``Optimizer`` will be a context manager. This enables easy setup of
cleanup actions to be performed by the object when it exits. For example,
some optimizers could make use of multiprocessing, the resources created by
the ``Optimizer`` could be released on __exit__.

``Function`` class is now separate to an Optimizer
--------------------------------------------------

A function and optimizer are fundamentally different things. The function
calculates a value, the other knows how to minimize the function. Their different
nature is ideally represented by different classes.

Currently the optimizers are responsible for calculating gradients and hessians
themselves (if specific functions are not provided) via numerical differentiation.
The default approach for scalar minimizers is to use the same absolute step size
for numerical differentiation for the duration of the optimization. However,
the fd-step size should really be relative to parameter value as it changes.

It is not easy change in the current implementation without placing the onus
on the user to write their own grad function (this is the job of the library),
or without rewriting all the locations in the existing optimizers where numerical
differentiation is required, with concomitant introduction of yet more keywords.
Crucially, if the ``Function`` is separated out into its own class it can then
be responsible for its own differentiation, and the ``Optimizer`` class becomes
agnostic of any differentiation requirements.

The new Function object will offer more options for numerical differentiation
(absolute step, relative step, 2-point/3-point/complex step, bounds). Of course,
the user can still provide their own gradient implementation if preferred.
Use of the ``Function`` class means that ``Optimizer`` will no longer need to have
the `args` parameters. These can be provided directly to the ``Function``
constructor. The constructor will also be able to take a `kwargs` parameter,
which will provide keyword arguments to the user provided callables.

Normally `func`, `grad` and `hess` callables are used to construct a ``Function``
instance. However, the ``Function`` class can be inherited by a user who wishes
to have extra functionality. Because of the object oriented nature of framework,
and the extended ``Function`` will be usable by any ``Optimizer``.

Approximations to class based ``Function`` are already in use in the SciPy
benchmarks in `test_functions.py`_.

.. _test_functions.py: https://github.com/scipy/scipy/blob/895a7741b12c2c3f816bfd27e5249468bea64a26/benchmarks/benchmarks/test_functions.py

This is also the approach being taken in a constrained trust region minimizer in
"ENH: optimize: ``trust-constr`` optimization algorithms [GSoC 2017]" under
`PR#8328`_, in which scalar functions are being described by a class object.

.. _PR#8328: https://github.com/scipy/scipy/pull/8328

Existing interface has room for improvement
-------------------------------------------

Whilst additional method options can be passed to ``minimize``, the details of
those options are harder to access. Moreover, some of the underlying methods
currently lack arguments to change some key behaviour in the way they operate.
We have seen this an issue with:

* expensive functions time-wise
* the ``callback`` argument

Expensive functions time-wise
"""""""""""""""""""""""""""""

If function evaluation is expensive time-wise, there may be some more
optimizations required based on low level function calls. Currently, this
requires rewriting the function and all the functions that call the function
desired to be changed.

A good example is at scikit-learn, where they've rewritten the Newton-CG method
for evaluating expensive functions at `sklearn/utils/optimize.py`_ because they
saw issues with expensive time-wise functions. They provide a modification to
get the function value and gradient with one function call. The proposed class
framework would make this simple -- it could be implemented in the ``Function``
class.

.. _sklearn/utils/optimize.py: https://github.com/scikit-learn/scikit-learn/blob/931fae8753ad0d9cef1c923ba38932074a8d8027/sklearn/utils/optimize.py

``callback`` improvements
"""""""""""""""""""""""""

The ``callback`` argument is not adequate for advanced use. It only sends the
current estimate ``x``, not the values of ``func``, ``grad`` or ``hess``, etc.
The new implementation will send an intermediate ``OptimizeResult`` to the
callback, with each ``Optimizer`` adding method specific attributes (such as
constraints) if required. There will be an additional attribute, 'wall-time',
which will allow the user to halt on elapsed optimization time.
This PR has not considered how a richer callback approach will be exposed to
``minimize``, that should be the focus of a separate PR.

Related issues: `#7425 <https://github.com/scipy/scipy/pull/7425>`_, `#6907
<https://github.com/scipy/scipy/pull/6907>`_, `#4384
<https://github.com/scipy/scipy/pull/4384>`_.



Open issues
-----------

The following open issues/PRs would be significantly easier to be addressed (or
tackled by the user themselves) with subclassing of an Optimizer base class.
That there are many signifies the level of difficulty implementing a coherent
solution across the multiplicity of scipy.optimize minimizer functions.

Issues resolved by easier testing
"""""""""""""""""""""""""""""""""

* `PR#7819 <https://github.com/scipy/scipy/pull/7819>`_ WIP: Basin hopping
  improvements.

  * This PR discusses behavior of how a minimizer should signify
      success/failure, e.g. if a constraint is violated

* `PR#8375 <https://github.com/scipy/scipy/pull/8375>`_: optimize - check that maxiter is not exceeded

  * Correct implementation is inherited by all Optimizers. Testing is simple for all Optimizers

* `#7854 <https://github.com/scipy/scipy/issues/7854>`_: "BUG: L-BFGS-B does two more iterations than specified in maxiter"

  * Again, inheriting would help resolve this, in testing and implementation.

* `#6019 <https://github.com/scipy/scipy/issues/6019>`_: "minimize_scalar doesn't honor disp option".

  * Again, inheriting would help resolve this.

Issues resolved by inheritance
""""""""""""""""""""""""""""""

* `#8419 <https://github.com/scipy/scipy/issues/8419>`_: "some optimize.minimize methods modify the parameter vector
  in-place"

  * This could be inherited by every instance of ``Optimizer``

* `#6673 <https://github.com/scipy/scipy/issues/6673>`_, "return value of scipy.optimize.minimize not consistent for 1D"

  * Again, inheriting would help resolve this, in testing and implementation.

* `6019 <https://github.com/scipy/scipy/issues/6019>`_ minimize_scalar doesn't seem to honor "disp" option

* `5161 <https://github.com/scipy/scipy/issues/5161>`_ "Optimizers reporting success when the minimum is NaN."

  * This would be standardized to make success False

* `4921 <https://github.com/scipy/scipy/issues/4921>`_ "scipy.optimize maxiter option not working as expected"

  * Optimizer.solve standardises for all subclasses

* `3816 <https://github.com/scipy/scipy/issues/3816>`_ wrap_function seems not to be working when wrapper_args is a one element list

  * fix in Optimizer, fix in all subclasses

* `PR#7425 <https://github.com/scipy/scipy/pull/7425>`_ ENH: optimize: more complete callback signature.

  * ``Optimizer.solve`` calls the callback with an intermediate Optimizer result, all Optimizer subclasses inherit.

* `PR#6907 <https://github.com/scipy/scipy/pull/6907>`_ differential_evolution: improve callback

  * ``Optimizer`` base or sub class calls the callback with an intermediate
    Optimizer result

* `PR#4384 <https://github.com/scipy/scipy/pull/4384>`_: ENH: optimize, returning True from callback function halts minimization

  * User code (either in callback, or in objective function) can raise StopIteration as a way of halting minimization.

Issues related to class interaction
"""""""""""""""""""""""""""""""""""

* `#7306 <https://github.com/scipy/scipy/issues/7306>`_ "any way of stopping
  optimization?".

  * Quote: "Beyond a pre-specified iteration limit, I always wanted some way of
    gracefully terminating an optimization routine during execution. I was
    working on problems that took a very long time to solve and sometimes I
    wanted to see what was going on when the algorithm seemed close to a
    solution but never seemed to achieve the termination conditions."

  * User code (either in callback, or in objective function) can raise
    StopIteration as a way of halting minimization.  Alternatively the user can
    move through the optimization stepwise, and make their own decision on when
    to stop.

* `6878 <https://github.com/scipy/scipy/issues/6878>`_ differential_evolution:
  make callback receive fun(xk)

  * Straightforward with stepwise interaction

* `6026 <https://github.com/scipy/scipy/issues/6026>`_ "Replace approx_grad
  with _numdiff.approx_derivative in scipy.optimize"

  * All numerical diff done in ``Function``, fix needed in one place.
    ``Optimizer`` don't need to know.

Documentation issues
""""""""""""""""""""

* `#5832 <https://github.com/scipy/scipy/issues/5832>`_ grad.T should be
  returned but not documented

* `#8031 <https://github.com/scipy/scipy/issues/8031>`_: "Scipy
  optimize.minimize maxfun has confusing behavior".

  * Documentation in one place will make things clear

* `#8373 <https://github.com/scipy/scipy/issues/8373>`_ "scipy.optimize has
  broken my trust."

  * A quote: "This has cost me dozens of hours of debugging time, only to learn
    that something is wrong with L-BFGS-B that causes the output of
    options={'disp': True} to not be the the cost function at a given parameter
    vector."

  * Again, inheriting would help resolve this.

* `5481 <https://github.com/scipy/scipy/issues/5481>`_ "1D root-finding
  interface and documentation could be improved" (asking for a standardised
  approach to root finding).

  * Inheriting Optimizer class for root finding can standardise behaviour?

Existing work
=============

Minimization has been adopted by libraries including SciPy and many related
libraries (e.g., scikit-learn). Optimization has received significant attention
from industry as well -- Google, Facebook, Amazon, Microsoft and Preferred
Networks have developed Keras/Tensorflow, PyTorch, MXNet, CNTK and Chainer
respectively, all of which are for deep learning and heavily involve
optimization. These libraries have Python bindings and are open source. Of
these, PyTorch and Chainer are the most Pythonic libraries.

Chainer and PyTorch both define classes for their optimizers and a similar
interface to the one above. They also provide a function to provide one
optimization step (``step`` in PyTorch, ``update`` in Chainer). Both libraries
have subclassed their base optimizer to provide different optimization
algorithms suited for deep learning (SGD, Adam, RMSprop, etc).

Keras provides a separation of optimizers and functions: any loss function can
be optimized by any of their optimizers. By contrast, scikit-learn takes a
different approach. They hide all optimization from the user, and instead
provide a wealth of different models (which has been a personal frustration),
though they do support the ``fit`` and ``predict``. There's no requirement that
particular algorithms optimize certain functions, though they have an
``SGDRegressor``. Also of note, scikit-learn has rewritten SciPy's optimize.py
at `utils/optimize.py <https://github.com/scikit-learn/scikit-learn/blob/6b5440a9964480ccb0fe1b59ab516d9228186571/sklearn/utils/optimize.py>`_
to specialize it for expensive functions.

scikit-optimize provides an Optimizer class (and is designed for Bayesian
optimization). They provide functions to perform one optimization step (``ask``
and ``tell``, ``run(..., niter=1)``) to perform a single optimization step.

.. note

    Projects related to sklearn: https://github.com/scikit-learn/scikit-learn/blob/4f710cdd088aa8851e8b049e4faafa03767fda10/doc/related_projects.rst

Issues with current codebase
============================
Maintainability/develop/test
----------------------------
- it is very difficult to introduce new features in a uniform manner
  across all scalar minimizers. New functionality has to be
  replicated in lots of places.
- the maintenance of the scipy.optimize codebase is complicated
  by the number of scalar minimizers, changes have to be made in
  lots of places to implement fixes.
- uniform testing of all the scalar minimizers is complicated by their
  number.

.. todo

    FIX ME - add your views on ``minimize`` issues>

Concerns
========

There are already two interfaces - minimize(method='nelder-mead') and fmin, why a third?
----------------------------------------------------------------------------------------
The function based optimizers (fmin, etc) were the original API for scalar
minimizers in scipy. They had slightly different meanings for a lot of their
options (e.g. tol vs ftol/xtol). A more unified way of using those optimizers
were introduced at a later stage with the ``minimize`` function, where the
meanings of each of the keywords was made more uniform, and the selection of
the optimization method was done with a string option. At this stage both
``fmin`` and ``minimize`` call a base ``_minimize_neldermead`` function
to execute the core optimization functionality.

As such ``minimize`` and ``fmin`` do a fair bit of pre-processing to massage
arguments into the form that should be expected by ``_minimize_neldermead``.
There is also checking of arguments, for example ``_minimize_neldermead``
cannot use bounds, so should reject a `bounds` argument.

The ``Optimizer`` objects will provide the core functionality of the
optimization, extracted from ``_minimize_neldermead``. Because backwards
compatibility is required both ``minimize`` and ``fmin`` need to be retained.
Initially the ``Optimizer`` objects will be private, so only the current API's
will be visible. However, full functionality (viz. iteration, etc) of the
framework is only realised when the ``Optimizer`` become public. At this
point there will be three API's. However, we believe the added features of
the classes are compelling. One way of simplification could be to deprecate
``fmin``-like functions for subsequent removal - a decision to be made at a
future point. The constructor for each of the ``Optimizers`` will be consistent
with parameter naming from ``minimize``.

Why not just enhance the callback with the existing minimizers?
---------------------------------------------------------------
Ironically the easiest way to enhance the callback for all optimizers is to
have a class based system where each Optimizer inherits a standardised way of
providing a user callback with more detailed information (write once, all
benefit). To  implement this for existing optimizers without a class based
framework will involve changes to many locations in the codebase, and is harder
to uniformly implement and test. In the new framework the callback will receive
an intermediate ``OptimizeResult`` with extra attributes (such as solver state)
being added if needed.

Will allowing the user to access Optimizer attributes lead to problems?
-----------------------------------------------------------------------
Currently there is no way of accessing and changing optimization hyper-parameters
as optimization is proceeding. As outlined above, the ability to change those
attributes are of interest to the advanced user. However, the ability to inspect
instance attributes may tempt the user to change attributes that would break
the optimizers operation. Class attributes that shouldn't be modified will be
given names conforming to usual Python practice (i.e. prefixed by _ or __).
Solver hyperparameters will be made available in the dict returned by the
hyper_parameter attribute.

Won't the code be more complex to maintain, because it inverts structure of control flow?
-----------------------------------------------------------------------------------------

The class based nature means that if the base class is improved, then
all inheriting classes are improved without extra work, and the tests
automatically cover all the inheriting classes. Consider, for example,
the ability for a solver to halt on demand (#4384, #7306). This was
implemented by four lines of code in the ``Optimizer.solve`` method. The
user simply raises a StopIteration in their objective function, or in
their callback. Moreover the solver state is not lost, it's still
possible to obtain the current solution, and optimization can continue
if required.

To implement the same in each scalar minimizer would be very difficult.
Each function and callback location would have to be wrapped with a
try/except clause (multiple locations in each optimizer), it gets harder
to return the current solution, and all optimization state is lost.
Implementing the class based structure is relatively fast, and there
is pre-existing testing from the current test suite to check correctness.
The main task is to override the __next__ method with the core iteration
process.

It's worth noting that the ``DifferentialEvolutionSolver`` and
``BasinHoppingRunner`` are already written in class form. Indeed,
``DifferentialEvolutionSolver`` has naturally moved (evolved) towards
this proposed implementation since it was first implemented (so we
have experience that the gross design works).

Will existing foreign functionality in C/Fortran be moved into Python?
-----------------------------------------------------------------------
Calculation by foreign functions can be split into two:

 - Iteration occurs in Python, but the foreign function suggests a new
   step (LBFGSB)
 - All looping takes place in C/Fortran, with that code calling back into
   Python for objective evaluation (e.g. ``leastsq``)

The former can simply be translated to the class based framework, with
the external calls still taking place (already completed for LBFGSB).
The second cannot (without significant rewrite). It's not intended at the
moment to rewrite the latter. It's not mandatory that all optimizers are
converted to use the class based code, just those that are possible.

Do people actually care about class based Optimizers?
-----------------------------------------------------

Depends on the user. For basic use, the user won't care: they only want a
functional interface. However, for anything more advanced it's sorely lacking.
A class based API can provide enhanced interaction, which will be useful to the
intermediate and advanced user.  We still believe the class based approach is
easier to maintain, even if it's advanced features aren't used.

Will API changes be made?
-------------------------
The ``minimizer`` and individual scalar minimizer functions (such as
``fmin``) can continue to be used as before, see back-compatibility section
below. Once made public the  new Optimizer objects can be used by themselves
to access richer functionality.

Won't performance suffer?
-------------------------
The benchmarking of the new approach is outlined in the Implementation
section below.

Why not just implement a scikit?
--------------------------------
The major reason why many people/projects use scipy is the optimize
module, it's a cornerstone of the project. It's worth implementing
``Optimizer`` on the maintainability argument alone, in which case
this question does not apply. If one is asking that question with
specific regard to the introduction of more advanced features, then
the question that would need asking in that case is:
"why have scipy at all, why not just have multiple scikits?"

Why not have a similar interface for root finding and linear programs?
----------------------------------------------------------------------

We have personally experienced issues with the ``minimize`` interface and API.
We suggest changing it because of this. We have not experienced similar issues
with root finding and linear programs. If other have, please come forward.

.. note

    * We have personal experience that makes minimize a problem. We are open to
      expanding this class interface but currently see no need to expand
      root/minimize_scalar/linprog.
    * `minimize` is similar to `solve_ivp` (see
      https://github.com/scipy/scipy/pull/8414#issuecomment-366372052) I said
      "minimize has been an issue to me". Can point to other examples.  and
      implementing classes could lower barrier to implementing new minimizers

Implementation
==============
An Optimizer and Function class will be created. Using two classes clearly
separates their functionality, for example, it shouldn't be necessary for a
minimizer to worry about how gradients are calculated.

Scope
-----
The scalar minimizers will be tackled first, with NelderMead, LBFGB, BFGS,
DifferentialEvolutionSolver in the initial PR.  Subsequent PR's will modify
BasinHopping, with the rest following. It should be noted here that
the modifications to BasinHopping and DifferentialEvolution are relatively
minor, given that they substantially already conform to the new class based
form.

Staggered modifications permit simpler code review, and allow fine tuning of
the Optimizer and Function classes as work continues.

Extension to root finding, least-squares, and minimize scalar, can be
performed at a later date. The ``Optimizer`` class will be written in
such a way that those optimizers are tacklable with this class based
framework.

Backwards compatibility
-----------------------
Both the ``minimize``, and ``fmin``, etc, functions will continue to work
with the current form. However, at their core calculation will be carried
out by the various ``Optimizer`` objects. Once the Optimizer classes are
exposed to the scipy public API the new objects can be used by themselves.
Perfect back compatibility will not be kept if the class implementation
discovers bugs/flaws in existing code. One example of this are the current
minimizers that perform more than `maxiter` iterations.

Testing
-------
The Optimizers will be exposed to the the existing test suite for the scalar
minimizers as they will constitute the core functionality for the existing
functions.

New tests will be written for each Optimizer and Function class. The tests
are parameterised across the Optimizer classes meaning that consistent
behaviour is achieved.

Performance
-----------
Initial benchmarking of ``Optimizer`` for the ``optimize`` asv benchmark suite
are listed below (Python 3.6, current numpy). 3cb29828 is current master, and
1bbf5b1a is from the proposed PR. Only the 'scalar minimizer + problem'
combinations that showed significant speedup/slowdown are tabulated. There is
a general slowdown for those methods (``nelder-mead``, ``bfgs``, ``l-bfgs-b``)
using an ``Optimizer`` in 1bbf5b1a. However, this slowdown must be considered
against a wider context:

- maintenance/develop/test for the class based optimizers is much simpler (c.f.
  Knuth's "Premature Optimization" quote)
- the evaluation time for objective functions posed by users are going to range
  over several orders of magnitude. Some will be sub microsecond, others will be
  multi-minute. Only those at the fastest end are affected by optimizer overhead.
- if a slowdown on the order of 2 ms to 3 ms is significant, then one may argue
  that other approaches are needed (such as compiled extensions).

There is some degree of uncertainty in the actual timings. For example,
'rosenbrock_tight + CG', or 'LJ + COBYLA' have a slowdown factor of ~1.23, yet
neither of those solvers use an ``Optimizer``. This may indicate that more
repeats are necessary to obtain more accurate timings.

+---+------------+------------+-------+------------------------+----------------+-----------+
|   | before     | after      | ratio |  problem               | method         | resource  |
|   | [3cb29828] | [1bbf5b1a] |       |                        |                |           |
+===+============+============+=======+========================+================+===========+
| + | 544.76μs   | 1.02ms     | 1.87  | sin1d                  | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 178.85μs   | 283.61μs   | 1.59  | booth                  | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 163.45μs   | 256.32μs   | 1.57  | sin1d                  | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 132.17μs   | 201.15μs   | 1.52  | asymmetric_quadratic   | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 12.88ms    | 19.27ms    | 1.50  | rosenbrock_slow        | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 133.09μs   | 197.68μs   | 1.49  | simple_quadratic       | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 158.79μs   | 234.93μs   | 1.48  | simple_quadratic       | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 1.75ms     | 2.57ms     | 1.46  | rosenbrock             | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 224.64μs   | 325.62μs   | 1.45  | sin1d                  | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 2.12ms     | 3.03ms     | 1.43  | simple_quadratic       | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 6.71ms     | 9.41ms     | 1.40  | rosenbrock_nograd      | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 2.06ms     | 2.85ms     | 1.39  | rosenbrock_nograd      | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 2.30ms     | 3.18ms     | 1.38  | beale                  | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 3.70ms     | 5.08ms     | 1.37  | rosenbrock             | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 36.47ms    | 49.28ms    | 1.35  | rosenbrock_slow        | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 2.32ms     | 3.09ms     | 1.33  | asymmetric_quadratic   | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 1.82ms     | 2.37ms     | 1.30  | booth                  | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 1.82ms     | 2.30ms     | 1.27  | rosenbrock_tight       | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 8.94ms     | 11.18ms    | 1.25  | rosenbrock_tight       | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 5.59ms     | 6.86ms     | 1.23  | rosenbrock_tight       | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 800.57μs   | 979.82μs   | 1.22  | beale                  | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 28.85ms    | 35.24ms    | 1.22  | rosenbrock_slow        | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 32.87ms    | 40.01ms    | 1.22  | LJ                     | COBYLA         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 4.55ms     | 5.52ms     | 1.21  | rosenbrock_nograd      | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 7.36ms     | 8.60ms     | 1.17  | rosenbrock             | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 209.12μs   | 240.41μs   | 1.15  | asymmetric_quadratic   | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 4.20ms     | 4.82ms     | 1.15  | rosenbrock             | trust-exact    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 6.19ms     | 7.06ms     | 1.14  | rosenbrock_tight       | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 4.99ms     | 5.50ms     | 1.10  | rosenbrock             | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 124.39ms   | 134.51ms   | 1.08  | LJ                     | nelder-mead    | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 89.75ms    | 96.86ms    | 1.08  | LJ                     | Powell         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 121.43μs   | 130.99μs   | 1.08  | simple_quadratic       | COBYLA         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 185.56μs   | 200.15μs   | 1.08  | asymmetric_quadratic   | Newton-CG      | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 101.36ms   | 107.99ms   | 1.07  | rosenbrock_slow        | COBYLA         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 13.63ms    | 14.52ms    | 1.06  | LJ                     | L-BFGS-B       | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 1.40ms     | 1.49ms     | 1.06  | asymmetric_quadratic   | Powell         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 262.67μs   | 277.57μs   | 1.06  | simple_quadratic       | dogleg         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| + | 935.46μs   | 984.38μs   | 1.05  | beale                  | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 20.26ms    | 19.20ms    | 0.95  | beale                  | Powell         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 9.00s      | 8.50s      | 0.94  | simple_quadratic       | BFGS           | mean_nfev |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 223.78μs   | 209.74μs   | 0.94  | booth                  | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 120.76μs   | 113.15μs   | 0.94  | simple_quadratic       | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 11.10ms    | 10.36ms    | 0.93  | rosenbrock_tight       | COBYLA         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 97.75μs    | 91.05μs    | 0.93  | booth                  | TNC            | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 8.15ms     | 7.59ms     | 0.93  | rosenbrock_tight       | Newton-CG      | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 1.60ms     | 1.47ms     | 0.92  | beale                  | COBYLA         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 101.40μs   | 92.60μs    | 0.91  | asymmetric_quadratic   | TNC            | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 173.31μs   | 155.38μs   | 0.90  | sin1d                  | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 821.45μs   | 720.21μs   | 0.88  | booth                  | Powell         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 25.50ms    | 21.98ms    | 0.86  | LJ                     | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 44.70ms    | 37.89ms    | 0.85  | LJ                     | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 7.67ms     | 6.49ms     | 0.85  | LJ                     | SLSQP          | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 17.89ms    | 14.72ms    | 0.82  | LJ                     | TNC            | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 299.45μs   | 246.45μs   | 0.82  | booth                  | SLSQP          | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 103.64μs   | 83.85μs    | 0.81  | sin1d                  | COBYLA         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 434.90μs   | 340.53μs   | 0.78  | sin1d                  | Powell         | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 14.16ms    | 10.42ms    | 0.74  | beale                  | BFGS           | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 145.10μs   | 94.44μs    | 0.65  | sin1d                  | TNC            | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 183.18μs   | 113.39μs   | 0.62  | asymmetric_quadratic   | CG             | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+
| - | 313.59μs   | 193.45μs   | 0.62  | sin1d                  | SLSQP          | mean_time |
+---+------------+------------+-------+------------------------+----------------+-----------+

``Optimizer``: common methods and attributes
--------------------------------------------

Methods
"""""""

* ``__init__(self, Function, <minimizer specific options>)``: initialization of
  an individual optimizer
* ``__next__(self)``:
			the core of the optimization functionality, performs a single
			optimization iteration.
			returns solution and function: x, f, at each iteration.
* ``__call__(self, iterations, callback=None, maxfun=np.inf, allvecs=False)``
			performs a set number of iterations, with a max number of function evaluations
			calls user callback at end of each iteration with intermediate OptimizeResult
* ``solve(self, maxiter=np.inf, callback=None, maxfun=np.inf, allvecs=False)``:
			performs a set number of iterations, with a max number of function evaluations
* ``_finish_up(self)``:  defines cleanup actions when solve finishes.
* ``__enter__(self)``:
			performs setup if Optimizer is used as context manager
* ``__exit__(self)``:
			performs cleanup if Optimizer is used as context manager
			e.g. closing multiprocessing pools.

Attributes
""""""""""

* ``converged``: truth as to whether Optimizer has converged
* ``N``: dimensionality of problem
* ``x0``: initial guess
* ``x``: current best solution
* ``fun``: current function evaluation at current solution
* ``hyper_parameters``: dict containing optimizer hyper parameters
* ``result``: current OptimizerResult (returned by solver)
* ``nit``: total number of iterations peformed by Optimizer
* ``nfev``: number of function evaluations used by Optimizer
* ``njev``: number of gradient evaluations used by Optimizer
* ``nhev``: number of Hessian evaluations used by Optimizer


``Function``: methods and attributes
-------------------------------------

The Function class is responsible for evaluating its function, its gradient, and
its Hessian. Minimization of scalar functions and vector functions will require
separate implementations, but will have the same methods.

* ``__init__(self, func=None, grad=None, hess=None, fd_method='3-point', step=None args=(), kwargs=None)``:
			args and kwargs are passed to func, grad, and hess.
* ``func(self)``: implementation of the function
* ``grad(self, f0=None)``: implementation of the gradient
* ``hess(self)``: implementation of the Hessian

There will be different ways of creating a function. Either the Function can be
initialised with `func`, `grad`, `hess` callables, or a Function may be
subclassed. If the Function is not subclassed then it must be initialised with
a `func` callable. If `grad` and `hess` are not provided, or not overridden,
then the gradient and hessian will be numerically estimated with finite
differences. The finite differences will either be absolute or relative step
(approx_fprime or approx_derivative), and controlled by the `fd_method` or
`step` keywords.


Example usage
-------------

.. code-block:: python

    import time
    from scipy.optimize import Function, BFGS, minimize

    def func(x, *args):
        return x**2 + args[0]
    def grad(x, *args):
        return 2 * x

    def callback(res): print(res)

    x0 = [2.0]

    # existing call has lots of parameters, mixing optimizer args with func args
    # it might be nice to have **kwds as well, but not possible with current approach
    result = minimize(func, x0, args=(2,), jac=grad, method='BFGS', maxiter=10, callback=callback)

    # proposed
    function = Function(func=func, args=(2,), grad=grad)
    opt = BFGS(function, x0)

    # do 10 iterations
    result = opt.solve(maxiter=10, callback=callback)

    # alternatively control how iteration occurs/halts
    d = opt.hyper_parameters
	start_time = time.time()
    for i, v in enumerate(opt):
        x, f = v
        print(i, f, x)
        d['my_hyper_parameter'] = np.inf
	    if time.time() - start_time > 100:
	  	    break

    # use of Function Class offers the potential for more sophisticated calculation.

    class Quad(Function):
        def __init__(self, bkg):
            super(Quad, self).__init__(self)
            self.bkg = bkg

        def func(self, x):
            return (x**2 + args[0])

        def grad(self, x):
            return 2*x

        def hess(self, x):
            return 2

    opt = BFGS(Quad, x0).solve(maxiter=10)

This is an example of machine learning. A function (``L2Loss``) is defined and
needs to be minimized over different training examples.

.. code-block:: python

    from scipy.optimize import Function, Optimizer

	# create your own GradientDescent optimizer, minimizing L2norm
    class L2Loss(Function):
        def __init__(self, A, y, *args, **kwargs):
            self.A = A
            self.y = y
            super().__init__(self, *args, **kwargs)

        def func(x):
            return LA.norm(self.A@x - self.y)**2

        def grad(x):
            return 2 * self.A.T @ (self.A@x - self.y)

    class GradientDescent(Optimizer):
        def __init__(self, *args, step_size=1e-3, **kwargs):
            self.step_size = step_size
            super().__init__(*arg, **kwargs)

        def __next__(self):
            self.x -= self.step_size*self.grad(x)

	n, d = 100, 10
	A = np.random.randn(n, d)
	x_star = np.random.randn(d)
	y = np.sign(A @ x_star)

	loss = L2Loss(A, y)
	opt = GradientDescent(loss)

	for k, _ in enumerate(opt):  # Optimizer.__next__ implement minimization
		if k % 100 == 0:
			compute_stats(opt, loss)


Context managers offer the chance for cleanup actions, for example multiprocessing:

.. code-block:: python
   
   with DifferentialEvolutionSolver(function, bounds, workers=2) as opt:
        # the __entry__ and __exit__ in the solver could create and close
        # multiprocessing pools.
        res = opt.solve()

Other possible scipy.optimize improvements
------------------------------------------

These suggestions aren't directly related to the PEP, but are worth
mentioning to provide a picture of other possible ``optimize``
improvements.

1. Cythonise scalar minimizers, and allow ``LowLevelCallable`` as the
   objective function. ``approx_derivative`` and ``approx_fprime`` (and
   other support functions) would have to be ported. The real speedup
   would not be Cythonisation per-se, but the use of `LowLevelCallable`
   from Cython, avoiding Python call overhead. Experience gained during
   this project might aid application in ``leastsq``, if
   ``LowLevelCallables`` could be integrated with its Fortran nature.
2. Implement (optional) parallelisation in calculation of finite
   difference gradients/Hessians (e.g. ``approx_derivative``).
3. Use parallelisation in `optimize.brute`
4. Use parallelisation in NelderMead (10.1007/s10614-007-9094-2, http://www.econ.nyu.edu/user/wiswall/research/lee_wiswall_parallel_simplex_edit_2_8_2007.pdf).

Copyright
=========
This document has been placed in the public domain, Feb 2018.
