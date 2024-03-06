Random variables
----------------

There are two general distribution classes that have been implemented
for encapsulating :ref:`continuous random variables
<continuous-random-variables>` and :ref:`discrete random variables
<discrete-random-variables>`. Over 80 continuous random variables
(RVs) and 10 discrete random variables have been implemented using
these classes. Besides this, new routines and distributions can be
easily added by the end user. (If you create one, please contribute it.)

All of the statistics functions are located in the sub-package
:mod:`scipy.stats` and a fairly complete listing of these functions
can be obtained using ``info(stats)``. The list of the random
variables available can also be obtained from the docstring for the
stats sub-package.

In the discussion below, we mostly focus on continuous RVs. Nearly everything
also applies to discrete variables, but we point out some differences
here: :ref:`discrete_points_label`.

In the code samples below, we assume that the :mod:`scipy.stats` package
is imported as

    >>> from scipy import stats

and in some cases we assume that individual objects are imported as

    >>> from scipy.stats import norm

Getting help
^^^^^^^^^^^^

First of all, all distributions are accompanied with help
functions. To obtain just some basic information, we print the relevant
docstring: ``print(stats.norm.__doc__)``.

To find the support, i.e., upper and lower bounds of the distribution,
call:

    >>> print('bounds of distribution lower: %s, upper: %s' % norm.support())
    bounds of distribution lower: -inf, upper: inf

We can list all methods and properties of the distribution with
``dir(norm)``. As it turns out, some of the methods are private,
although they are not named as such (their names do not start
with a leading underscore), for example ``veccdf``, are only available
for internal calculation (those methods will give warnings when one tries to
use them, and will be removed at some point).

To obtain the *real* main methods, we list the methods of the frozen
distribution. (We explain the meaning of a `frozen` distribution
below).

    >>> rv = norm()
    >>> dir(rv)  # reformatted
    ['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__',
     '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__',
     '__init__', '__le__', '__lt__', '__module__', '__ne__', '__new__',
     '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
     '__str__', '__subclasshook__', '__weakref__', 'a', 'args', 'b', 'cdf',
     'dist', 'entropy', 'expect', 'interval', 'isf', 'kwds', 'logcdf',
     'logpdf', 'logpmf', 'logsf', 'mean', 'median', 'moment', 'pdf', 'pmf',
     'ppf', 'random_state', 'rvs', 'sf', 'stats', 'std', 'var']

Finally, we can obtain the list of available distribution through
introspection:

    >>> dist_continu = [d for d in dir(stats) if
    ...                 isinstance(getattr(stats, d), stats.rv_continuous)]
    >>> dist_discrete = [d for d in dir(stats) if
    ...                  isinstance(getattr(stats, d), stats.rv_discrete)]
    >>> print('number of continuous distributions: %d' % len(dist_continu))
    number of continuous distributions: 108
    >>> print('number of discrete distributions:   %d' % len(dist_discrete))
    number of discrete distributions:   20

Common methods
^^^^^^^^^^^^^^

The main public methods for continuous  RVs are:

* rvs:   Random Variates
* pdf:   Probability Density Function
* cdf:   Cumulative Distribution Function
* sf:    Survival Function (1-CDF)
* ppf:   Percent Point Function (Inverse of CDF)
* isf:   Inverse Survival Function (Inverse of SF)
* stats: Return mean, variance, (Fisher's) skew, or (Fisher's) kurtosis
* moment: non-central moments of the distribution


Let's take a normal RV as an example.

    >>> norm.cdf(0)
    0.5

To compute the ``cdf`` at a number of points, we can pass a list or a numpy array.

    >>> norm.cdf([-1., 0, 1])
    array([ 0.15865525,  0.5,  0.84134475])
    >>> import numpy as np
    >>> norm.cdf(np.array([-1., 0, 1]))
    array([ 0.15865525,  0.5,  0.84134475])

Thus, the basic methods, such as `pdf`, `cdf`, and so on, are vectorized.

Other generally useful methods are supported too:

    >>> norm.mean(), norm.std(), norm.var()
    (0.0, 1.0, 1.0)
    >>> norm.stats(moments="mv")
    (array(0.0), array(1.0))

To find the median of a distribution, we can use the percent point
function ``ppf``, which is the inverse of the ``cdf``:

    >>> norm.ppf(0.5)
    0.0

To generate a sequence of random variates, use the ``size`` keyword
argument:

    >>> norm.rvs(size=3)
    array([-0.35687759,  1.34347647, -0.11710531])   # random

Don't think that ``norm.rvs(5)`` generates 5 variates:

    >>> norm.rvs(5)
    5.471435163732493  # random

Here, ``5`` with no keyword is being interpreted as the first possible
keyword argument, ``loc``, which is the first of a pair of keyword arguments
taken by all continuous distributions.
This brings us to the topic of the next subsection.

Random number generation
^^^^^^^^^^^^^^^^^^^^^^^^

Drawing random numbers relies on generators from `numpy.random` package.
In the examples above, the specific stream of
random numbers is not reproducible across runs. To achieve reproducibility,
you can explicitly *seed* a random number generator. In NumPy, a generator
is an instance of `numpy.random.Generator`. Here is the canonical way to create
a generator:

    >>> from numpy.random import default_rng
    >>> rng = default_rng()

And fixing the seed can be done like this:

    >>> # do NOT copy this value
    >>> rng = default_rng(301439351238479871608357552876690613766)

.. warning:: Do not use this number or common values such as 0. Using just a
             small set of seeds to instantiate larger state spaces means that
             there are some initial states that are impossible to reach. This
             creates some biases if everyone uses such values. A good way to
             get a seed is to use a `numpy.random.SeedSequence`:

             >>> from numpy.random import SeedSequence
             >>> print(SeedSequence().entropy)
             301439351238479871608357552876690613766  # random

The `random_state` parameter in distributions accepts an instance of
`numpy.random.Generator` class, or an integer, which is then used to
seed an internal ``Generator`` object:

    >>> norm.rvs(size=5, random_state=rng)
    array([ 0.47143516, -1.19097569,  1.43270697, -0.3126519 , -0.72058873])  # random

For further info, see `NumPy's documentation
<https://numpy.org/doc/stable/reference/random/index.html>`__.

To learn more about the random number samplers implemented in SciPy, see
:ref:`non-uniform random number sampling tutorial
<non-uniform-random-number-sampling>` and :ref:`quasi monte carlo tutorial
<quasi-monte-carlo>`

Shifting and scaling
^^^^^^^^^^^^^^^^^^^^

All continuous distributions take ``loc`` and ``scale`` as keyword
parameters to adjust the location and scale of the distribution,
e.g., for the standard normal distribution, the location is the mean and
the scale is the standard deviation.

    >>> norm.stats(loc=3, scale=4, moments="mv")
    (array(3.0), array(16.0))

In many cases, the standardized distribution for a random variable ``X``
is obtained through the transformation ``(X - loc) / scale``. The
default values are ``loc = 0`` and ``scale = 1``.

Smart use of ``loc`` and ``scale`` can help modify the standard
distributions in many ways. To illustrate the scaling further, the
``cdf`` of an exponentially distributed RV with mean :math:`1/\lambda`
is given by

.. math::

    F(x) = 1 - \exp(-\lambda x)

By applying the scaling rule above, it can be seen that by
taking ``scale  = 1./lambda`` we get the proper scale.

    >>> from scipy.stats import expon
    >>> expon.mean(scale=3.)
    3.0

.. note:: Distributions that take shape parameters may
   require more than simple application of ``loc`` and/or
   ``scale`` to achieve the desired form. For example, the
   distribution of 2-D vector lengths given a constant vector
   of length :math:`R` perturbed by independent N(0, :math:`\sigma^2`)
   deviations in each component is
   rice(:math:`R/\sigma`, scale= :math:`\sigma`). The first argument
   is a shape parameter that needs to be scaled along with :math:`x`.

The uniform distribution is also interesting:

    >>> from scipy.stats import uniform
    >>> uniform.cdf([0, 1, 2, 3, 4, 5], loc=1, scale=4)
    array([ 0.  ,  0.  ,  0.25,  0.5 ,  0.75,  1.  ])


Finally, recall from the previous paragraph that we are left with the
problem of the meaning of ``norm.rvs(5)``. As it turns out, calling a
distribution like this, the first argument, i.e., the 5, gets passed
to set the ``loc`` parameter. Let's see:

    >>> np.mean(norm.rvs(5, size=500))
    5.0098355106969992  # random

Thus, to explain the output of the example of the last section:
``norm.rvs(5)`` generates a single normally distributed random variate with
mean ``loc=5``, because of the default ``size=1``.

We recommend that you set ``loc`` and ``scale`` parameters explicitly, by
passing the values as keywords rather than as arguments. Repetition
can be minimized when calling more than one method of a given RV by
using the technique of `Freezing a Distribution`_, as explained below.


Shape parameters
^^^^^^^^^^^^^^^^

While a general continuous random variable can be shifted and scaled
with the ``loc`` and ``scale`` parameters, some distributions require
additional shape parameters. For instance, the gamma distribution with density

.. math::

    \gamma(x, a) = \frac{\lambda (\lambda x)^{a-1}}{\Gamma(a)} e^{-\lambda x}\;,

requires the shape parameter :math:`a`. Observe that setting
:math:`\lambda` can be obtained by setting the ``scale`` keyword to
:math:`1/\lambda`.

Let's check the number and name of the shape parameters of the gamma
distribution. (We know from the above that this should be 1.)

    >>> from scipy.stats import gamma
    >>> gamma.numargs
    1
    >>> gamma.shapes
    'a'

Now, we set the value of the shape variable to 1 to obtain the
exponential distribution, so that we compare easily whether we get the
results we expect.

    >>> gamma(1, scale=2.).stats(moments="mv")
    (array(2.0), array(4.0))

Notice that we can also specify shape parameters as keywords:

   >>> gamma(a=1, scale=2.).stats(moments="mv")
   (array(2.0), array(4.0))


Freezing a distribution
^^^^^^^^^^^^^^^^^^^^^^^

Passing the ``loc`` and ``scale`` keywords time and again can become
quite bothersome. The concept of `freezing` a RV is used to
solve such problems.

    >>> rv = gamma(1, scale=2.)

By using ``rv`` we no longer have to include the scale or the shape
parameters anymore. Thus, distributions can be used in one of two
ways, either by passing all distribution parameters to each method
call (such as we did earlier) or by freezing the parameters for the
instance of the distribution. Let us check this:

    >>> rv.mean(), rv.std()
    (2.0, 2.0)

This is, indeed, what we should get.


Broadcasting
^^^^^^^^^^^^

The basic methods ``pdf``, and so on, satisfy the usual numpy broadcasting rules. For
example, we can calculate the critical values for the upper tail of
the t distribution for different probabilities and degrees of freedom.

    >>> stats.t.isf([0.1, 0.05, 0.01], [[10], [11]])
    array([[ 1.37218364,  1.81246112,  2.76376946],
           [ 1.36343032,  1.79588482,  2.71807918]])

Here, the first row contains the critical values for 10 degrees of freedom
and the second row for 11 degrees of freedom (d.o.f.). Thus, the
broadcasting rules give the same result of calling ``isf`` twice:

    >>> stats.t.isf([0.1, 0.05, 0.01], 10)
    array([ 1.37218364,  1.81246112,  2.76376946])
    >>> stats.t.isf([0.1, 0.05, 0.01], 11)
    array([ 1.36343032,  1.79588482,  2.71807918])

If the array with probabilities, i.e., ``[0.1, 0.05, 0.01]`` and the
array of degrees of freedom i.e., ``[10, 11, 12]``, have the same
array shape, then element-wise matching is used. As an example, we can
obtain the 10% tail for 10 d.o.f., the 5% tail for 11 d.o.f. and the
1% tail for 12 d.o.f. by calling

    >>> stats.t.isf([0.1, 0.05, 0.01], [10, 11, 12])
    array([ 1.37218364,  1.79588482,  2.68099799])


.. _discrete_points_label:

Specific points for discrete distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Discrete distributions have mostly the same basic methods as the
continuous distributions. However ``pdf`` is replaced by the probability
mass function ``pmf``, no estimation methods, such as fit, are
available, and ``scale`` is not a valid keyword parameter. The
location parameter, keyword ``loc``, can still be used to shift the
distribution.

The computation of the cdf requires some extra attention. In the case
of continuous distribution, the cumulative distribution function is, in
most standard cases, strictly monotonic increasing in the bounds (a,b)
and has, therefore, a unique inverse. The cdf of a discrete
distribution, however, is a step function, hence the inverse cdf,
i.e., the percent point function, requires a different definition:

::

    ppf(q) = min{x : cdf(x) >= q, x integer}

For further info, see the docs :ref:`here<discrete-ppf>`.


We can look at the hypergeometric distribution as an example

    >>> from scipy.stats import hypergeom
    >>> [M, n, N] = [20, 7, 12]

If we use the cdf at some integer points and then evaluate the ppf at those
cdf values, we get the initial integers back, for example

    >>> x = np.arange(4) * 2
    >>> x
    array([0, 2, 4, 6])
    >>> prb = hypergeom.cdf(x, M, n, N)
    >>> prb
    array([  1.03199174e-04,   5.21155831e-02,   6.08359133e-01,
             9.89783282e-01])
    >>> hypergeom.ppf(prb, M, n, N)
    array([ 0.,  2.,  4.,  6.])

If we use values that are not at the kinks of the cdf step function, we get
the next higher integer back:

    >>> hypergeom.ppf(prb + 1e-8, M, n, N)
    array([ 1.,  3.,  5.,  7.])
    >>> hypergeom.ppf(prb - 1e-8, M, n, N)
    array([ 0.,  2.,  4.,  6.])


Fitting distributions
^^^^^^^^^^^^^^^^^^^^^

The main additional methods of the not frozen distribution are related
to the estimation of distribution parameters:

* fit:   maximum likelihood estimation of distribution parameters, including location
         and scale
* fit_loc_scale: estimation of location and scale when shape parameters are given
* nnlf:  negative log likelihood function
* expect: calculate the expectation of a function against the pdf or pmf


.. _performance_issues_label:

Performance issues and cautionary remarks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The performance of the individual methods, in terms of speed, varies
widely by distribution and method. The results of a method are
obtained in one of two ways: either by explicit calculation, or by a
generic algorithm that is independent of the specific distribution.

Explicit calculation, on the one hand, requires that the method is
directly specified for the given distribution, either through analytic
formulas or through special functions in ``scipy.special`` or
``numpy.random`` for ``rvs``. These are usually relatively fast
calculations.

The generic methods, on the other hand, are used if the distribution
does not specify any explicit calculation. To define a distribution,
only one of pdf or cdf is necessary; all other methods can be derived
using numeric integration and root finding. However, these indirect
methods can be `very` slow. As an example, ``rgh =
stats.gausshyper.rvs(0.5, 2, 2, 2, size=100)`` creates random
variables in a very indirect way and takes about 19 seconds for 100
random variables on my computer, while one million random variables
from the standard normal or from the t distribution take just above
one second.


Remaining issues
^^^^^^^^^^^^^^^^

The distributions in ``scipy.stats`` have recently been corrected and improved
and gained a considerable test suite; however, a few issues remain:

* The distributions have been tested over some range of parameters;
  however, in some corner ranges, a few incorrect results may remain.
* The maximum likelihood estimation in `fit` does not work with
  default starting parameters for all distributions and the user
  needs to supply good starting parameters. Also, for some
  distribution using a maximum likelihood estimator might
  inherently not be the best choice.
  
