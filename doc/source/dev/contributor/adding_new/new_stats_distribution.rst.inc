Adding A New Statistics Distribution
------------------------------------

For hundreds of years statisticians, mathematicians and scientists have needed
to understand, analyze and model data.
This has led to a plethora of statisics distributions, many of which are
related to others.
Modeling of new types of data continues to give rise to new distributions,
as does theoretical considerations being applied to new disciplines.
SciPy models about a dozen discrete distributions
:ref:`discrete-random-variables` and 100 continuous distributions
:ref:`continuous-random-variables`.

To add a new distribution, a good reference is needed.
Scipy typically uses [JKB]_ as its gold standard, with WikipediaDistributions_
articles often providing some extra details and/or graphical plots.

How to create a new continuous distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a few steps to be done to add continuous distribution XXXXX to SciPy.
(Adding a discrete distribution is similar).


Before Implementation
.....................

1. See if it already been implemented (that saves a lot of effort!)

  - It may have been implemented with a different name.
  - It may have been implemented with a different parameterization (shape parameters).
  - It may be a specialization of a more general family of distributions.

It is very common for multiple disciplines to discover/rediscover a
distribution (or a specialization or different parameterization).
There are a few existing SciPy distributions which are specializations
of other distributions.  E.g. The :class:`scipy.stats.arcsine` distribution
is a specialization of the :class:`scipy.stats.beta` distribution.
These duplications exist for (very!) historical and widespread usage reasons.
At this time, adding new specializations/reparametrizations of existing
distributions to SciPy is not supported,  mainly due to the increase in user
confusion resulting from such additions.

2. Create a SciPy Issue on github listing the distribution, references and
   reasons for its inclusion.


Implementation
..............

#. Find an already existing distribution similar to ``XXXXX``.
   Use its code as a template for ``XXXXX``.
#. Read the docstring for class ``rv_continuous`` in
   ``scipy/stats/_distn_infrastructure.py``.
#. Write the new code for class ``xxxx_gen`` and insert it into
   ``scipy/stats/_continuous_distns.py``,  which is in (mostly) alphabetical
   order by distribution name.
#. Does the distribution have infinite support? If not, left and/or right
   endpoints ``a``, ``b`` need to be specified in the call to ``xxxx_gen()``.
#. If the support depends upon the shape parameters,
   ``xxxx_gen._get_support()`` needs to be implemented.
#. The default ``_argcheck()`` implementation checks that the shape parameters
   are positive.  Create a more appropriate implementation.
#. If ``xxxx_gen.ppf()`` is expensive to compute relative to
   ``xxxx_gen.pdf()``, consider setting the ``momtype`` in the call
   to ``xxxx_gen()``.
#. If ``xxxx_gen.rvs()`` is expensive to compute, consider implementing a
   specific ``xxxx_gen._rvs()``.
#. Add the name to the listing in the docstring of
   ``scipy/stats/__init__.py``.
#. Add the name and a good set of example shape parameters to the ``distcont``
   list in ``scipy/stats/_distr_params.py``. These shape parameters are used
   for testing and automatic documentation generation.
#. Add a ``TestXXXXX` class and any specific tests to
   ``scipy/stats/tests/test_distributions.py``.
#. Run and pass(!) the tests.

After Implementation
....................

#. Add a tutorial ``doc/source/tutorial/stats/continuous_XXXXX.rst``
#. Add an entry to ``doc/source/tutorial/stats/continuous.rst``.
#. Update the ``number of continuous distributions`` in the example code in ``doc/source/tutorial/stats.rst``.
#. Build the documentation successfully.
#. Submit a PR.

References
^^^^^^^^^^

.. [JKB] Johnson, Kotz, and Balakrishnan, "Continuous Univariate Distributions, Volume 1", Second Edition, John Wiley and Sons,
           p. 173 (1994).

.. _WikipediaDistributions: https://en.wikipedia.org/wiki/List_of_probability_distributions