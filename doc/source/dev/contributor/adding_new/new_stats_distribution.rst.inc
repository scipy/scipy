Adding A New Statistics Distribution
------------------------------------

For hundreds of years statisticians, mathematicians and scientists have needed
to understand, analyze and model data.
This has led to a plethora of statistics distributions, many of which are
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

There are a few steps to be done to add a continuous distribution to SciPy.
(Adding a discrete distribution is similar).
We'll use the fictitious "Squirrel" distribution in the instructions below.


Before Implementation
.....................

1. See if ``Squirrel`` has already been implemented--that saves a lot
   of effort!

   - It may have been implemented with a different name.
   - It may have been implemented with a different parameterization
     (shape parameters).
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

2. Create a `SciPy Issue on GitHub <https://github.com/scipy/scipy/issues>`_,
   listing the distribution, references and reasons for its inclusion.


Implementation
..............

#. Find an already existing distribution similar to ``Squirrel``.
   Use its code as a template for ``Squirrel``.
#. Read the docstring for class ``rv_continuous`` in
   `scipy/stats/_distn_infrastructure.py <https://github.com/scipy/scipy/blob/main/scipy/stats/_distn_infrastructure.py#L1378>`_.
#. Write the new code for class ``squirrel_gen`` and insert it into
   `scipy/stats/_continuous_distns.py <https://github.com/scipy/scipy/blob/main/scipy/stats/_continuous_distns.py>`_,
   which is in (mostly) alphabetical order by distribution name.
#. Does the distribution have infinite support? If not, left and/or right
   endpoints ``a``, ``b`` need to be specified in the call to ``squirrel_gen(name='squirrel', a=?, b=?)``.
#. If the support depends upon the shape parameters,
   ``squirrel_gen._get_support()`` needs to be implemented.
#. The default inherited ``_argcheck()`` implementation checks that the shape parameters
   are positive.  Create a more appropriate implementation.
#. If ``squirrel_gen.ppf()`` is expensive to compute relative to
   ``squirrel_gen.pdf()``, consider setting the ``momtype`` in the call
   to ``squirrel_gen()``.
#. If ``squirrel_gen.rvs()`` is expensive to compute, consider implementing a
   specific ``squirrel_gen._rvs()``.
#. Add the name to the listing in the docstring of
   `scipy/stats/__init__.py <https://github.com/scipy/scipy/blob/main/scipy/stats/__init__.py>`_.
#. Add the name and a good set of example shape parameters to the ``distcont``
   list in
   `scipy/stats/_distr_params.py <https://github.com/scipy/scipy/blob/main/scipy/stats/_distr_params.py#L5>`_.
   These shape parameters are used both for testing and automatic documentation generation.
#. Add the name and an _invalid_ set of example shape parameters to the
   list in ``invdistcont``, also in
   `_distr_params.py <https://github.com/scipy/scipy/blob/main/scipy/stats/_distr_params.py>`_.
   These shape parameters are also used for testing.
#. Add a ``TestSquirrel`` class and any specific tests to
   `scipy/stats/tests/test_distributions.py <https://github.com/scipy/scipy/blob/main/scipy/stats/tests/test_distributions.py>`_.
#. Run and pass(!) the tests.

After Implementation
....................

#. Add a tutorial ``doc/source/tutorial/stats/continuous_squirrel.rst``
#. Add it to the listing of continuous distributions in
   `doc/source/tutorial/stats/continuous.rst <https://github.com/scipy/scipy/blob/main/doc/source/tutorial/stats/continuous.rst>`_.
#. Update the ``number of continuous distributions`` in the example code in
   `doc/source/tutorial/stats.rst <https://github.com/scipy/scipy/blob/main/doc/source/tutorial/stats.rst>`_.
#. Build the documentation successfully.
#. Submit a PR.

References
^^^^^^^^^^

.. [JKB] Johnson, Kotz, and Balakrishnan, "Continuous Univariate Distributions, Volume 1", Second Edition, John Wiley and Sons,
           p. 173 (1994).

.. _WikipediaDistributions: https://en.wikipedia.org/wiki/List_of_probability_distributions
