
.. _continuous-halflogistic:

Half-Logistic Distribution
==========================

In the limit as :math:`c\rightarrow\infty` for the generalized half-logistic we have the half-logistic defined
over :math:`x\geq0.` Also, the distribution of :math:`\left|X\right|` where :math:`X` has logistic distribution.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{2e^{-x}}{\left(1+e^{-x}\right)^{2}}=\frac{1}{2}\mathrm{sech}^{2}\left(\frac{x}{2}\right)\\ F\left(x\right) & = & \frac{1-e^{-x}}{1+e^{-x}}=\tanh\left(\frac{x}{2}\right)\\ G\left(q\right) & = & \log\left(\frac{1+q}{1-q}\right)=2\mathrm{arctanh}\left(q\right)\end{eqnarray*}

.. math::

     M\left(t\right)=1-t\psi_{0}\left(\frac{1}{2}-\frac{t}{2}\right)+t\psi_{0}\left(1-\frac{t}{2}\right)

.. math::

     \mu_{n}^{\prime}=2\left(1-2^{1-n}\right)n!\zeta\left(n\right)\quad n\neq1

.. math::
   :nowrap:

    \begin{eqnarray*} \mu_{1}^{\prime} & = & 2\log\left(2\right)\\ \mu_{2}^{\prime} & = & 2\zeta\left(2\right)=\frac{\pi^{2}}{3}\\ \mu_{3}^{\prime} & = & 9\zeta\left(3\right)\\ \mu_{4}^{\prime} & = & 42\zeta\left(4\right)=\frac{7\pi^{4}}{15}\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & 2-\log\left(2\right)\\  & \approx & 1.3068528194400546906.\end{eqnarray*}

Implementation: `scipy.stats.halflogistic`
