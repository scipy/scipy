
.. _continuous-laplace:

Laplace (Double Exponential, Bilateral Exponential) Distribution
================================================================

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{2}e^{-\left|x\right|}\\ F\left(x\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}e^{x} &  & x\leq0\\ 1-\frac{1}{2}e^{-x} &  & x>0\end{array}\right.\\ G\left(q\right) & = & \left\{ \begin{array}{ccc} \log\left(2q\right) &  & q\leq\frac{1}{2}\\ -\log\left(2-2q\right) &  & q>\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\ \mu_{2} & = & 2\\ \gamma_{1} & = & 0\\ \gamma_{2} & = & 3\end{eqnarray*}

The ML estimator of the location parameter is

.. math::

     \hat{L}=\mathrm{median}\left(X_{i}\right)

where :math:`X_{i}` is a sequence of :math:`N` mutually independent Laplace RV's and the median is some number
between the :math:`\frac{1}{2}N\mathrm{th}` and the :math:`(N/2+1)\mathrm{th}` order statistic ( *e.g.* take the average of these two) when :math:`N` is even. Also,

.. math::

     \hat{S}=\frac{1}{N}\sum_{j=1}^{N}\left|X_{j}-\hat{L}\right|.

Replace :math:`\hat{L}` with :math:`L` if it is known. If :math:`L` is known then this estimator is distributed as :math:`\left(2N\right)^{-1}S\cdot\chi_{2N}^{2}` .

.. math::
   :nowrap:

    \begin{eqnarray*} h\left[X\right] & = & \log\left(2e\right)\\  & \approx & 1.6931471805599453094.\end{eqnarray*}

Implementation: `scipy.stats.laplace`
