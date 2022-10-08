
.. _continuous-alpha:

Alpha Distribution
==================

One shape parameter :math:`\alpha>0` (parameter :math:`\beta` in DATAPLOT
is a scale-parameter). The support for the standard form is :math:`x>0`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\alpha\right) & = & \frac{1}{x^{2}\Phi\left(\alpha\right)\sqrt{2\pi}}\exp\left(-\frac{1}{2}\left(\alpha-\frac{1}{x}\right)^{2}\right)\\ F\left(x;\alpha\right) & = & \frac{\Phi\left(\alpha-\frac{1}{x}\right)}{\Phi\left(\alpha\right)}\\ G\left(q;\alpha\right) & = & \left[\alpha-\Phi^{-1}\left(q\Phi\left(\alpha\right)\right)\right]^{-1}\end{eqnarray*}

.. math::

     M\left(t\right)=\frac{1}{\Phi\left(a\right)\sqrt{2\pi}}\int_{0}^{\infty}\frac{e^{xt}}{x^{2}}\exp\left(-\frac{1}{2}\left(\alpha-\frac{1}{x}\right)^{2}\right)dx

No moments?

.. math::

     l_{\mathbf{x}}\left(\alpha\right)=N\log\left[\Phi\left(\alpha\right)\sqrt{2\pi}\right]+2N\overline{\log\mathbf{x}}+\frac{N}{2}\alpha^{2}-\alpha\overline{\mathbf{x}^{-1}}+\frac{1}{2}\overline{\mathbf{x}^{-2}}

Implementation: `scipy.stats.alpha`
