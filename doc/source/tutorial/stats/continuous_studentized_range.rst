.. _continuous-studentized_range:

Studentized Range Distribution
==============================
This distribution has two shape parameters :math:`\nu>0.` and :math:`k>0` and support is :math:`q \geq 0.`

.. math::
   :nowrap:

    \begin{eqnarray*} f(q; k, \nu) = \frac{k(k-1)\nu^{\nu/2}}{\Gamma(\nu/2)2^{\nu/2-1}}\int_{0}^{\infty}\int_{-\infty}^{\infty}s^{\nu-1}e^{\nu s^2/2}s\phi(z)\phi(sq + z)[\Phi(sq + z) - \Phi(z)]^{k-2}\mathrm{d}z\mathrm{d}s \end{eqnarray*}



Implementation: `scipy.stats.studentized_range`
