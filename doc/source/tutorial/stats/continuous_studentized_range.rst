
.. _continuous-studentized_range:

Studentized Range Distribution
==============================
This distribution has two shape parameters :math:`\nu>0` and :math:`k>0` and support is :math:`q\geq0`.

.. math::
   :nowrap:

    \begin{eqnarray*}
        f(q;k,\nu) = \frac{\nu^{\nu/2}}{\Gamma(\nu/2)\,2^{\nu/2-1}}\,  \int_0^\infty\int_{-\infty}^\infty  s^{\nu-1}\,e^{-\nu s^{2}/2}  \,sk(k-1) \varphi(z) \varphi(sq+z) [\Phi(sq+z)-\Phi(z)]^{k-2} \,\mathrm{d}z \, \mathrm{d}s\\

        F(q;k,\nu) = \frac{\nu^{\nu/2}}{\Gamma(\nu/2)\,2^{\nu/2-1}}\, \int_0^\infty\int_{-\infty}^\infty s^{\nu-1}\,e^{-\nu s^{2}/2}  k\varphi(z) [\Phi(sq+z)-\Phi(z)]^{k-1} \mathrm{d}z  \, \mathrm{d}s\\

    \end{eqnarray*}

Implementation: `scipy.stats.studentized_range`
