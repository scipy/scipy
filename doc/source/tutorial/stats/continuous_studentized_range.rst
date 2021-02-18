
.. _continuous-studentized_range:

Studentized Range Distribution
==============================
This distribution has two shape parameters :math:`\nu>0` and :math:`k>0`.

.. math::
   :nowrap:

    \begin{eqnarray*}
        f(q;k,\nu)=\frac{\,\sqrt{2\pi\,}\,k\,(k-1)\,\nu^{\nu/2}\,}{\Gamma(\nu/2)\,2^{\left(\nu/2-1\right)}}  \int_0^\infty\int_{-\infty}^\infty s^\nu\,\varphi(\sqrt{\nu\,}\,s)\,  \varphi(z+q\,s)\,\varphi(z)\, \left[\Phi(z+q\,s)-\Phi(z)\right]^{k-2} \, \mathrm{d}z \, \mathrm{d}s\\

        F(q;k,\nu) = \frac{\sqrt{2\pi\,}\,k\,\nu^{\nu/2}}{\,\Gamma(\nu/2)\,2^{(\nu/2-1)}\,} \int_0^\infty\int_{-\infty}^\infty s^{\nu-1} \varphi(\sqrt{\nu\,}\,s)  \varphi(z) \left[\Phi(z+q\,s)-\Phi(z)\right]^{k-1} \, \mathrm{d}z  \, \mathrm{d}s\\

    \end{eqnarray*}

Note that in :math:`f(q;k,\nu)` the substition :math:`\varphi(\sqrt{\nu\,}\,s) \, \sqrt{2\pi\,} = e^{-\left(\nu\, s^2/2\right)}` was used to remove an exponential term.

Implementation: `scipy.stats.studentized_range`
