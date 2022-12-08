.. _continuous-studentized_range:

Studentized Range Distribution
==============================
This distribution has two shape parameters, :math:`k>1` and :math:`\nu>0`, and the support is :math:`x \geq 0`.

.. math::
   :nowrap:

    \begin{eqnarray*}
    f(x; k, \nu) = \frac{k(k-1)\nu^{\nu/2}}{\Gamma(\nu/2)2^{\nu/2-1}}
    \int_{0}^{\infty} \int_{-\infty}^{\infty} s^{\nu} e^{-\nu s^2/2} \phi(z) \phi(sx + z)
    [\Phi(sx + z) - \Phi(z)]^{k-2} \,dz \,ds
    \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*}
    F(q; k, \nu) = \frac{k\nu^{\nu/2}}{\Gamma(\nu/2)2^{\nu/2-1}}
    \int_{0}^{\infty} \int_{-\infty}^{\infty} s^{\nu-1} e^{-\nu s^2/2} \phi(z)
    [\Phi(sq + z) - \Phi(z)]^{k-1} \,dz \,ds
    \end{eqnarray*}

Note: :math:`\phi(z)` and :math:`\Phi(z)` represent the normal PDF and normal CDF, respectively.

When :math:`\nu` exceeds 100,000, the asymptotic approximation of :math:`F(x; k, \nu=\infty)` or :math:`f(x; k, \nu=\infty)` is used:

.. math::
   :nowrap:

    \begin{eqnarray*}
    F(x; k, \nu=\infty) = k \int_{-\infty}^{\infty} \phi(z)
    [\Phi(x + z) - \Phi(z)]^{k-1} \,dz
    \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*}
    f(x; k, \nu=\infty) = k(k-1) \int_{-\infty}^{\infty} \phi(z)\phi(x + z)
    [\Phi(x + z) - \Phi(z)]^{k-2} \,dz
    \end{eqnarray*}


Implementation: `scipy.stats.studentized_range`
