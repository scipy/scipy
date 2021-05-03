.. _continuous-studentized_range:

Studentized Range Distribution
==============================
This distribution has two shape parameters :math:`\nu>0.` and :math:`k>0` and support is :math:`q \geq 0.`

.. math::
   :nowrap:

    \begin{eqnarray*}
    f(q; k, \nu) = \frac{k(k-1)\nu^{\nu/2}}{\Gamma(\nu/2)2^{\nu/2-1}}
    \int_{0}^{\infty} \int_{-\infty}^{\infty} s e^{-\nu s^2/2} \phi(z) \phi(sq + z)
    [\Phi(sq + z) - \Phi(z)]^{k-2} \,dz \,ds
    \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*}
    F(q; k, \nu) = \frac{k\nu^{\nu/2}}{\Gamma(\nu/2)2^{\nu/2-1}}
    \int_{0}^{\infty} \int_{-\infty}^{\infty} s^{\nu-1} e^{-\nu s^2/2} \phi(z)
    [\Phi(sq + z) - \Phi(z)]^{k-1} \,dz \,ds
    \end{eqnarray*}

Note: :math:`\phi(z)` and :math:`\Phi(z)` represent the normal PDF and normal CDF, respectively.

When :math:`\nu` is sufficiently large, the asymptopic approximation of :math:`F(q; k, \nu)` is used:

.. math::
   :nowrap:

    \begin{eqnarray*}
    F(q; k, \nu) = k \int_{-\infty}^{\infty} \phi(z)
    [\Phi(q + z) - \Phi(z)]^{k-1} \,dz
    \end{eqnarray*}


Implementation: `scipy.stats.studentized_range`
