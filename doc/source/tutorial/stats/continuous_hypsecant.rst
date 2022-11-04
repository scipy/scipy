
.. _continuous-hypsecant:

Hyperbolic Secant Distribution
==============================

Related to the logistic distribution and used in lifetime analysis.
Standard form is (defined over all :math:`x` )

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x\right) & = & \frac{1}{\pi}\mathrm{sech}\left(x\right)\\
    F\left(x\right) & = & \frac{2}{\pi}\arctan\left(e^{x}\right)\\
    G\left(q\right) & = & \log\left(\tan\left(\frac{\pi}{2}q\right)\right)\end{eqnarray*}

.. math::

     M\left(t\right)=\sec\left(\frac{\pi}{2}t\right)

.. math::
   :nowrap:

   \begin{eqnarray*} \mu_{n}^{\prime} & = & \frac{1+\left(-1\right)^{n}}{2\pi2^{2n}}n!\left[\zeta\left(n+1,\frac{1}{4}\right)-\zeta\left(n+1,\frac{3}{4}\right)\right]\\
    & = & \left\{
      \begin{array}{cc}
        0 & n \text{ odd}\\
        C_{n/2}\frac{\pi^{n}}{2^{n}} & n \text{ even}
      \end{array}
    \right.\end{eqnarray*}

where :math:`C_{m}` is an integer given by

.. math::
   :nowrap:

    \begin{eqnarray*} C_{m} & = & \frac{\left(2m\right)!\left[\zeta\left(2m+1,\frac{1}{4}\right)-\zeta\left(2m+1,\frac{3}{4}\right)\right]}{\pi^{2m+1}2^{2m}}\\
     & = & 4\left(-1\right)^{m-1}\frac{16^{m}}{2m+1}B_{2m+1}\left(\frac{1}{4}\right)\end{eqnarray*}

where :math:`B_{2m+1}\left(\frac{1}{4}\right)` is the Bernoulli polynomial of order :math:`2m+1` evaluated at :math:`1/4.` Thus

.. math::

     \mu_{n}^{\prime}=\left\{
       \begin{array}{cc}
        0 & n \text{ odd}\\
        4\left(-1\right)^{n/2-1}\frac{\left(2\pi\right)^{n}}{n+1}B_{n+1}\left(\frac{1}{4}\right) & n \text{ even}
      \end{array}
      \right.

.. math::
   :nowrap:

    \begin{eqnarray*} m_{d}=m_{n}=\mu & = & 0\\
    \mu_{2} & = & \frac{\pi^{2}}{4}\\
    \gamma_{1} & = & 0\\
    \gamma_{2} & = & 2\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(2\pi\right).

Implementation: `scipy.stats.hypsecant`
