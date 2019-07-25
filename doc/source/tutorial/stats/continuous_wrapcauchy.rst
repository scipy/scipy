
.. _continuous-wrapcauchy:

Wrapped Cauchy Distribution
===========================

Thre is one shape parameter :math:`c\in\left(0,1\right)` with support :math:`x\in\left[0,2\pi\right]`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c\right) & = & \frac{1-c^{2}}{2\pi\left(1+c^{2}-2c\cos x\right)}\\
    g_{c}\left(x\right) & = & \frac{1}{\pi}\arctan\left(\frac{1+c}{1-c}\tan\left(\frac{x}{2}\right)\right)\\
    r_{c}\left(q\right) & = & 2\arctan\left(\frac{1-c}{1+c}\tan\left(\pi q\right)\right)\\
    F\left(x;c\right) & = & \left\{
          \begin{array}{ccc}
            g_{c}\left(x\right) &  & 0\leq x<\pi\\
            1-g_{c}\left(2\pi-x\right) &  & \pi\leq x\leq2\pi
          \end{array}
          \right.\\
   G\left(q;c\right) & = & \left\{
          \begin{array}{ccc}
            r_{c}\left(q\right) &  & 0\leq q<\frac{1}{2}\\
            2\pi-r_{c}\left(1-q\right) &  & \frac{1}{2}\leq q\leq1
          \end{array}
          \right.\end{eqnarray*}

.. math::

     h\left[X\right]=\log\left(2\pi\left(1-c^{2}\right)\right).

Implementation: `scipy.stats.wrapcauchy`
