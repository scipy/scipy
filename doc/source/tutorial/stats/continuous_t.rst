
.. _continuous-t:

Student t Distribution
======================

There is one shape parameter :math:`\nu>0` and the support is :math:`x\in\mathbb{R}`.

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu\right) & = & \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\pi\nu}\Gamma\left(\frac{\nu}{2}\right)\left[1+\frac{x^{2}}{\nu}\right]^{\frac{\nu+1}{2}}}\\
    F\left(x;\nu\right) & = &
      \left\{
        \begin{array}{ccc}
          \frac{1}{2}I\left(\frac{\nu}{\nu+x^{2}}; \frac{\nu}{2},\frac{1}{2}\right) &  & x\leq0\\
          1-\frac{1}{2}I\left(\frac{\nu}{\nu+x^{2}}; \frac{\nu}{2},\frac{1}{2}\right) &  & x\geq0
        \end{array}
      \right.\\
    G\left(q;\nu\right) & = & \left\{
      \begin{array}{ccc}
        -\sqrt{\frac{\nu}{I^{-1}\left(2q; \frac{\nu}{2},\frac{1}{2}\right)}-\nu} &  & q\leq\frac{1}{2}\\
        \sqrt{\frac{\nu}{I^{-1}\left(2-2q; \frac{\nu}{2},\frac{1}{2}\right)}-\nu} &  & q\geq\frac{1}{2}
      \end{array}
      \right. \end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{n}=m_{d}=\mu & = & 0\\
    \mu_{2} & = & \frac{\nu}{\nu-2}\quad\nu>2\\
    \gamma_{1} & = & 0\quad\nu>3\\
    \gamma_{2} & = & \frac{6}{\nu-4}\quad\nu>4\end{eqnarray*}

where :math:`I\left(x; a,b\right)` is the incomplete beta integral and :math:`I^{-1}\left(I\left(x; a,b\right); a,b\right)=x`.
As :math:`\nu\rightarrow\infty,` this distribution approaches the standard normal distribution.

.. math::

     h\left[X\right]=\frac{\nu+1}{2} \left[\psi \left(\frac{1+\nu}{2} \right) -\psi \left(\frac{\nu}{2} \right) \right] + \ln \left[ \sqrt{\nu} B \left( \frac{\nu}{2}, \frac{1}{2} \right) \right]

where :math:`\psi(x)` is the digamma function and :math:`B(x, y)` is the
beta function.

References
----------

- "Student's t-distribution", Wikipedia, https://en.wikipedia.org/wiki/Student%27s_t-distribution

Implementation: `scipy.stats.t`
