
.. _continuous-truncnorm:

Truncated Normal Distribution
=============================

A normal distribution restricted to lie within a certain range given
by two parameters :math:`A` and :math:`B` . Notice that this :math:`A` and :math:`B` correspond to the bounds on :math:`x` in standard form. For :math:`x\in\left[A,B\right]` we get

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;A,B\right) & = & \frac{\phi\left(x\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\\
    F\left(x;A,B\right) & = & \frac{\Phi\left(x\right)-\Phi\left(A\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\\
    G\left(q;A,B\right) & = & \Phi^{-1}\left(q\Phi\left(B\right)+\Phi\left(A\right)\left(1-q\right)\right)\end{eqnarray*}

where

.. math::
   :nowrap:

    \begin{eqnarray*} \phi\left(x\right) & = & \frac{1}{\sqrt{2\pi}}e^{-x^{2}/2}\\
    \Phi\left(x\right) & = & \int_{-\infty}^{x}\phi\left(u\right)du.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} \mu & = & \frac{\phi\left(A\right)-\phi\left(B\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\\
    \mu_{2} & = & 1+\frac{A\phi\left(A\right)-B\phi\left(B\right)}{\Phi\left(B\right)-\Phi\left(A\right)}-\left(\frac{\phi\left(A\right)-\phi\left(B\right)}{\Phi\left(B\right)-\Phi\left(A\right)}\right)^{2}\end{eqnarray*}

Implementation: `scipy.stats.truncnorm`
