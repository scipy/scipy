
.. _continuous-bradford:

Bradford Distribution
=====================

There is one shape parameter, :math:`c>0`, and the support is :math:`x\in [0,1]`.

.. math::
   :nowrap:

    \begin{eqnarray*} \textrm{Let } k & = & \log\left(1+c\right)\\
    \textrm{Then}\\
    f\left(x;c\right) & = & \frac{c}{k\left(1+cx\right)}\\
    F\left(x;c\right) & = & \frac{\log\left(1+cx\right)}{k}\\
    G\left(q\; c\right) & = & \frac{\left(1+c\right)^{q}-1}{c}\\
    M\left(t\right) & = & \frac{1}{k}e^{-t/c}\left[\mathrm{Ei}\left(t+\frac{t}{c}\right)-\mathrm{Ei}\left(\frac{t}{c}\right)\right]\\
    \mu & = & \frac{c-k}{ck}\\
    \mu_{2} & = & \frac{\left(c+2\right)k-2c}{2ck^{2}}\\
    \gamma_{1} & = & \frac{\sqrt{2}\left(12c^{2}-9kc\left(c+2\right)+2k^{2}\left(c\left(c+3\right)+3\right)\right)}{\sqrt{c\left(c\left(k-2\right)+2k\right)}\left(3c\left(k-2\right)+6k\right)}\\
    \gamma_{2} & = & \frac{c^{3}\left(k-3\right)\left(k\left(3k-16\right)+24\right)+12kc^{2}\left(k-4\right)\left(k-3\right)+6ck^{2}\left(3k-14\right)+12k^{3}}{3c\left(c\left(k-2\right)+2k\right)^{2}}\\
    m_{d} & = & 0\\
    m_{n} & = & \sqrt{1+c}-1\\
    h\left[X\right]& = & \frac{1}{2}\log\left(1+c\right)-\log\left(\frac{c}{\log\left(1+c\right)}\right)\end{eqnarray*}

where :math:`\mathrm{Ei}\left(\mathrm{z}\right)` is the exponential integral function.

Implementation: `scipy.stats.bradford`
