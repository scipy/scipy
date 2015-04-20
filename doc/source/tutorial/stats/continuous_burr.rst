
.. _continuous-burr:

Burr Distribution
=================

.. math::
   :nowrap:

    \begin{eqnarray*} c & > & 0\\ d & > & 0\\ k & = & \Gamma\left(d\right)\Gamma\left(1-\frac{2}{c}\right)\Gamma\left(\frac{2}{c}+d\right)-\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+d\right)\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;c,d\right) & = & \frac{cd}{x^{c+1}\left(1+x^{-c}\right)^{d+1}}I_{\left(0,\infty\right)}\left(x\right)\\ F\left(x;c,d\right) & = & \left(1+x^{-c}\right)^{-d}\\ G\left(\alpha;c,d\right) & = & \left(\alpha^{-1/d}-1\right)^{-1/c}\\ \mu & = & \frac{\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+d\right)}{\Gamma\left(d\right)}\\ \mu_{2} & = & \frac{k}{\Gamma^{2}\left(d\right)}\\ \gamma_{1} & = & \frac{1}{\sqrt{k^{3}}}\left[2\Gamma^{3}\left(1-\frac{1}{c}\right)\Gamma^{3}\left(\frac{1}{c}+d\right)+\Gamma^{2}\left(d\right)\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(\frac{3}{c}+d\right)\right.\\  &  & \left.-3\Gamma\left(d\right)\Gamma\left(1-\frac{2}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+d\right)\Gamma\left(\frac{2}{c}+d\right)\right]\\ \gamma_{2} & = & -3+\frac{1}{k^{2}}\left[6\Gamma\left(d\right)\Gamma\left(1-\frac{2}{c}\right)\Gamma^{2}\left(1-\frac{1}{c}\right)\Gamma^{2}\left(\frac{1}{c}+d\right)\Gamma\left(\frac{2}{c}+d\right)\right.\\  &  & -3\Gamma^{4}\left(1-\frac{1}{c}\right)\Gamma^{4}\left(\frac{1}{c}+d\right)+\Gamma^{3}\left(d\right)\Gamma\left(1-\frac{4}{c}\right)\Gamma\left(\frac{4}{c}+d\right)\\  &  & \left.-4\Gamma^{2}\left(d\right)\Gamma\left(1-\frac{3}{c}\right)\Gamma\left(1-\frac{1}{c}\right)\Gamma\left(\frac{1}{c}+d\right)\Gamma\left(\frac{3}{c}+d\right)\right]\\ m_{d} & = & \left(\frac{cd-1}{c+1}\right)^{1/c}\,\mathrm{if }cd>1\,\mathrm{otherwise }0\\ m_{n} & = & \left(2^{1/d}-1\right)^{-1/c}\end{eqnarray*}

Implementation: `scipy.stats.burr`
