
.. _continuous-nct:

Noncentral t Distribution
=========================

The distribution of the ratio

.. math::

     \frac{U+\lambda}{\chi_{\nu}/\sqrt{\nu}}

where :math:`U` and :math:`\chi_{\nu}` are independent and distributed as a standard normal and chi with :math:`\nu` degrees of freedom. Note :math:`\lambda>0` and :math:`\nu>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda,\nu\right) & = & \frac{\nu^{\nu/2}\Gamma\left(\nu+1\right)}{2^{\nu}e^{\lambda^{2}/2}\left(\nu+x^{2}\right)^{\nu/2}\Gamma\left(\nu/2\right)}\\  &  & \times\left\{ \frac{\sqrt{2}\lambda x\,_{1}F_{1}\left(\frac{\nu}{2}+1;\frac{3}{2};\frac{\lambda^{2}x^{2}}{2\left(\nu+x^{2}\right)}\right)}{\left(\nu+x^{2}\right)\Gamma\left(\frac{\nu+1}{2}\right)}\right.\\  &  & -\left.\frac{\,_{1}F_{1}\left(\frac{\nu+1}{2};\frac{1}{2};\frac{\lambda^{2}x^{2}}{2\left(\nu+x^{2}\right)}\right)}{\sqrt{\nu+x^{2}}\Gamma\left(\frac{\nu}{2}+1\right)}\right\} \\  & = & \frac{\Gamma\left(\nu+1\right)}{2^{\left(\nu-1\right)/2}\sqrt{\pi\nu}\Gamma\left(\nu/2\right)}\exp\left[-\frac{\nu\lambda^{2}}{\nu+x^{2}}\right]\\  &  & \times\left(\frac{\nu}{\nu+x^{2}}\right)^{\left(\nu-1\right)/2}Hh_{\nu}\left(-\frac{\lambda x}{\sqrt{\nu+x^{2}}}\right)\\ F\left(x;\lambda,\nu\right) & =\end{eqnarray*}

Implementation: `scipy.stats.nct`
