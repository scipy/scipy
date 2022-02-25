
.. _continuous-nct:

Noncentral t Distribution
=========================

The distribution of the ratio

.. math::

     \frac{U+\lambda}{\chi_{\nu}/\sqrt{\nu}}

where :math:`U` and :math:`\chi_{\nu}` are independent and distributed as a standard normal and chi with :math:`\nu` degrees of freedom. Note :math:`\lambda>0` and :math:`\nu>0` .

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\lambda,\nu\right) & = & \frac{\nu^{\nu/2}\Gamma\left(\nu+1\right)}{2^{\nu}e^{\lambda^{2}/2}\left(\nu+x^{2}\right)^{\nu/2}\Gamma\left(\nu/2\right)}\\
     &  & \times\left\{ \frac{\sqrt{2}\lambda x\,_{1}F_{1}\left(\frac{\nu}{2}+1;\frac{3}{2};\frac{\lambda^{2}x^{2}}{2\left(\nu+x^{2}\right)}\right)}{\left(\nu+x^{2}\right)\Gamma\left(\frac{\nu+1}{2}\right)}\right.\\
     &  & -\left.\frac{\,_{1}F_{1}\left(\frac{\nu+1}{2};\frac{1}{2};\frac{\lambda^{2}x^{2}}{2\left(\nu+x^{2}\right)}\right)}{\sqrt{\nu+x^{2}}\Gamma\left(\frac{\nu}{2}+1\right)}\right\} \\
     & = & \frac{\Gamma\left(\nu+1\right)}{2^{\left(\nu-1\right)/2}\sqrt{\pi\nu}\Gamma\left(\nu/2\right)}\exp\left[-\frac{\nu\lambda^{2}}{\nu+x^{2}}\right]\\
     &  & \times\left(\frac{\nu}{\nu+x^{2}}\right)^{\left(\nu-1\right)/2}Hh_{\nu}\left(-\frac{\lambda x}{\sqrt{\nu+x^{2}}}\right)\\
     F\left(x;\lambda,\nu\right) & = & \left\{
                                  \begin{array}{cc}
                                    {\tilde{F}}_{{\nu ,\mu }}(x) & x\geq0 \\
                                    1 - {\tilde{F}}_{{\nu ,-\mu }}(x) & x<0
                                    \end{array}
                                 \right. \\
    \text{where} \\
     {\tilde{F}}_{{\nu ,\mu }}(x) & = & \Phi (-\mu )+{\frac{1}{2}}\sum _{{j=0}}^{\infty }\left[p_{j}I_{y}\left(j+{\frac{1}{2}},{\frac{\nu }{2}}\right)+q_{j}I_{y}\left(j+1,{\frac{\nu }{2}}\right)\right]\\
     y & = & \frac{x^2}{x^2+\nu}\\
     p_{j} & = & \frac{e^{\left( -\frac{\mu^2}{2} \right)} }{j!} \left(\frac{\mu^2}{2}\right)^{j}\\
     q_{j} & = & {\frac{\mu e^{\left( -\frac{\mu^2}{2} \right)} } {\sqrt{2}\Gamma(j+3/2)}} \left({\frac{\mu^2}{2}}\right)^{j} \end{eqnarray*}

where :math:`I_{y}(a,b)` is the regularized incomplete beta function and
Airy's Hh function is :math:`Hh_{\nu}(x)=\frac{1}{\Gamma(\nu+1)}\int_0^\infty t^\nu e^{\frac{-(t+x)^2}{2}}dt`.

Implementation: `scipy.stats.nct`
