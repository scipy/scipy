
.. _continuous-t:

Student t Distribution
======================

Shape parameter :math:`\nu>0.` :math:`I\left(a,b,x\right)` is the incomplete beta integral and :math:`I^{-1}\left(a,b,I\left(a,b,x\right)\right)=x`

.. math::
   :nowrap:

    \begin{eqnarray*} f\left(x;\nu\right) & = & \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\sqrt{\pi\nu}\Gamma\left(\frac{\nu}{2}\right)\left[1+\frac{x^{2}}{\nu}\right]^{\frac{\nu+1}{2}}}\\ F\left(x;\nu\right) & = & \left\{ \begin{array}{ccc} \frac{1}{2}I\left(\frac{\nu}{2},\frac{1}{2},\frac{\nu}{\nu+x^{2}}\right) &  & x\leq0\\ 1-\frac{1}{2}I\left(\frac{\nu}{2},\frac{1}{2},\frac{\nu}{\nu+x^{2}}\right) &  & x\geq0\end{array}\right.\\ G\left(q;\nu\right) & = & \left\{ \begin{array}{ccc} -\sqrt{\frac{\nu}{I^{-1}\left(\frac{\nu}{2},\frac{1}{2},2q\right)}-\nu} &  & q\leq\frac{1}{2}\\ \sqrt{\frac{\nu}{I^{-1}\left(\frac{\nu}{2},\frac{1}{2},2-2q\right)}-\nu} &  & q\geq\frac{1}{2}\end{array}\right.\end{eqnarray*}

.. math::
   :nowrap:

    \begin{eqnarray*} m_{n}=m_{d}=\mu & = & 0\\ \mu_{2} & = & \frac{\nu}{\nu-2}\quad\nu>2\\ \gamma_{1} & = & 0\quad\nu>3\\ \gamma_{2} & = & \frac{6}{\nu-4}\quad\nu>4\end{eqnarray*}

As :math:`\nu\rightarrow\infty,` this distribution approaches the standard normal distribution.

.. math::

     h\left[X\right]=\frac{1}{4}\log\left(\frac{\pi c\Gamma^{2}\left(\frac{c}{2}\right)}{\Gamma^{2}\left(\frac{c+1}{2}\right)}\right)-\frac{\left(c+1\right)}{4}\left[\Psi\left(\frac{c}{2}\right)-cZ\left(c\right)+\pi\tan\left(\frac{\pi c}{2}\right)+\gamma+2\log2\right]

where

.. math::

     Z\left(c\right)=\,_{3}F_{2}\left(1,1,1+\frac{c}{2};\frac{3}{2},2;1\right)=\sum_{k=0}^{\infty}\frac{k!}{k+1}\frac{\Gamma\left(\frac{c}{2}+1+k\right)}{\Gamma\left(\frac{c}{2}+1\right)}\frac{\Gamma\left(\frac{3}{2}\right)}{\Gamma\left(\frac{3}{2}+k\right)}

Implementation: `scipy.stats.t`
