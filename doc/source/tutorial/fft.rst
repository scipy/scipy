.. _tutorial_FFT:

**********************************************
Discrete Fourier Transforms (:mod:`scipy.fft`)
**********************************************
.. sectionauthor:: SciPy Developers

.. currentmodule:: scipy.fft

.. include:: _LaTeXMathMacros.rst

Fourier analysis is a method for expressing a function as a sum of periodic
components, and for recovering the signal from those components. When both
the function and its Fourier transform are replaced with discretized
counterparts, it is called the discrete Fourier transform (DFT). The DFT has
become a mainstay of numerical computing in part because of a very fast
algorithm for computing it, called the Fast Fourier Transform (FFT), which was
known to Gauss (1805) and was brought to light in its current form by Cooley
and Tukey [#CT65]_. Press et al. [#NR07]_ provide an accessible introduction to
Fourier analysis and its applications.

This chapter is structured as follows: The first section defines the discrete Fourier
transform and presents its key properties. The second section discusses the
relationship to continuous Fourier series and introduces various forms of its Fourier
transform. Then practical implications, which mainly results from the DFT's periodicity
properties, are discussed in the :ref:`tutorial_FFT_Caveats` section. The remaining
sections give definitions of higher-dimensional DFTs, discrete sine and cosine
transforms as well as the Hankel transform.

.. contents::


.. _tutorial_FFT_DFT:

The Discrete Fourier Transform
==============================
The discrete Fourier Transform [#WIKI_DFT]_ is a one-to-one mapping between two
complex-valued sequences of length :math:`N`. The mappings between :math:`x[k],
X[l]\in\IC`, with :math:`k, l \in\{0, \ldots, N-1\}`, are defined by

.. math::
    :label: eq_FFT_DFT

    X[l] = \frac{1}{\gamma} \sum_{k=0}^{N-1} x[n] \e^{-\jj 2 \pi l k / N}\,,\qquad
    x[k] = \frac{\gamma}{N}\sum_{l=0}^{N-1} X[l] \e^{\jj 2 \pi k l / N}\,,

which are implemented by the :func:`~scipy.fft.fft` / :func:`~scipy.fft.ifft`
functions. :math:`\gamma` is a normalization constant, which is set by the `norm`
parameter. It may have the values :math:`1` (``norm='backward'``, default),
:math:`N` (``norm='forward'``), or :math:`\sqrt{N}` (``norm='ortho'``).


Note that :math:`X[l]` can be continued periodically outside its domain :math:`\{0,
\ldots, N-1\}`. I.e., :math:`X[l+pN] = X[l]` for :math:`p\in\IZ`  when :math:`X[l]` is
calculated with Eq. :math:numref:`eq_FFT_DFT`. This allows to define the DFT over
shifted coefficients, e.g., using :math:`l\in\{-N//2, \ldots, (N-1)-N//2\}`. The same
holds for :math:`x[k]`: When determined with Eq. :math:numref:`eq_FFT_DFT`, :math:`x[k
+ pN] = x[k]`, :math:`p \in\IZ` is true.

This DFT can be rewritten in vector-matrix form as

.. math::

    \underbrace{\begin{bmatrix} X[0]\\ X[1]\\ \vdots\\ X[N-1] \end{bmatrix}}_{=:\vb{X}}
      = \frac{1}{\gamma}\underbrace{\begin{bmatrix}
         \e^{-\jj 2 \pi \cdot 0 \cdot 0 / N} & \cdots & \e^{-\jj 2 \pi \cdot 0 \cdot (N-1) / N}\\
         \vdots & \ddots & \vdots\\
         \e^{-\jj 2 \pi (N-1) \cdot 0 / N} & \cdots & \e^{-\jj 2 \pi (N-1) (N-1) / N}
       \end{bmatrix}}_{=:\mF}
    \underbrace{\begin{bmatrix} x[0]\\ x[1]\\ \vdots\\ x[N-1] \end{bmatrix}}_{=:\vb{x}}\,,

with the matrix :math:`\mF\in\IC^{N\times N}` being symmetric, i.e., :math:`\mF^T=\mF`.
It has the elements :math:`F[l, k] = \exp(-\jj2\pi l k / N)` and allows the DFT and its
inverse to be expressed as

.. math::

    \vb{X} = \frac{1}{\gamma}\mF\,\vb{x}\,, \qquad
    \vb{x} = \frac{\gamma}{N}\conj{\mF}\,\vb{X}
           = \frac{\gamma}{N}\conjT{\mF}\,\vb{X} \,,


with :math:`\conjT{\mF}` being the conjugate transpose of :math:`\mF`. It can be shown
that :math:`\conjT{\mF}\mF = N\mI` holds and thus :math:`\mF^{-1} = \conjT{\mF}/N`,
which makes it straightforward to verify that the two relations of Eq.
:math:numref:`eq_FFT_DFT` are inverses of each other. The :mod:`~scipy.linalg` module
provides the function :func:`~scipy.linalg.dft` to calculate the matrix :math:`\mF /
\gamma`:

>>> import numpy as np
>>> from scipy.linalg import dft
>>> from scipy.fft import fft
>>> np.set_printoptions(precision=2, suppress=True)  # for compact output
...
>>> dft(4)
array([[ 1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j],
       [ 1.+0.j,  0.-1.j, -1.-0.j, -0.+1.j],
       [ 1.+0.j, -1.-0.j,  1.+0.j, -1.-0.j],
       [ 1.+0.j, -0.+1.j, -1.-0.j,  0.-1.j]])
>>> fft(np.eye(4)) # produces the same matrix
array([[ 1.-0.j,  1.+0.j,  1.-0.j,  1.-0.j],
       [ 1.-0.j,  0.-1.j, -1.-0.j,  0.+1.j],
       [ 1.-0.j, -1.+0.j,  1.-0.j, -1.-0.j],
       [ 1.-0.j,  0.+1.j, -1.-0.j,  0.-1.j]])

.. A comment for the sake of completeness: Proving that `FF := conj(F) @ F` = N I holds
   amounts to noting that each entry FF[p, q] = sum_k( exp(2jπ k (p-q) / N) ) is a
   finite geometric series, which is straightforward to evaluate.

This notation makes it also easy to show that the DFT changes the inner product (or
scalar product) only by the scaling factor :math:`N/\gamma^2`. I.e., expressing the
inner product of :math:`\vb{X}` and :math:`\vb{Y}` in terms of :math:`\vb{x}` and
:math:`\vb{y}` gives

.. math::

    \langle \vb{X}, \vb{Y}\rangle := \conjT{\vb{X}}\vb{Y} =
     \conjT{\big(\frac{1}{\gamma}\mF\,\vb{x}\big)}
            \big(\frac{1}{\gamma}\mF\,\vb{y}\big) =
    \frac{1}{\gamma^2}\conjT{\vb{x}} \conjT{\mF} \mF \vb{y} =
    \frac{N}{\gamma^2}\conjT{\vb{x}} \vb{y} =
    \frac{N}{\gamma^2} \langle \vb{x}, \vb{y}\rangle\,.

This illustrates that the DFT is a unitary transform for :math:`\gamma = \sqrt{N}`.
The sums of the absolute squares are a special case, i.e.,

.. math::

    \sum_{l=0}^{N-1} \abs{X[l]}^2 =
    \langle \vb{X}, \vb{X}\rangle =
    \frac{N}{\gamma^2} \langle \vb{x}, \vb{x}\rangle =
    \frac{N}{\gamma^2} \sum_{k=0}^{N-1} \abs{x[k]}^2\,.

Note that the sum of the unsquared values can be retrieved from the zeroth DFT
component, i.e.,

.. math::

    X[0] = \frac{1}{\gamma} \sum_{k=0}^{N-1} x[n]\,,\qquad
    x[0] = \frac{\gamma}{N}\sum_{l=0}^{N-1} X[l]\,.

The circular convolution of two sequences of length :math:`N`, i.e.,

.. math::  z[k] = x[k] \ast y[k] := \sum_{q=0}^{N-1} x[q]\, y\big[(k-q) \bmod N\big]

can be expressed as a multiplication of their DFTs and vice versa, i.e.,

.. math::
    :label: eq_FFT_convolution_relation

    \begin{align}
       z[k] &= x[k] \ast y[k] & \Leftrightarrow &&
       Z[l] &= \gamma\, X[l]\, Y[l] \,,\\
       Z[l] &= X[l] \ast Y[l]  & \Leftrightarrow &&
       z[l] &= \frac{N}{\gamma} x[k] \,y[k] \,.
    \end{align}

Here, the circular convolution is defined over finite length sequences, whereas in
literature it is often defined over infinite length sequences. For sufficiently large
:math:`N` (roughly: :math:`N > 500`), it is typically computationally more efficient to
use this multiplication together with the :func:`~scipy.fft.fft` /
:func:`~scipy.fft.ifft` functions than to calculate the convolution directly. Note that
:func:`~scipy.signal.convolve`, :func:`~scipy.signal.fftconvolve` and
:func:`numpy.convolve` implement the non-circular convolution by using appropriate
zero-padding before applying the FFT.

.. Note: Source for the "N > 500" statement is the `scipy.signal.fftconvolve` docstr.


Fourier Series and the DFT
==========================
A Fourier series

.. math::
    :label: eq_FFT_FSeries

    z(t) = \frac{\gamma}{N}
             \sum_{l=-N_2}^{N_2} c_l\, \e^{\jj2\pi l t/\tau}\,, \qquad t\in[0, \tau]\,,

with :math:`N=2N_2+1` complex-valued Fourier coefficients :math:`c_l` represents a
closed curve :math:`z: [0, \tau] \mapsto \IC`. Determining the values of :math:`N`
equidistant samples with sampling interval :math:`T := \tau / N`, i.e.,

.. math::

    z[k] := z(kT) = \frac{\gamma}{N}
        \sum_{l=-N_2}^{N_2} c_l\, \e^{\jj2\pi l k / N}\,, \qquad
        k \in \{0,\ldots,N-1\}\,,

may be interpreted as a shifted version the inverse DFT of Eq.
:math:numref:`eq_FFT_DFT`. Since the (inverse) DFT is a one-to-one mapping, the
:math:`N` samples are an alternative representation of the :math:`N` Fourier
coefficients. Here, the scaling factor :math:`\gamma/N` exists only to ensure
notational consistency with Eq. :math:numref:`eq_FFT_DFT`.

The following example produces a closed curve :math:`z(t)` resembling the letter ω by
utilizing a Fourier series: The Fourier coefficients are calculated by the
:func:`~scipy.fft.fft` function from 31 provided samples. In the plot, the dots
represent the samples :math:`z[k]` and the coloration of the line the value of the
parameter :math:`t\in[0, 1]`. Since the :func:`~scipy.fft.fft` function expects indexes
:math:`\{0,\ldots,N-1\}` but :math:`z(t)` is defined on the indexes
:math:`\{-N_2,\ldots,N_2\}`, the :func:`~scipy.fft.fftshift` function is utilized.

.. plot:: tutorial/examples/fft_ParametricClosedCurve.py

.. Note: The provided samples were determined by:
   1. In Inkscape a bitmap image was traced into a path and saved as a SVG.
   2. The path points were loaded into a numpy array by using the
      package `svgpathtools`.
   3. Downsampling to 31 samples was performed by FFT-based resampling aka
      `signal.resample`.

Remarks:

* That the plot is a closed curve is due to the periodic nature of Fourier series,
  i.e., :math:`z(t+p\tau) = z(t)\,,\ p\in\IZ`. As a consequence, a given closed curve
  can be represented by various Fourier series and with that by various sets of
  samples. E.g., increasing the number of samples can be achieved by adding zero-valued
  Fourier coefficients to their DFT and then calculating the corresponding inverse DFT.
  This technique is used in the :func:`~scipy.signal.resample` function.
* A Fourier series is a decomposition of a complex-valued closed curve into circular
  closed curves with various frequencies: The :math:`l`-th term :math:`z_l(t) :=
  \tfrac{\gamma}{N} c_l\, \exp(\jj2\pi l t / \tau)` may be interpreted as a circular
  curve centered at the origin with radius :math:`|\gamma c_l / N|` and starting point
  :math:`z_l(0)= \gamma c_l / N`. For :math:`t\in[0,\tau]`, it cycles the circle
  :math:`l` times in counterclockwise direction if :math:`l > 0` or clockwise direction
  if :math:`l < 0`. The number of cycles per unit interval, i.e, :math:`f_l := l /
  \tau`, is usually called the frequency of :math:`z_l(t)` (having the same sign as
  :math:`l`). The :func:`~scipy.fft.fftfreq` function determines those corresponding
  frequencies for the Fourier coefficients.
* Especially in signal processing, :math:`t` frequently represents time. Hence, the
  sample values :math:`z[k]` with their associated sample times :math:`t[k] = k T` are
  commonly referred to as being in the so-called "time-domain", whereas coefficients
  :math:`c_l` with their associated frequency values :math:`f_l = l \Delta f` with
  :math:`\Delta f := 1/\tau = 1 / (NT)`, are referred as to being in the so-called
  "frequency domain" or "Fourier domain".
* If only real-valued curves are of interest, then
  :math:`\Im\{z(t)\} = \big(z(t) - \conj{z(t)}\big)/2\jj \equiv 0` holds. This is
  equivalent to :math:`z(t) = \conj{z(t)}` as well as to :math:`c_l = \conj{c_{-l}}`.
  As a consequence, for DFTs with real-valued inputs only the non-negative DFT
  coefficients need to be calculated. This is implemented by the
  :func:`~scipy.fft.rfft` / :func:`~scipy.fft.irfft` / :func:`~scipy.fft.rfftfreq`
  functions.

Besides the DFT, other variations for calculating a Fourier transform of a Fourier
series exist. The important case of extending the domain from the fixed interval
:math:`[0, \tau]` to the real axis :math:`\IR` is discussed in the following: When
assuming :math:`z(t)\equiv 0` holds for :math:`t \not\in[0, \tau]`, it is
straightforward to apply the continuous Fourier transform to the Fourier Series of Eq.
:math:numref:`eq_FFT_FSeries`:

.. math::

    Z(f) &:= \int_\IR z(t)\, \e^{-\jj2\pi t f} \dd t
      = \frac{\gamma}{N}\sum_{l=-N_2}^{N_2}
               c_l \int_0^\tau \exp\!\big(\jj2\pi (l - f\tau) t / \tau\big)\,  \dd t \\
     &=  \frac{\gamma}{N} \sum_{l=-N_2}^{N_2}
                     c_l\frac{\tau (\e^{\jj2\pi (l - f\tau)} - 1)}{\jj2\pi (l - f\tau)}
      = \tau \frac{\gamma}{N}\sum_{l=-N_2}^{N_2}
                                   c_l\, \sinc(l - f\tau)\, \e^{\jj\pi (l - f\tau)} \,,

with :math:`\sinc(f) := \sin(\pi f) / (\pi f)`. Note that :math:`Z(f)` is continuous
and its support extends over the complete real axis. The DFT may be interpreted as
sampling the continuous Fourier transform due to

.. math::

    Z(l \Delta f) = \tau\frac{\gamma}{N} c_l = T \gamma c_l
                                                         \quad\text{ for }\quad l\in\IZ

being true. When continuing :math:`z(t)` periodically outside the interval :math:`[0,
\tau]`, i.e.,

.. math::

    z_p(t) = \sum_{p\in\IZ} z(t + p\tau) \,,

then the continuous Fourier transform can be expressed as

.. math::

    Z_p(f) = \frac{\gamma}{N}\sum_{l=-N_2}^{N_2} c_l\, \delta(f- l \Delta f) \,,\qquad
             \Delta f = 1 / \tau = 1/ (NT) \,,

where :math:`\delta(.)` denotes the Dirac delta function. Note that :math:`Z_p(f)` is
only non-zero at multiples of :math:`\Delta f`. Furthermore, :math:`Z_p(f)` is
bandlimited, i.e., :math:`Z_p(f)\equiv 0` for :math:`|f| > N \Delta f / 2 = 1 / (2T)`.
In signal processing, :math:`f_S := 1/T` is usually referred to as sampling frequency
and :math:`f_\text{Ny} := f_S / 2 = 1 / (2T)` is called the Nyquist frequency.

The so-called "discrete-time Fourier transform" (DTFT) [#WIKI_DTFT]_ of the sampled
signal :math:`z_p[k] := z_p(kT)` is defined as

.. math::

    Z_p^\text{DTFT}(f) := \sum_{k\in\IZ}  z_p[k] \e^{-2\jj\pi f k T}
    \quad\text{ or }\quad
    Z_p^\text{DTFT}(\jj\Omega) := \sum_{k\in\IZ}  z_p[k] \e^{-\jj k \Omega}\,,

with :math:`\Omega:= 2\pi f/f_S = 2\pi f T`. The convention of using the normalized
frequency :math:`\Omega` instead of :math:`f` is common in signal processing.
:math:`Z_p^\text{DTFT}(f)` is a continuous function in :math:`f` with a periodicity of
:math:`f_S`. Due to :math:`z_p[k]` being periodic in :math:`N`, it can be expressed in
terms of the Fourier coefficients as

.. math::

    Z_p^\text{DTFT}(f) = \frac{\gamma}{N T} \sum_{l=-N_2}^{N_2} \sum_{p\in\IZ}
                                      c_l\, \delta(f - \frac{l}{N T} - \frac{p}{T})
                       = \frac{1}{T} \sum_{p\in\IZ}  Z_p(f - \frac{p}{T}) \,.

Hence, the DTFT may be interpreted as a periodic continuation of :math:`Z_p(f)`, i.e.,
as the periodic continuation in :math:`f` of the Fourier transform of the
continuous-time periodic signal.

A brief introduction about analyzing the spectral properties of sampled signals can be
found in the :ref:`tutorial_SpectralAnalysis` section.

.. _tutorial_FFT_Caveats:

Caveats
=======
This subsection briefly discusses practical implications of the periodicity properties
of Fourier transforms and provides recommendations for improving computational
efficiency.



Periodicity in the Time Domain
------------------------------
When working with sampled signals, it is typically assumed that the underlying
continuous signals are represented by Fourier series. The circumstance that to
accurately model signals with discontinuities requires many Fourier coefficients is
known as Gibbs phenomenon. Since Fourier series are inherently periodic, representing
non-periodic signals introduces artificial discontinuities. A common technique to make
non-periodic signals "more" periodic is subtracting a best-fit polynomial.
:func:`~scipy.signal.detrend` can be used to subtract a signal's mean (0\ :sup:`th`
order polynomial) or linear trend (1\ :sup:`st` order polynomial). NumPy's
:class:`~numpy.polynomial.polynomial.Polynomial` class provides means to
:meth:`~numpy.polynomial.polynomial.Polynomial.fit` higher order polynomials.

The following example illustrates the effect of removing the linear trend from a noisy
signal on its spectrum: The signal :math:`x(t)` is represented by :math:`N=1000`
samples with a sampling interval of :math:`T=0.1\,`\ ms. The first row depicts
:math:`x(t)` (with linear trend) and :math:`x_d(t)` (with removed linear trend). To
determine the periodic components, a magnitude spectrum of :math:`x(t)` and of
:math:`x_d(t)` is calculated and drawn is second and third row.

.. plot:: tutorial/examples/fft_remove_trend.py
    :include-source: True

The signal's trend manifests as a decreasing slope in the spectrum's low-frequency part.
Removing the linear trend steepens this slope significantly, making the peaks at
:math:`100`, :math:`200`, and :math:`300\,`\ Hz much better distinguishable. These peaks, which
correspond to the sinusoidal components of :math:`x(t)`, would have a magnitude of
:math:`0.5` in a noise- and trend-free signal.

Note that detrending a signal does not necessarily improve its spectrum. I.e., periodic
components whose frequencies are not integer multiples of the frequency resolution
:math:`\Delta f` produce a trend that should not be removed. More information on
calculating and interpreting spectra can be found in the
:ref:`tutorial_SpectralAnalysis` section.


Periodicity in the Frequency Domain
-----------------------------------
The following plot compares :func:`~scipy.fft.fft` and :func:`~scipy.fft.rfft`
representations of two signals made up of an odd number (:math:`N=5`, left column) and
an even number of samples (:math:`N=4`, right column). The real-valued input signals
were chosen in such a way so that their Fourier transforms are real-valued, i.e.,
``x[k] = x[N-k]``. The first row shows the DFTs of Eq. :math:numref:`eq_FFT_DFT` with
its non-negative frequencies. The gray stems represent the periodic continuation to the
left and to the right. The center row shows the :func:`~scipy.fft.fft` with the shifted
frequencies produced by :func:`~scipy.fft.fftfreq`, where the two highest frequencies
are shifted to the negative side. The last row depicts the one-sided
:func:`~scipy.fft.rfft`, which only produces only the non-negative frequencies of the
(two-sided) :func:`~scipy.fft.fft`.

.. plot:: tutorial/examples/fft_compare_DFT_fft_rfft.py
    :include-source: True

This plot illustrates the following two aspects:

* For an odd number of samples (here: :math:`N=5`), :func:`~scipy.fft.fftfreq` produces
  an equal number of positive and negative frequencies, i.e. ``[-2, -1, 0, 1, 2]``. For
  an even number of samples (here: :math:`N=4`), it would be valid to have either the
  frequencies,  ``[-2, -1, 0, 1]`` or the frequencies ``[-1, 0, 1, 2]`` due to the FFT
  values at :math:`f = \pm2\,`\ Hz being identical. By convention,
  :func:`~scipy.fft.fftfreq` places the Nyquist frequency
  :math:`f_\text{Ny} = \Delta f\, N/2 = 2\,`\ Hz on the negative side.
* Here, the :func:`~scipy.fft.rfft` values for input lengths of :math:`N=4,5` are equal
  in spite of the input signals being different. This illustrates that the number
  of signal samples `N` need to be passed to :func:`~scipy.fft.irfft` to be able
  reconstruct the original signal, e.g.:

  >>> import numpy as np
  >>> from scipy.fft import fft, irfft, rfft
  ...
  >>> x = np.array([4, 1, 0, 1])  # n = 4
  >>> X = rfft(x)
  >>> y0, y1 = irfft(X, n=4), irfft(X, n=5)
  >>> np.allclose(y0, x)
  True
  >>> len(y1) == len(x)  # signals differ
  False
  >>> np.allclose(rfft(y1), rfft(x))  # one-sided FFTs are equal
  True


Speed Considerations
--------------------
The efficiency of the FFT algorithm depends on how well the number of input samples
:math:`N` can be factored into small prime factors. I.e., it is slowest if :math:`N` is
prime and fastest if :math:`N` is a power of 2. The :mod:`~scipy.fft` module provides
the functions :func:`~scipy.fft.prev_fast_len` and :func:`~scipy.fft.next_fast_len` for
determining adequate input lengths for fast execution.

Furthermore, the :mod:`~scipy.fft` module provides also the context manager
:func:`~scipy.fft.set_workers` to set the maximum number of allowed threads when
calculating the FFT. Note that some third-party FFT libraries provide SciPy backends,
which can be activated by utilizing :func:`~scipy.fft.set_backend`.


2- and N-D discrete Fourier transforms
======================================

The functions :func:`fft2` and :func:`ifft2` provide 2-D FFT and
IFFT, respectively. Similarly, :func:`fftn` and :func:`ifftn` provide
N-D FFT, and IFFT, respectively.

For real-input signals, similarly to :func:`rfft`, we have the functions
:func:`rfft2` and :func:`irfft2` for 2-D real transforms;
:func:`rfftn` and :func:`irfftn` for N-D real transforms.

The example below demonstrates a 2-D IFFT and plots the resulting
(2-D) time-domain signals.

.. plot::
    :alt: "This code generates six heatmaps arranged in a 2x3 grid. The top row shows mostly blank canvases with the exception of two tiny red peaks on each image. The bottom row shows the real-part of the inverse FFT of each image above it. The first column has two dots arranged horizontally in the top image and in the bottom image a smooth grayscale plot of 5 black vertical stripes representing the 2-D time domain signal. The second column has two dots arranged vertically in the top image and in the bottom image a smooth grayscale plot of 5 horizontal black stripes representing the 2-D time domain signal. In the last column the top image has two dots diagonally located; the corresponding image below has perhaps 20 black stripes at a 60 degree angle."

    >>> from scipy.fft import ifftn
    >>> import matplotlib.pyplot as plt
    >>> import matplotlib.cm as cm
    >>> import numpy as np
    >>> N = 30
    >>> f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    >>> xf = np.zeros((N,N))
    >>> xf[0, 5] = 1
    >>> xf[0, N-5] = 1
    >>> Z = ifftn(xf)
    >>> ax1.imshow(xf, cmap=cm.Reds)
    >>> ax4.imshow(np.real(Z), cmap=cm.gray)
    >>> xf = np.zeros((N, N))
    >>> xf[5, 0] = 1
    >>> xf[N-5, 0] = 1
    >>> Z = ifftn(xf)
    >>> ax2.imshow(xf, cmap=cm.Reds)
    >>> ax5.imshow(np.real(Z), cmap=cm.gray)
    >>> xf = np.zeros((N, N))
    >>> xf[5, 10] = 1
    >>> xf[N-5, N-10] = 1
    >>> Z = ifftn(xf)
    >>> ax3.imshow(xf, cmap=cm.Reds)
    >>> ax6.imshow(np.real(Z), cmap=cm.gray)
    >>> plt.show()


Discrete Cosine Transforms
==========================

SciPy provides a DCT with the function :func:`dct` and a corresponding IDCT
with the function :func:`idct`. There are 8 types of the DCT [#WIKI_DCT]_ [#Mak]_;
however, only the first 4 types are implemented in SciPy. "The" DCT generally
refers to DCT type 2, and "the" Inverse DCT generally refers to DCT type 3. In
addition, the DCT coefficients can be normalized differently (for most types,
scipy provides ``None`` and ``ortho``). Two parameters of the dct/idct
function calls allow setting the DCT type and coefficient normalization.

For a single dimension array `x`, ``dct(x, norm='ortho')`` is equal to
MATLAB's "dct(x)".


Type I DCT
----------

SciPy uses the following definition of the unnormalized DCT-I
(``norm=None``):

.. math::

    y[k] = x_0 + (-1)^k x_{N-1} + 2\sum_{n=1}^{N-2} x[n]
    \cos\left(\frac{\pi nk}{N-1}\right),
    \qquad 0 \le k < N.

Note that the DCT-I is only supported for input size > 1.

Type II DCT
-----------

SciPy uses the following definition of the unnormalized DCT-II
(``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n] \cos \left({\pi(2n+1)k \over 2N} \right)
    \qquad 0 \le k < N.

In case of the normalized DCT (``norm='ortho'``), the DCT coefficients
:math:`y[k]` are multiplied by a scaling factor `f`:

.. math::

    f = \begin{cases} \sqrt{1/(4N)}, & \text{if $k = 0$} \\    \sqrt{1/(2N)},
    & \text{otherwise} \end{cases} \, .

In this case, the DCT "base functions" :math:`\phi_k[n] = 2 f \cos
\left({\pi(2n+1)k \over 2N} \right)` become orthonormal:

.. math::

   \sum_{n=0}^{N-1} \phi_k[n] \phi_l[n] = \delta_{lk}.


Type III DCT
------------

SciPy uses the following definition of the unnormalized DCT-III
(``norm=None``):

.. math::

    y[k] = x_0 + 2 \sum_{n=1}^{N-1} x[n] \cos\left({\pi n(2k+1) \over 2N}\right)
    \qquad 0 \le k < N,

or, for ``norm='ortho'``:

.. math::

    y[k] = {x_0\over\sqrt{N}} + {2\over\sqrt{N}} \sum_{n=1}^{N-1} x[n]
    \cos\left({\pi n(2k+1) \over 2N}\right) \qquad 0 \le k < N.


Type IV DCT
-----------

SciPy uses the following definition of the unnormalized DCT-IV
(``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n] \cos\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N,

or, for ``norm='ortho'``:

.. math::

    y[k] = \sqrt{2\over N}\sum_{n=0}^{N-1} x[n] \cos\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N


DCT and IDCT
------------


The (unnormalized) DCT-III is the inverse of the (unnormalized) DCT-II, up to a
factor of ``2N``. The orthonormalized DCT-III is exactly the inverse of the
orthonormalized DCT- II. The function :func:`idct` performs the mappings between
the DCT and IDCT types, as well as the correct normalization.

The following example shows the relation between DCT and IDCT for different
types and normalizations.

>>> from scipy.fft import dct, idct
>>> x = np.array([1.0, 2.0, 1.0, -1.0, 1.5])

The DCT-II and DCT-III are each other's inverses, so for an orthonormal transform
we return back to the original signal.

>>> dct(dct(x, type=2, norm='ortho'), type=3, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Doing the same under default normalization, however, we pick up an extra scaling
factor of :math:`2N=10` since the forward transform is unnormalized.

>>> dct(dct(x, type=2), type=3)
array([ 10.,  20.,  10., -10.,  15.])

For this reason, we should use the function `idct` using the same type for both,
giving a correctly normalized result.

>>> # Normalized inverse: no scaling factor
>>> idct(dct(x, type=2), type=2)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Analogous results can be seen for the DCT-I, which is its own inverse up to a
factor of :math:`2(N-1)`.

>>> dct(dct(x, type=1, norm='ortho'), type=1, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>> # Unnormalized round-trip via DCT-I: scaling factor 2*(N-1) = 8
>>> dct(dct(x, type=1), type=1)
array([ 8. ,  16.,  8. , -8. ,  12.])
>>> # Normalized inverse: no scaling factor
>>> idct(dct(x, type=1), type=1)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

And for the DCT-IV, which is also its own inverse up to a factor of :math:`2N`.

>>> dct(dct(x, type=4, norm='ortho'), type=4, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>> # Unnormalized round-trip via DCT-IV: scaling factor 2*N = 10
>>> dct(dct(x, type=4), type=4)
array([ 10.,  20.,  10., -10.,  15.])
>>> # Normalized inverse: no scaling factor
>>> idct(dct(x, type=4), type=4)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Example
-------

The DCT exhibits the "energy compaction property", meaning that for many
signals only the first few DCT coefficients have significant magnitude.
Zeroing out the other coefficients leads to a small reconstruction error, a
fact which is exploited in lossy signal compression (e.g. JPEG compression).

The example below shows a signal x and two reconstructions (:math:`x_{20}` and
:math:`x_{15}`) from the signal's DCT coefficients. The signal :math:`x_{20}`
is reconstructed from the first 20 DCT coefficients, :math:`x_{15}` is
reconstructed from the first 15 DCT coefficients. It can be seen that the
relative error of using 20 coefficients is still very small (~0.1%), but
provides a five-fold compression rate.


.. plot::
    :alt: "This code generates an X-Y plot showing amplitude on the Y axis and time on the X axis. The first blue trace is the original signal and starts at amplitude 1 and oscillates down to 0 amplitude over the duration of the plot resembling a frequency chirp. The second red trace is the x_20 reconstruction using the DCT and closely follows the original signal in the high amplitude region but it is unclear to the right side of the plot. The third green trace is the x_15 reconstruction using the DCT and is less precise than the x_20 reconstruction but still similar to x."

    >>> from scipy.fft import dct, idct
    >>> import matplotlib.pyplot as plt
    >>> N = 100
    >>> t = np.linspace(0,20,N, endpoint=False)
    >>> x = np.exp(-t/3)*np.cos(2*t)
    >>> y = dct(x, norm='ortho')
    >>> window = np.zeros(N)
    >>> window[:20] = 1
    >>> yr = idct(y*window, norm='ortho')
    >>> sum(abs(x-yr)**2) / sum(abs(x)**2)
    0.0009872817275276098
    >>> plt.plot(t, x, '-bx')
    >>> plt.plot(t, yr, 'ro')
    >>> window = np.zeros(N)
    >>> window[:15] = 1
    >>> yr = idct(y*window, norm='ortho')
    >>> sum(abs(x-yr)**2) / sum(abs(x)**2)
    0.06196643004256714
    >>> plt.plot(t, yr, 'g+')
    >>> plt.legend(['x', '$x_{20}$', '$x_{15}$'])
    >>> plt.grid()
    >>> plt.show()

Discrete Sine Transforms
========================

SciPy provides a DST [#Mak]_ with the function :func:`dst` and a corresponding IDST
with the function :func:`idst`.

There are, theoretically, 8 types of the DST for different combinations of
even/odd boundary conditions and boundary offsets [#WIKI_DST]_, only the first 4
types are implemented in SciPy.

Type I DST
----------

DST-I assumes the input is odd around n=-1 and n=N. SciPy uses the following
definition of the unnormalized DST-I (``norm=None``):

.. math::

    y[k] = 2\sum_{n=0}^{N-1} x[n]  \sin\left( \pi {(n+1) (k+1)}\over{N+1}
    \right), \qquad 0 \le k < N.

Note also that the DST-I is only supported for input size > 1. The
(unnormalized) DST-I is its own inverse, up to a factor of ``2(N+1)``.

Type II DST
-----------

DST-II assumes the input is odd around n=-1/2 and even around n=N. SciPy uses
the following definition of the unnormalized DST-II (``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n]  \sin\left( {\pi (n+1/2)(k+1)} \over N
    \right), \qquad 0 \le k < N.

Type III DST
------------

DST-III assumes the input is odd around n=-1 and even around n=N-1. SciPy uses
the following definition of the unnormalized DST-III (``norm=None``):

.. math::

    y[k] = (-1)^k x[N-1] + 2 \sum_{n=0}^{N-2} x[n] \sin \left( {\pi
    (n+1)(k+1/2)} \over N \right), \qquad 0 \le k < N.

Type IV DST
-----------

SciPy uses the following definition of the unnormalized DST-IV
(``norm=None``):

.. math::

    y[k] = 2 \sum_{n=0}^{N-1} x[n] \sin\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N,

or, for ``norm='ortho'``:

.. math::

    y[k] = \sqrt{2\over N}\sum_{n=0}^{N-1} x[n] \sin\left({\pi (2n+1)(2k+1) \over 4N}\right)
    \qquad 0 \le k < N,


DST and IDST
------------


The following example shows the relation between DST and IDST for
different types and normalizations.

>>> from scipy.fft import dst, idst
>>> x = np.array([1.0, 2.0, 1.0, -1.0, 1.5])

The DST-II and DST-III are each other's inverses, so for an orthonormal transform
we return back to the original signal.

>>> dst(dst(x, type=2, norm='ortho'), type=3, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Doing the same under default normalization, however, we pick up an extra scaling
factor of :math:`2N=10` since the forward transform is unnormalized.

>>> dst(dst(x, type=2), type=3)
array([ 10.,  20.,  10., -10.,  15.])

For this reason, we should use the function `idst` using the same type for both,
giving a correctly normalized result.

>>> idst(dst(x, type=2), type=2)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

Analogous results can be seen for the DST-I, which is its own inverse up to a
factor of :math:`2(N-1)`.

>>> dst(dst(x, type=1, norm='ortho'), type=1, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>>  # scaling factor 2*(N+1) = 12
>>> dst(dst(x, type=1), type=1)
array([ 12.,  24.,  12., -12.,  18.])
>>>  # no scaling factor
>>> idst(dst(x, type=1), type=1)
array([ 1. ,  2. ,  1. , -1. ,  1.5])

And for the DST-IV, which is also its own inverse up to a factor of :math:`2N`.

>>> dst(dst(x, type=4, norm='ortho'), type=4, norm='ortho')
array([ 1. ,  2. ,  1. , -1. ,  1.5])
>>>  # scaling factor 2*N = 10
>>> dst(dst(x, type=4), type=4)
array([ 10.,  20.,  10., -10.,  15.])
>>>  # no scaling factor
>>> idst(dst(x, type=4), type=4)
array([ 1. ,  2. ,  1. , -1. ,  1.5])


Fast Hankel Transform
=====================

SciPy provides the functions ``fht`` and ``ifht`` to perform the Fast
Hankel Transform (FHT) and its inverse (IFHT) on logarithmically-spaced input
arrays.

The FHT is the discretized version of the continuous Hankel transform defined
by [#Ham00]_

.. math::

    A(k) = \int_{0}^{\infty} \! a(r) \, J_{\mu}(kr) \, k \, dr \;,

with :math:`J_{\mu}` the Bessel function of order :math:`\mu`. Under a change
of variables :math:`r \to \log r`, :math:`k \to \log k`, this becomes

.. math::

    A(e^{\log k})
    = \int_{0}^{\infty} \! a(e^{\log r}) \, J_{\mu}(e^{\log k + \log r})
                                        \, e^{\log k + \log r} \, d{\log r}

which is a convolution in logarithmic space. The FHT algorithm uses the FFT
to perform this convolution on discrete input data.

Care must be taken to minimise numerical ringing due to the circular nature
of FFT convolution. To ensure that the low-ringing condition [#Ham00]_ holds,
the output array can be slightly shifted by an offset computed using the
``fhtoffset`` function.


References
==========

.. [#CT65] Cooley, James W., and John W. Tukey, 1965, "An algorithm for the
        machine calculation of complex Fourier series," *Math. Comput.*
        19: 297-301.

.. [#NR07] Press, W., Teukolsky, S., Vetterline, W.T., and Flannery, B.P.,
        2007, *Numerical Recipes: The Art of Scientific Computing*, ch.
        12-13.  Cambridge Univ. Press, Cambridge, UK.

.. [#WIKI_DFT] "Discrete Fourier transform", Wikipedia,
        https://en.wikipedia.org/wiki/Discrete_Fourier_transform

.. [#WIKI_DTFT] "Discrete-time Fourier transform", Wikipedia,
        https://en.wikipedia.org/wiki/Discrete-time_Fourier_transform

.. [#WIKI_DCT] "Discrete cosine transform", Wikipedia,
        https://en.wikipedia.org/wiki/Discrete_cosine_transform

.. [#Mak] J. Makhoul, 1980, 'A Fast Cosine Transform in One and Two Dimensions',
        `IEEE Transactions on acoustics, speech and signal processing`
        vol. 28(1), pp. 27-34, :doi:`10.1109/TASSP.1980.1163351`

.. [#WIKI_DST] "Discrete sine transform", Wikipedia,
        https://en.wikipedia.org/wiki/Discrete_sine_transform

.. [#Ham00] A. J. S. Hamilton, 2000, "Uncorrelated modes of the non-linear power
        spectrum", *MNRAS*, 312, 257. :doi:`10.1046/j.1365-8711.2000.03071.x`





