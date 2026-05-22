.. This following environment provides LaTeX macros for math environments.
   To use these macros add the following line to the applicable `.rst` file:

   .. include:: MathMacros.rst

   The macros are defined after this line to the end of the file.
   Currently, this file is only included by the file `fft.rst`.


.. math::

    \newcommand{\IC}{{\mathbb{C}}}  % set of complex numbers
    \newcommand{\IN}{{\mathbb{N}}}  % set of natural numbers
    \newcommand{\IR}{{\mathbb{R}}}  % set of real numbers
    \newcommand{\IZ}{{\mathbb{Z}}}  % set of integers
    \newcommand{\jj}{{\mathbb{j}}}  % imaginary unit
    \newcommand{\e}{\operatorname{e}}  % Euler's number
    \newcommand{\dd}{\operatorname{d}} % infinitesimal operator
    \newcommand{\abs}[1]{\left|#1\right|} % absolute value
    \newcommand{\conj}[1]{\overline{#1}} % complex conjugate
    \newcommand{\conjT}[1]{\overline{#1^T}} % transposed complex conjugate
    \newcommand{\inv}[1]{\left(#1\right)^{\!-1}} % inverse
    \newcommand{\rect}{\operatorname{rect}}  % rect or boxcar function
    \newcommand{\sinc}{\operatorname{sinc}}  % sinc(t) := sin(pi*t) / (pi*t)
    % overwrite macros:
    \renewcommand{\Re}{\operatorname{Re}}  % real part of a complex number
    \renewcommand{\Im}{\operatorname{Im}}  % imaginary part of a complex number
    % Since the physics package is not loaded, we define the macros ourselves:
    \newcommand{\vb}[1]{\mathbf{#1}} % vectors  are bold
    \newcommand{\mat}[1]{\mathrm{#1}} % matrices are upright
    %
    % Abbreviations of matrices and vectors :
    \newcommand{\mF}{\mat{F}}  % matrix with DFT coefficients
    \newcommand{\mI}{\mat{I}}  % unity matrix

