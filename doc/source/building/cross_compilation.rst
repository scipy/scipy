Cross compilation
=================

Cross compilation is a complex topic, we only add some hopefully helpful hints
here (for now). As of May 2023, cross-compilation based on ``crossenv`` is
known to work, as used (for example) in conda-forge. Cross-compilation without
``crossenv`` requires some manual overrides, as done in (for example) Void Linux.

Please see `Meson's documentation on cross compilation
<https://mesonbuild.com/Cross-compilation.html>`__
for details on Meson's support for cross-compilation.

One common hiccup is that ``numpy`` and ``pythran`` require
running Python code in order to obtain their include directories. This tends to
not work well, either accidentally picking up the packages from the build
(native) Python rather than the host (cross) Python or requiring ``crossenv``
or QEMU to run the host Python. To avoid this problem, specify the paths to the
relevant directories in your *cross file*:

.. code:: ini

    [constants]
    sitepkg = '/abspath/to/host-pythons/site-packages/'

    [properties]
    numpy-include-dir = sitepkg + 'numpy/core/include'
    pythran-include-dir = sitepkg + 'pythran'

For more details and the current status around cross compilation, see:

- Tracking issue for SciPy cross-compilation needs and issues:
  `scipy#14812 <https://github.com/scipy/scipy/issues/14812>`__
- The state of cross compilation in Python:
  `pypackaging-native key issue page <https://pypackaging-native.github.io/key-issues/cross_compilation/>`__
