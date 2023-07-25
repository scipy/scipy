Cross compilation
=================

Cross compilation is a complex topic, we only add some hopefully helpful hints
here (for now). As of May 2023, cross-compilation based on ``crossenv`` is
known to work, as used (for example) in conda-forge. Cross-compilation without
``crossenv`` requires some manual overrides. You instruct these overrides by
passing options to ``meson setup`` via `meson-python`_.

.. _meson-python: https://meson-python.readthedocs.io/en/latest/how-to-guides/meson-args.html

All distributions that are known to successfully cross compile SciPy are using
``python -m build`` (``pypa/build``), but using ``pip`` for that should be
possible as well. Here are links to the SciPy's "build recipes" on those
distros:

- `Void Linux <https://github.com/void-linux/void-packages/blob/master/srcpkgs/python3-scipy/template>`_
- `Nix <https://github.com/nixos/nixpkgs/blob/master/pkgs/development/python-modules/scipy/default.nix>`_
- `Conda-forge <https://github.com/conda-forge/scipy-feedstock/blob/main/recipe/build.sh>`_

See also `Meson's documentation on cross compilation
<https://mesonbuild.com/Cross-compilation.html>`__ to learn what options you
may need to pass to Meson to successfully cross compile.

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
