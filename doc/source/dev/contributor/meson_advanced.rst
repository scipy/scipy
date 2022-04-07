.. _meson-advanced:

Build with Meson Optimizations
==============================

Meson provides different optimizations levels while configuring the project. You can see
the available options for optimizations at
`meson documentation <https://mesonbuild.com/Builtin-options.html#core-options>`.

Assuming that you are building from scratch(do ``git clean -xdf`` if needed), you can
configure the build as following::

    meson setup build --optimization g --prefix=$PWD/build-install

Now, you can use the ``dev.py`` interface for further building, installing and testing SciPy::

    python dev.py -s linalg

This will work because after initial configuration, Meson will remember the config options.
