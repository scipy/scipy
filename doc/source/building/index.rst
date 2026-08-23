.. _building-from-source:

Building from source
====================

Building SciPy from source requires setting up system-level dependencies
(compilers, BLAS/LAPACK libraries, etc.) first, and then invoking a build. The
build may be done in order to install SciPy for local usage, develop SciPy
itself, or build redistributable binary packages.
See the guides below for key steps to follow:

.. grid:: 1 1 2 2
    :gutter: 2 2 2 2

    .. grid-item-card::
        :text-align: center

        **Building from source to use SciPy**
        ^^^

        .. button-ref:: user
            :color: secondary
            :click-parent:

            To the user building guide

    .. grid-item-card::
        :text-align: center

        **Building from source to develop SciPy**
        ^^^

        .. button-ref:: developer
            :color: secondary
            :click-parent:

            To the developer building guide

.. toctree::
   :hidden:

   User building guide <user>
   Developer building guide <developer>


Customizing builds
------------------

It may be desired to customize aspects of how the build is done.
See the following pages for guidance on various customization options.

.. toctree::
   :maxdepth: 1

   compilers_and_options
   blas_lapack
   cross_compilation
   redistributable_binaries


Background information
----------------------

See the following pages for background information on how the SciPy build
works, and links to up-to-date guides for generic Python build & packaging
documentation that is relevant.

.. toctree::
    :maxdepth: 1

    understanding_meson
    introspecting_a_build
    distutils_equivalents
