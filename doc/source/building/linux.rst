.. _build-linux:

Building from source on Linux
=============================

.. tabs::

   .. tab:: Debian/Ubuntu Linux

      .. note::

         These instructions have been tested on Ubuntu Linux 16.04, 18.04, and
         20.04.

      Python should be available in your system via the ``python3`` command. To
      install the remaining system-level dependencies, run::

         sudo apt install -y gcc g++ gfortran libopenblas-dev liblapack-dev pkg-config
         sudo apt install -y python3-pip

      Alternatively, you can do::

         sudo apt build-dep scipy

      This command installs whatever is needed to build SciPy, with the
      advantage that new dependencies or updates to required versions are
      handled by the package managers.

      All further work should proceed in a virtual environment. Popular options
      include the standard library ``venv`` module or a separate ``virtualenv``
      package. See also :ref:`ubuntu-guide`.

      **Custom BLAS**

      To customize which BLAS is used, you can set up a `site.cfg` file. See
      the `site.cfg.example` file in the numpy source for the options you
      can set.

      Note that Debian and Ubuntu package optimized BLAS libraries in an
      exchangeable way. You can install libraries, such as ATLAS or OpenBLAS
      and change the default one used via the alternatives mechanism:

      ::

         $ sudo apt-get install libopenblas-base libatlas3-base
         $ update-alternatives --list libblas.so.3
         /usr/lib/atlas-base/atlas/libblas.so.3
         /usr/lib/libblas/libblas.so.3
         /usr/lib/openblas-base/libopenblas.so.0

         $ sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libopenblas.so.0

      See ``/usr/share/doc/libatlas3-base/README.Debian`` for instructions on
      how to build optimized ATLAS packages for your specific CPU. The
      packaged OpenBLAS chooses the optimal code at runtime so it does not
      need recompiling unless the packaged version does not yet support the
      used CPU.

      You can also use a library you built yourself by preloading it. This does
      not require administrator rights.

      ::

         LD_PRELOAD=/path/to/libatlas.so.3 ./my-application

   .. tab:: Fedora

      To install scipy build requirements, you can do::

         sudo dnf install gcc-gfortran python3-devel openblas-devel lapack-devel

      Alternatively, you can do::

         sudo dnf builddep scipy

      This command installs whatever is needed to build SciPy, with the
      advantage that new dependencies or updates to required versions are
      handled by the package managers.

   .. tab:: CentOS/RHEL

      To install scipy build requirements, you can do::

         sudo yum install gcc-gfortran python3-devel openblas-devel lapack-devel

      Alternatively, you can do::

         sudo yum-builddep scipy

      This command installs whatever is needed to build SciPy, with the
      advantage that new dependencies or updates to required versions are
      handled by the package managers.

   .. tab:: Arch

      To install scipy build requirements, you can do::

         sudo pacman -S gcc-gfortran python-devel openblas-devel lapack-devel
