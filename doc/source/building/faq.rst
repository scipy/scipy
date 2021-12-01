.. _build-install-faq:

=================
Build/Install FAQ
=================

*How do I set up a development version of SciPy in parallel to a released
version that I use to do my job/research?*

If you use the ``conda`` package manager, this is covered in the
:ref:`quickstart-mac`.

Another simple way to achieve this is to install the released version in
site-packages, by using a binary installer or pip, for example, and set
up the development version in a virtualenv. First, install
`virtualenv`_ (optionally, use `virtualenvwrapper`_), then create your
virtualenv (named scipy-dev here) with::

    $ virtualenv scipy-dev

Now, whenever you want to switch to the virtual environment, you can use the
command ``source scipy-dev/bin/activate``, and ``deactivate`` to exit from the
virtual environment and back to your previous shell. With scipy-dev
activated, first install Scipy's dependencies::

    $ pip install numpy pytest cython pybind11

After that, you can install a development version of Scipy, for example, via::

    $ python setup.py install

The installation goes to the virtual environment.


*How do I set up an in-place build for development?*

For development, you can set up an in-place build so that changes made to
``.py`` files have effect without rebuild. First, run::

    $ python setup.py build_ext -i

Then you need to point your PYTHONPATH environment variable to this directory.
Some IDEs (`Spyder`_, for example) have utilities to manage PYTHONPATH. On Linux
and OSX, you can run the command::

    $ export PYTHONPATH=$PWD

and on Windows::

    $ set PYTHONPATH=/path/to/scipy

Now, editing a Python source file in SciPy allows you to immediately
test and use your changes (in ``.py`` files), by simply restarting the
interpreter.

If you use the macOS, please find more information in the
:ref:`quickstart-mac`.

.. _virtualenv: https://virtualenv.pypa.io/

.. _virtualenvwrapper: https://bitbucket.org/dhellmann/virtualenvwrapper/

.. _Spyder: https://www.spyder-ide.org/
