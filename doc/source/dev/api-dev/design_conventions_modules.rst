.. The contents of this file previously resided in file `reference/inex.rst`

.. _design_conventions_modules:

Design Conventions for Modules
==============================
All SciPy modules should follow the following conventions. In the
following, a *SciPy module* is defined as a Python package, say
``yyy``, that is located in the scipy/ directory.

* Ideally, each SciPy module should be as self-contained as possible.
  That is, it should have minimal dependencies on other packages or
  modules. Even dependencies on other SciPy modules should be kept to
  a minimum. A dependency on NumPy is of course assumed.

* Directory ``yyy/`` contains:

  - A file ``meson.build`` with build configuration for the submodule.

  - A directory ``tests/`` that contains files ``test_<name>.py``
    corresponding to modules ``yyy/<name>{.py,.so,/}``.

  - An ``__init__.py`` file, which loads functionality from the other files in the
    directory as needed. Furthermore, it must list all public members, classes and
    attributes into its ´´__all__`` attribute. Following Python's conventions, names
    starting with an underscore character are considered private.

* Private modules should be prefixed with an underscore ``_``,
  for instance ``yyy/_some_module.py``.

* User-visible functions should have good documentation following
  the `NumPy documentation style`_.

* The ``__init__.py`` of the module should contain the main reference documentation in
  its docstring.  This is connected to the Sphinx
  documentation under ``doc/`` via Sphinx's automodule directive.

  The reference documentation should first give a categorized list of
  the contents of the module using ``autosummary::`` directives, and
  after that explain points essential for understanding the use of the
  module.

  Tutorial-style documentation with extensive examples should be
  separate and put under ``doc/source/tutorial/``.

See the existing SciPy submodules for guidance.

.. _NumPy documentation style: https://numpydoc.readthedocs.io/en/latest/format.html
