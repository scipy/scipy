Building SciPy from source for SciPy development
````````````````````````````````````````````````

Building SciPy from source for SciPy development can be split into
two main steps: setting up system-level dependencies, and building SciPy itself.

.. _system-level:

.. include:: system_level.rst

.. _the-spin-interface:

Building SciPy itself
^^^^^^^^^^^^^^^^^^^^^

First, clone the SciPy repository::

      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init


.. tip::

    Many of the steps described below can now be accomplished automatically
    with commands which execute tasks in SciPy's Pixi workspace,
    like ``pixi run build``.
    To use this workspace, `install Pixi <https://pixi.sh/latest/installation/>`__
    and execute ``pixi task list`` in a local clone of SciPy's source to see
    the various tasks available.
    
    This removes the need for developers to keep track of development environments
    and installed dependencies, as running a task automatically installs and uses
    a suitable environment.
    A future update to this guide will provide full details on using the Pixi
    workspace for SciPy development.


Then you want to do the following:

1. Create a dedicated development environment (virtual environment or conda
   environment),
2. Install all needed dependencies (*build*, and also *test*, *doc* and
   *optional* dependencies),
3. Build SciPy with our ``spin`` developer interface.

Step (3) is always the same, steps (1) and (2) are different between conda and
virtual environments:

.. tab-set::

  .. tab-item:: Conda env
    :sync: conda

    To create a ``scipy-dev`` development environment with every required and
    optional dependency installed, run::

        conda env create -f environment.yml
        conda activate scipy-dev

  .. tab-item:: Virtual env or system Python
    :sync: pip

    .. note::

       There are many tools to manage virtual environments, like ``venv``,
       ``virtualenv``/``virtualenvwrapper``, ``pyenv``/``pyenv-virtualenv``,
       Poetry, PDM, Hatch, and more. Here we use the basic ``venv`` tool that
       is part of the Python stdlib. You can use any other tool; all we need is
       an activated Python environment.

    Create and activate a virtual environment in a new directory named ``venv`` (
    note that the exact activation command may be different based on your OS and shell
    - see `"How venvs work" <https://docs.python.org/3/library/venv.html#how-venvs-work>`__
    in the ``venv`` docs).

    .. tab-set::

      .. tab-item:: Linux
        :sync: linux

        ::

          python -m venv venv
          source venv/bin/activate

      .. tab-item:: macOS
        :sync: macos

        ::

          python -m venv venv
          source venv/bin/activate

      .. tab-item:: Windows
        :sync: windows

        ::

          python -m venv venv
          venv\Scripts\Activate.ps1

    Then install the Python-level dependencies from PyPI with::

       python -m pip install --group dev

    If you want to install dependencies in a more granular fashion,
    see the ``doc``, ``test``, and other groups under ``[dependency-groups]``
    in ``pyproject.toml``.

To build SciPy in an activated development environment, run::

    spin build

This will install SciPy inside the repository (by default in a
``build-install`` directory). You can then run tests (``spin test``),
drop into IPython (``spin ipython``), or take other development steps
like build the html documentation or running benchmarks. The ``spin``
interface is self-documenting, so please see ``spin --help`` and
``spin <subcommand> --help`` for detailed guidance.


.. admonition:: IDE support & editable installs

    While the ``spin build`` interface is our recommended way of working on SciPy,
    it has one limitation: because of the custom install location, SciPy
    will not be recognized automatically within an
    IDE (e.g., for running a script via a "run" button, or setting breakpoints
    visually). This will work better with an *in-place build* (or "editable
    install").

    Editable installs are supported via ``spin install``.
    When making changes to SciPy code, including to compiled code, there is no
    need to manually rebuild or reinstall. However, should you need to run ``git clean -xdf``,
    which removes the built extension modules, remember to also uninstall SciPy
    with ``pip uninstall scipy``.

    See the meson-python_ documentation on editable installs for more details
    on how things work under the hood.

    Note that editable installations are unsuitable for some forms of development,
    such as working on sections of C/Cython API where tests are disabled for editable
    installations. They also tend to hit weird corner cases more frequently than
    regular installations, and have some known limitations like a lack of support
    for static typing.


.. _meson-python: https://mesonbuild.com/meson-python/


Installing static type stubs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would like to install static type stubs to aid your development of SciPy,
you can include the ``scipy-stubs`` package in your development environment.
It is available on PyPI and conda-forge - see the scipy-stubs_ installation guide.

.. _scipy-stubs: https://github.com/jorenham/scipy-stubs?tab=readme-ov-file#installation
