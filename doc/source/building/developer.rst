Building SciPy from source for SciPy development
````````````````````````````````````````````````

.. note::

    This page describes how to build and develop SciPy using its Pixi workspace,
    which provides the most accessible developer experience.

    For developers who do not want to use Pixi or who require customised environments,
    see `the advanced developer building guide <developer-advanced>`__.

Development of SciPy is made easy with [Pixi](https://pixi.prefix.dev/).
First, `clone the SciPy repository <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`__::

      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init

and `install Pixi <https://pixi.prefix.dev/latest/installation>`__.

SciPy can then be built with the single command::

    pixi run build

All other common development tasks are also available via ``pixi run``::

    # TODO: use actual SciPy tasks.
    pixi run tests      # run the tests
    pixi run open-docs  # build and preview the docs
    pixi run lint       # run the full lint suite
    pixi run ipython    # spawn an ipython prompt with array-api-extra installed
    pixi run hooks      # install pre-commit hooks

.. tip::

    Run `pixi task list` for a full list of available tasks.

.. tip::

    Run `pixi info` for a full list of environments and their tasks.
