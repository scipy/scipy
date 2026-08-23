Building SciPy for Contributors
```````````````````````````````

.. note::

    This page describes how to build and develop SciPy using its Pixi workspace,
    which provides the most accessible developer experience.

    For developers who do not want to use Pixi or who require customised environments,
    see :doc:`the advanced developer building guide <contributor-advanced>`.
    
Development of SciPy is made easy with `Pixi <https://pixi.prefix.dev>`__.
First, `clone the SciPy repository <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`__::

      git clone https://github.com/scipy/scipy.git
      cd scipy
      git submodule update --init

and `install Pixi <https://pixi.prefix.dev/latest/installation>`__.

SciPy can then be built with the single command::

    pixi run build

All other common development tasks are also available via ``pixi run``:

.. code-block:: console

    pixi run test       # run the tests
    pixi run open-docs  # build and preview the docs
    pixi run lint       # run main lint checks
    pixi run ipython    # spawn an ipython prompt with SciPy installed
    pixi run smoke-docs # run the doctests
    pixi run test-cpu   # run the tests with all cpu array backends
    pixi run bench      # run the benchmarks

.. tip::

    Run ``pixi task list`` for a full list of available tasks.

.. tip::

    Run ``pixi info`` for a full list of environments and their tasks.
