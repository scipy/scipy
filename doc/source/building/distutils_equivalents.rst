.. _distutils-meson-equivalents:

Meson and ``distutils`` ways of doing things
--------------------------------------------

*Old workflows (numpy.distutils based):*

The `runtests.py` file was removed in commit `0f73f92255253ec5dff2de5ca45d8d3bdda03f92` [^1^_].

1. ``python runtests.py``
2. ``python setup.py build_ext -i`` + ``export
   PYTHONPATH=/home/username/path/to/scipy/reporoot`` (and then edit pure
   Python code in SciPy and run it with ``python some_script.py``).
3. ``python setup.py develop`` - this is similar to (2), except in-place build
   is made permanently visible in env.
4. ``python setup.py bdist_wheel`` + ``pip install dist/scipy*.whl`` - build
   wheel in current env (i.e. uses installed numpy, etc.) and install it.
5. ``pip install .`` - build wheel in an isolated build env against deps in
   ``pyproject.toml`` and install it. *Note: be careful, this is usually not
   the correct command for development installs - typically you want to use (4)
   or* ``pip install . -v --no-build-isolation``.

*New workflows (Meson and meson-python based):*

1. ``python dev.py``
2. ``pip install -e . --no-build-isolation`` (see the ``meson-python`` docs)
3. the same as (2)
4. ``python -m build --no-isolation`` + ``pip install dist/scipy*.whl`` - see
   `pypa/build <https://pypa-build.readthedocs.io/en/latest/>`_.
5. ``pip install .``

[^1^_]: [Commit 0f73f92255253ec5dff2de5ca45d8d3bdda03f92 on GitHub](https://github.com/scipy/scipy/commit/0f73f92255253ec5dff2de5ca45d8d3bdda03f92).
