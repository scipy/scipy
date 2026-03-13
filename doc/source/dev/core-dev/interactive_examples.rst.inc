.. _interactive-docs:

Interactive examples in docstrings
==================================

Examples within docstrings can be made interactive using the
`jupyterlite-sphinx extension <https://jupyterlite-sphinx.readthedocs.io/en/stable/>`__,
which implements the ``.. try_examples::`` directive. When used in the
"Examples" section of a docstring, this directive will create a button that
allows users to open the example in an interactive JupyterLite window, or in a
new tab as a Jupyter notebook. For SciPy, this directive is automatically added
to all docstrings with an "Examples" section through the
``global_enable_try_examples`` configuration option in ``doc/source/conf.py``.

To hide the buttons that turn docstring examples into embedded notebooks, after
building the html documentation, edit the ``try_examples.json`` runtime
configuration file after ``scipy/doc/build/html/`` and add the url paths of the
files to be ignored to the ``ignore_patterns`` list. This list expects elements
to be JavaScript Regex patterns. The initial version of this file in a clean
documentation build is

.. code-block:: json

   {
       "global_min_height": "400px",
       "ignore_patterns": [".*"]
   }

The buttons for converting examples sections into embedded interactive
notebooks are hidden by default on clean docs builds.

See the documentation for the ``jupyterlite-sphinx``
`TryExamples directive <https://jupyterlite-sphinx.readthedocs.io/en/stable/directives/try_examples.html#>`__
for more information.

