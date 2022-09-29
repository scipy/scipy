:orphan:

.. _pep8-scipy:

==============
PEP8 and SciPy
==============

All SciPy Python code should adhere to `PEP8`_ style guidelines. It’s so
important that some continuous integration tests on GitHub will fail due
to certain PEP8 violations. Here are a few tips for ensuring PEP8
compliance before pushing your code:

-  Many integrated development environments (IDEs) have options that
   automatically check for PEP8 compliance. In Spyder, for example,
   `enable Real-time code style analysis`_ in Tools |rarr| Preferences |rarr|
   Editor |rarr| Code Introspection/Analysis and "Automatically remove
   trailing spaces when saving files" in in Tools |rarr| Preferences |rarr|
   Editor |rarr| Advanced Settings. This can help you fix PEP8 issues as you
   write your code.

-  SciPy makes use of special configuration files for linting with the
   |flake8|_ tool.

-  It is typically recommended to leave any existing style issues alone
   unless they are part of the code you're already modifying.
   This practice ensures that the codebase is gradually cleaned up
   without dedicating precious review time to style-only cleanups.
   Before sending a Pull Request, we suggest running the lint tests only
   for the changes you've made in your feature branch. This will mimic
   the continuous integration linting checks setup on GitHub.
   After installing ``flake8``, you can run the following check locally
   in the SciPy root directory to ensure your Pull Request doesn't
   break the Continuous Integration linting tests.

   ::

      python tools/lint_diff.py

   If you want to run the diff based lint tests only for specific files
   or directories, please consider using the ``--files`` option.

   ::

      python tools/lint_diff.py --files scipy/odr/models.py scipy/ndimage

-  If you have existing code with a lot of PEP8 issues, consider using
   |autopep8|_ to automatically fix them before incorporating the code into
   SciPy.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _enable Real-time code style analysis: https://stackoverflow.com/questions/51463223/how-to-use-pep8-module-using-spyder

.. |flake8| replace:: ``flake8``
.. _flake8: http://flake8.pycqa.org/en/latest/

.. |autopep8| replace:: ``autopep8``
.. _autopep8: https://pypi.org/project/autopep8/0.8/

.. include:: <isonum.txt>
