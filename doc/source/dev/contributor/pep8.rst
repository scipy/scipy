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

-  Note, however, that SciPy's linting configuration may not match
   that of your IDE exactly. See below on how to run the official
   checks.

-  It is recommended to leave existing style issues alone
   unless they exist in lines of code you are already modifying.
   This practice ensures that the codebase is gradually cleaned up
   without dedicating precious review time to style-only cleanups.

-  Before sending a Pull Request, run the linter on changes made in
   your feature branch. The checks will also be made during
   continuous integration, but it's quicker to catch them early.

   The easiest way to do so is to install our pre-commit hook (once)::

     cp tools/pre-commit-hook.py .git/hooks/pre-commit

   This will run linting checks before each commit is made.

   Alternatively, you can run the check manually from the SciPy root directory::

      python dev.py lint

   You can also run the linter on specific files, using the ``--files`` option::

      python tools/lint.py --files scipy/odr/models.py scipy/ndimage

-  If you have existing code with a lot of PEP8 issues, consider using
   |autopep8|_ to automatically fix them before incorporating the code into
   SciPy.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _enable Real-time code style analysis: https://stackoverflow.com/questions/51463223/how-to-use-pep8-module-using-spyder

.. |autopep8| replace:: ``autopep8``
.. _autopep8: https://pypi.org/project/autopep8/0.8/

.. include:: <isonum.txt>
