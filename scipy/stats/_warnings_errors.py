# This module exists only to allow Sphinx to generate docs
# for the result objects returned by some functions in stats
# _without_ adding them to the main stats documentation page.

"""
Errors
------

.. currentmodule:: scipy.stats._warnings_errors

.. autosummary::
   :toctree: generated/

   FitError

Warnings
--------

.. currentmodule:: scipy.stats._warnings_errors

.. autosummary::
   :toctree: generated/

   DegenerateDataWarning

"""

# Warnings


class DegenerateDataWarning(RuntimeWarning):
    """Warns when data is degenerate and results may not be reliable."""
    def __init__(self, msg=None):
        if msg is None:
            msg = ("Degenerate data encountered; results may not be reliable.")
        self.args = (msg,)


class ConstantInputWarning(DegenerateDataWarning):
    """Warns when all values in data are exactly equal."""
    def __init__(self, msg=None):
        if msg is None:
            msg = ("All values in data are exactly equal; "
                   "results may not be reliable.")
        self.args = (msg,)


class NearConstantInputWarning(DegenerateDataWarning):
    """Warns when all values in data are nearly equal."""
    def __init__(self, msg=None):
        if msg is None:
            msg = ("All values in data are nearly equal; "
                   "results may not be reliable.")
        self.args = (msg,)


# Errors


class FitError(RuntimeError):
    """Represents an error condition when fitting a distribution to data."""
    def __init__(self, msg=None):
        if msg is None:
            msg = ("An error occured when fitting a distribution to data.")
        self.args = (msg,)
