
# This is a hack to make sure we can import the shared libraries;
# All installed SOs should all be siblings to __init__.py.  The rpath
# of the SOs are set to be relative (=./) since there's no way to
# pull out the installed directory from distutils or setuptools
# before the installation happens
import pathlib
from os import chdir
chdir(str(pathlib.Path(__file__).parent))

# Wrappers
from pyHiGHS.highs_wrapper import highs_wrapper

# Constants
from pyHiGHS.highs_wrapper import (
    CONST_I_INF,
    CONST_INF,
)
