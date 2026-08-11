import os
import pytest
import sys

@pytest.fixture(scope='function')
def minimize_with_debugging():
    # This is a hack to force us to use DEBUGGING for a test, so that we get
    # code coverage. We definitely don't want to do this for the pycutest tests,
    # they are slow enough with USE_NAIVE_MATH set to True.
    os.environ['PRIMA_DEBUGGING'] = "True"
    if 'pyprima' in sys.modules:
        modules_to_delete = [m for m in sys.modules if 'pyprima' in m]
        for m in modules_to_delete:
            del sys.modules[m]
    # We have to reimport minimize here, if the tests try to use the minimize that
    # was imported at the top of the file, it will be the existing object which
    # captured DEBUGGING=False.
    from pyprima import minimize
    yield minimize
    del os.environ['PRIMA_DEBUGGING']