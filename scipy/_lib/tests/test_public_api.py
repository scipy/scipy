"""
This test script is adopted from:
    https://github.com/numpy/numpy/blob/main/numpy/tests/test_public_api.py
"""

import pkgutil
import types
import importlib
import warnings
from importlib import import_module

import pytest

import scipy

from scipy._lib._public_api import PUBLIC_MODULES
from scipy.conftest import xp_available_backends


def test_dir_testing():
    """Assert that output of dir has only one "testing/tester"
    attribute without duplicate"""
    assert len(dir(scipy)) == len(set(dir(scipy)))


# The PRIVATE_BUT_PRESENT_MODULES list contains modules that lacked underscores
# in their name and hence looked public, but weren't meant to be. All these
# namespace were deprecated in the 1.8.0 release - see "clear split between
# public and private API" in the 1.8.0 release notes.
# These private modules support will be removed in SciPy v2.0.0, as the
# deprecation messages emitted by each of these modules say.
PRIVATE_BUT_PRESENT_MODULES = [
    'scipy.fftpack.convolve',
    'scipy.misc',
    'scipy.misc.common',
    'scipy.misc.doccer',
    'scipy.special.cython_special',
]


def is_unexpected(name):
    """Check if this needs to be considered."""
    if '._' in name or '.tests' in name or '.setup' in name:
        return False

    if name in PUBLIC_MODULES:
        return False

    if name in PRIVATE_BUT_PRESENT_MODULES:
        return False

    return True


SKIP_LIST = [
    'scipy.conftest',
    'scipy.version',
    'scipy.special.libsf_error_state',
    'scipy.integrate.lsoda'
]


# XXX: this test does more than it says on the tin - in using `pkgutil.walk_packages`,
# it will raise if it encounters any exceptions which are not handled by `ignore_errors`
# while attempting to import each discovered package.
# For now, `ignore_errors` only ignores what is necessary, but this could be expanded -
# for example, to all errors from private modules or git subpackages - if desired.
@pytest.mark.thread_unsafe(
    reason=("crashes in pkgutil.walk_packages, see "
            "https://github.com/data-apis/array-api-compat/issues/343"))
def test_all_modules_are_expected():
    """
    Test that we don't add anything that looks like a new public module by
    accident.  Check is based on filenames.
    """

    def ignore_errors(name):
        # if versions of other array libraries are installed which are incompatible
        # with the installed NumPy version, there can be errors on importing
        # `array_api_compat`. This should only raise if SciPy is configured with
        # that library as an available backend.
        backends = {'cupy', 'torch', 'dask.array'}
        for backend in backends:
            path = f'array_api_compat.{backend}'
            if path in name and backend not in xp_available_backends:
                return
        raise

    modnames = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "scipy.misc", DeprecationWarning)
        for _, modname, _ in pkgutil.walk_packages(path=scipy.__path__,
                                                   prefix=scipy.__name__ + '.',
                                                   onerror=ignore_errors):
            if is_unexpected(modname) and modname not in SKIP_LIST:
                # We have a name that is new.  If that's on purpose, add it to
                # PUBLIC_MODULES.  We don't expect to have to add anything to
                # PRIVATE_BUT_PRESENT_MODULES.  Use an underscore in the name!
                modnames.append(modname)

    if modnames:
        raise AssertionError(f'Found unexpected modules: {modnames}')


# Stuff that clearly shouldn't be in the API and is detected by the next test
# below
SKIP_LIST_2 = [
    'scipy.char',
    'scipy.rec',
    'scipy.emath',
    'scipy.math',
    'scipy.random',
    'scipy.ctypeslib',
    'scipy.ma',
    'scipy.integrate.lsoda'
]


def test_all_modules_are_expected_2():
    """
    Method checking all objects. The pkgutil-based method in
    `test_all_modules_are_expected` does not catch imports into a namespace,
    only filenames.
    """

    def find_unexpected_members(mod_name):
        members = []
        module = importlib.import_module(mod_name)
        if hasattr(module, '__all__'):
            objnames = module.__all__
        else:
            objnames = dir(module)

        for objname in objnames:
            if not objname.startswith('_'):
                fullobjname = mod_name + '.' + objname
                if isinstance(getattr(module, objname), types.ModuleType):
                    if is_unexpected(fullobjname) and fullobjname not in SKIP_LIST_2:
                        members.append(fullobjname)

        return members
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",  "scipy.misc", DeprecationWarning)
        unexpected_members = find_unexpected_members("scipy")

    for modname in PUBLIC_MODULES:
        unexpected_members.extend(find_unexpected_members(modname))

    if unexpected_members:
        raise AssertionError("Found unexpected object(s) that look like "
                             f"modules: {unexpected_members}")


def test_api_importable():
    """
    Check that all submodules listed higher up in this file can be imported
    Note that if a PRIVATE_BUT_PRESENT_MODULES entry goes missing, it may
    simply need to be removed from the list (deprecation may or may not be
    needed - apply common sense).
    """
    def check_importable(module_name):
        try:
            importlib.import_module(module_name)
        except (ImportError, AttributeError):
            return False

        return True

    module_names = []
    for module_name in PUBLIC_MODULES:
        if not check_importable(module_name):
            module_names.append(module_name)

    if module_names:
        raise AssertionError("Modules in the public API that cannot be "
                             f"imported: {module_names}")

    with warnings.catch_warnings(record=True):
        warnings.simplefilter('always', category=DeprecationWarning)
        warnings.simplefilter('always', category=ImportWarning)
        for module_name in PRIVATE_BUT_PRESENT_MODULES:
            if not check_importable(module_name):
                module_names.append(module_name)

    if module_names:
        raise AssertionError("Modules that are not really public but looked "
                             "public and can not be imported: "
                             f"{module_names}")
