# prevent interference with KeyboardInterrupt on Windows
# due to Fortran libraries
# See stackoverflow for explanation:
# http://stackoverflow.com/questions/15457786/ctrl-c-crashes-python-after-importing-scipy-stats
import imp
import ctypes
import os

INSTALL = False

dirname = os.path.dirname(__file__)
config_file = os.path.join(dirname, '__config__.py')

if os.path.exists(config_file):
    with open(config_file, 'rb') as fid:
        text = fid.read()
    if 'mkl_blas' in text:
        INSTALL = True


def handler(sig):
    try:
        import _thread
    except ImportError:
        import thread as _thread
    _thread.interrupt_main()
    return 1

# load numpy  math and fortran libraries (but do not import numpy)
basepath = imp.find_module('numpy')[1]
ctypes.CDLL(os.path.join(basepath, 'core', 'libmmd.dll'))
ctypes.CDLL(os.path.join(basepath, 'core', 'libifcoremd.dll'))
# install handler
routine = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_uint)(handler)

if INSTALL:
    ctypes.windll.kernel32.SetConsoleCtrlHandler(routine, 1)
