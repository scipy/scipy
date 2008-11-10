import numpy as np

from scipy.lib.lapack import flapack, clapack

FUNCS_TP = {'ssygv' : np.float32, 
         'dsygv': np.float,
         'ssygvd' : np.float32,
         'dsygvd' : np.float,
         'ssyev' : np.float32, 
         'dsyev': np.float,
         'ssyevr' : np.float32,
         'dsyevr' : np.float,
         'ssyevr' : np.float32,
         'dsyevr' : np.float,
         'sgehrd' : np.float32,
         'dgehrd' : np.float,
         'sgebal' : np.float32, 
         'dgebal': np.float}

# Test FLAPACK if not empty
if hasattr(flapack, 'empty_module'):
    FLAPACK_IS_EMPTY = True
else:
    FLAPACK_IS_EMPTY = False

# Test CLAPACK if not empty and not the same as clapack
if hasattr(clapack, 'empty_module') or (clapack == flapack):
    CLAPACK_IS_EMPTY = True
else:
    CLAPACK_IS_EMPTY = False

if not FLAPACK_IS_EMPTY:
    FUNCS_FLAPACK = {'ssygv' : flapack.ssygv,
                     'dsygv': flapack.dsygv,
                     'ssygvd' : flapack.ssygvd,
                     'dsygvd' : flapack.dsygvd,
                     'ssyev' : flapack.ssyev,
                     'dsyev': flapack.dsyev,
                     'ssyevr' : flapack.ssyevr,
                     'dsyevr' : flapack.dsyevr,
                     'sgehrd' : flapack.sgehrd,
                     'dgehrd' : flapack.dgehrd,
                     'sgebal' : flapack.sgebal,
                     'dgebal': flapack.dgebal}
else:
    FUNCS_FLAPACK = None

if not CLAPACK_IS_EMPTY:
    FUNCS_CLAPACK = {'ssygv' : clapack.ssygv,
                     'dsygv': clapack.dsygv,
                     'ssygvd' : clapack.ssygvd,
                     'dsygvd' : clapack.dsygvd,
                     'ssyev' : clapack.ssyev,
                     'dsyev': clapack.dsyev,
                     'ssyevr' : clapack.ssyevr,
                     'dsyevr' : clapack.dsyevr,
                     'sgehrd' : flapack.sgehrd,
                     'dgehrd' : flapack.dgehrd,
                     'sgebal' : clapack.sgebal,
                     'dgebal': clapack.dgebal}
else:
    FUNCS_CLAPACK = None

PREC = {np.float32: 5, np.float: 12}

