"""
Checks the signature file for LAPACK entries in the actual
local library.

It does so by reading the cython_lapack_signatures.txt
file and tries to link all signatures to check their 
existance.
"""
from __future__ import division, print_function, absolute_import

from _library_check import get_linking_signatures


def run():
    from numpy.distutils.system_info import get_info, NotFoundError

    lapack_opt = get_info('lapack_opt')
    print('The LAPACK library option:')
    print(lapack_opt)

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    sigs = get_linking_signatures('cython_lapack_signatures.txt',lapack_opt)

    with open('cython_lapack_signatures_actual.txt','w') as f:
        f.writelines(sigs)

    
if __name__ == '__main__':
    run()
