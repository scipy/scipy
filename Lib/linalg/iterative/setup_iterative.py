#!/usr/bin/env python
#
# Sparse Iterative Package for SciPy
# Hongze Liu, Travis E. Oliphant,
# Brigham Young University
# 2004
#

import os
from glob import glob
from scipy_distutils.core import Extension
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
from scipy_distutils.system_info import get_info,dict_append,NotFoundError


def configuration(parent_package='',parent_path=None):

    package = 'iterative'
    
    config = default_config_dict(package,parent_package)
    lapack_opt = get_info('lapack_opt')


#    sources = ['getbreak.f']
    BiCG = ['BiCG.f.src','BiCGREVCOM.f.src']
    BiCGSTAB = ['BiCGSTAB.f.src','BiCGSTABREVCOM.f.src']
    CG = ['CG.f.src','CGREVCOM.f.src']
    CGS = ['CGS.f.src','CGSREVCOM.f.src']
    Cheby = ['Cheby.f.src','ChebyREVCOM.f.src']
#    GMRES = ['GMRES.f.src','GMRESREVCOM.f.src']
    GMRES = ['GMRES.f.src','GMRESREVCOM.f']
#    Jacobi = ['Jacobi.f.src','JacobiREVCOM.f.src']
    QMR = ['QMR.f.src','QMRREVCOM.f.src']
    SOR = ['SOR.f.src','SORREVCOM.f.src']
    Util = ['STOPTEST2.f.src','getbreak.f.src']
    sources = Util + BiCG + BiCGSTAB + CG + CGS + QMR + ['iterative.pyf.src']
    ext_args = {
        'name': dot_join(parent_package,package,'iterative'),
        'sources': sources,
        }
    dict_append(ext_args, **lapack_opt)
    ext = Extension(**ext_args)
    config['ext_modules'].append(ext)

    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
