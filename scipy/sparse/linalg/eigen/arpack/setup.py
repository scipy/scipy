#!/usr/bin/env python
import re
from os.path import join

def needs_veclib_wrapper(info):
    """Returns true if needs special veclib wrapper."""
    import re
    r_accel = re.compile("Accelerate")
    r_vec = re.compile("vecLib")
    res = False
    try:
        tmpstr = info['extra_link_args']
        for i in tmpstr:
            if r_accel.search(i) or r_vec.search(i):
                res = True
    except KeyError:
        pass

    return res

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, NotFoundError
    from numpy.distutils.misc_util import Configuration

    config = Configuration('arpack',parent_package,top_path)

    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError('no lapack/blas resources found')

    config = Configuration('arpack', parent_package, top_path)

    arpack_sources=[join('ARPACK','SRC', '*.f')]
    arpack_sources.extend([join('ARPACK','UTIL', '*.f')])
    arpack_sources.extend([join('ARPACK','LAPACK', '*.f')])

    if needs_veclib_wrapper(lapack_opt):
        arpack_sources += [join('ARPACK', 'FWRAPPERS', 'veclib_cabi_f.f'),
                           join('ARPACK', 'FWRAPPERS', 'veclib_cabi_c.c')]
    else:
        arpack_sources += [join('ARPACK', 'FWRAPPERS', 'dummy.f')]

    config.add_library('arpack_scipy', sources=arpack_sources,
                       include_dirs=[join('ARPACK', 'SRC')],
                       depends = [join('ARPACK', 'FWRAPPERS',
                                       'veclib_cabi_f.f'),
                                  join('ARPACK', 'FWRAPPERS',
                                       'veclib_cabi_c.c'),
                                  join('ARPACK', 'FWRAPPERS',
                                        'dummy.f')])


    config.add_extension('_arpack',
                         sources='arpack.pyf.src',
                         libraries=['arpack_scipy'],
                         extra_info = lapack_opt
                        )

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
