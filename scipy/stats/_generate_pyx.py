import pathlib
from shutil import copyfile


def isNPY_OLD():
    '''
    A new random C API was added in 1.18 and became stable in 1.19.
    Prefer the new random C API when building with recent numpy.
    '''
    import numpy as np
    ver = tuple(int(num) for num in np.__version__.split('.')[:2])
    return ver < (1, 19)


def make_biasedurn():
    '''Substitute True/False values for NPY_OLD Cython build variable.'''
    biasedurn_base = (pathlib.Path(__file__).parent / 'biasedurn').absolute()
    with open(biasedurn_base.with_suffix('.pyx.templ'), 'r') as src:
        contents = src.read()
    with open(biasedurn_base.with_suffix('.pyx'), 'w') as dest:
        dest.write(contents.format(NPY_OLD=str(bool(isNPY_OLD()))))


def make_boost():
    # create target directory
    (pathlib.Path(__file__).parent / '_boost/src').mkdir(exist_ok=True,
                                                         parents=True)
    src_dir = pathlib.Path(__file__).parent / '_boost/src'

    # copy contents of include into directory to satisfy Cython
    # PXD include conditions
    inc_dir = pathlib.Path(__file__).parent / '_boost/include'
    src = 'templated_pyufunc.pxd'
    copyfile(inc_dir / src, src_dir / src)

    # generate the PXD and PYX wrappers
    from _boost.include.gen_func_defs_pxd import (  # type: ignore
        _gen_func_defs_pxd)
    from _boost.include.code_gen import _ufunc_gen  # type: ignore
    from _boost._info import (  # type: ignore
        _x_funcs, _no_x_funcs, _klass_mapper)
    _gen_func_defs_pxd(
        f'{src_dir}/func_defs.pxd',
        x_funcs=_x_funcs,
        no_x_funcs=_no_x_funcs)
    for b, s in _klass_mapper.items():
        _ufunc_gen(
            scipy_dist=s.scipy_name,
            types=['NPY_FLOAT', 'NPY_DOUBLE', 'NPY_LONGDOUBLE'],
            ctor_args=s.ctor_args,
            filename=f'{src_dir}/{s.scipy_name}_ufunc.pyx',
            boost_dist=f'{b}_distribution',
            x_funcs=_x_funcs,
            no_x_funcs=_no_x_funcs,
        )


if __name__ == '__main__':
    make_biasedurn()
    make_boost()
