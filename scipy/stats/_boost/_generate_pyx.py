import pathlib
from shutil import copyfile
from typing import NamedTuple

class _KlassMap(NamedTuple):
    scipy_name: str
    ctor_args: tuple


# map boost stats classes to scipy class names and
# constructor arguments; b -> (s, ('ctor', 'args', ...))
_klass_mapper = {
    'beta': _KlassMap('beta', ('a', 'b')),
    'binomial': _KlassMap('binom', ('n', 'p')),
    'negative_binomial': _KlassMap('nbinom', ('n', 'p')),
}

# functions that take ctor params and parameter "x"
_x_funcs = ('pdf', 'cdf', 'sf', 'ppf', 'isf')

# functions that take only ctor params
_no_x_funcs = ('mean', 'variance', 'skewness', 'kurtosis_excess')


if __name__ == '__main__':
    # create target directory
    (pathlib.Path(__file__).parent / 'src').mkdir(exist_ok=True, parents=True)
    src_dir = pathlib.Path(__file__).parent / 'src'

    # copy contents of include into directory to satisfy Cython
    # PXD include conditions
    inc_dir = pathlib.Path(__file__).parent / 'include'
    src = 'templated_pyufunc.pxd'
    copyfile(inc_dir / src, src_dir / src)

    # generate the PXD and PYX wrappers
    from include.gen_func_defs_pxd import _gen_func_defs_pxd
    from include.code_gen import _ufunc_gen
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
