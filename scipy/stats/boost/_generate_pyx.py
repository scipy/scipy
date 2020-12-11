import pathlib
import inspect
from shutil import copyfile

import scipy.stats

# map boost stats classes to scipy classes; b -> s
klass_mapper = {
    'bernoulli': scipy.stats.bernoulli,
    'beta': scipy.stats.beta,
    'binomial': scipy.stats.binom,
    'negative_binomial': scipy.stats.nbinom,
    'non_central_chi_squared': scipy.stats.ncx2,
}


if __name__ == '__main__':
    # create target directory
    (pathlib.Path(__file__).parent / 'src').mkdir(exist_ok=True, parents=True)
    src_dir = pathlib.Path(__file__).parent / 'src'

    # copy contents of include into directory
    inc_dir = pathlib.Path(__file__).parent / 'include'
    src = 'templated_pyufunc.pxd'
    copyfile(inc_dir / src, src_dir / src)

    # generate the PXD and PYX wrappers
    from include.gen_func_defs_pxd import gen_func_defs_pxd
    from include.code_gen import ufunc_gen
    gen_func_defs_pxd(f'{src_dir}/func_defs.pxd')
    for b, s in klass_mapper.items():
        if s is None:
            print(f'{b} has no scipy equivalent! Skipping!')
            continue
        # cheat for now by pulling constructors from scipy classes
        scipy_name = s.__class__.__name__.split('_gen')[0]
        ctor_args = [p for p in inspect.signature(s._stats).parameters if p != 'moments']
        ufunc_gen(
            wrapper_prefix=scipy_name,
            types=['NPY_DOUBLE', 'NPY_FLOAT'],
            num_ctor_args=len(ctor_args),
            filename=f'{src_dir}/{scipy_name}_ufunc.pyx',
            boost_dist=f'{b}_distribution',
        )
