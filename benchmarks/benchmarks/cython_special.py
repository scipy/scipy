import re
import numpy as np
from scipy import special

from .common import with_attributes, safe_import

with safe_import():
    from scipy.special import cython_special


FUNC_ARGS = {
    'airy_d': (1,),
    'airy_D': (1,),
    'beta_dd': (0.25, 0.75),
    'erf_d': (1,),
    'erf_D': (1+1j,),
    'exprel_d': (1e-6,),
    'gamma_d': (100,),
    'gamma_D': (100+100j,),
    'jv_dd': (1, 1),
    'jv_dD': (1, (1+1j)),
    'loggamma_D': (20,),
    'logit_d': (0.5,),
    'psi_d': (1,),
    'psi_D': (1,),
}


class _CythonSpecialMeta(type):
    """
    Add time_* benchmarks corresponding to cython_special._bench_*_cy
    """

    def __new__(cls, cls_name, bases, dct):
        params = [(10, 100, 1000), ('python', 'numpy', 'cython')]
        param_names = ['N', 'api']

        def get_time_func(name, args):

            @with_attributes(params=[(name,), (args,)] + params,
                             param_names=['name', 'argument'] + param_names)
            def func(self, name, args, N, api):
                if api == 'python':
                    self.py_func(N, *args)
                elif api == 'numpy':
                    self.np_func(*self.obj)
                else:
                    self.cy_func(N, *args)

            func.__name__ = 'time_' + name
            return func

        for name in FUNC_ARGS.keys():
            func = get_time_func(name, FUNC_ARGS[name])
            dct[func.__name__] = func

        return type.__new__(cls, cls_name, bases, dct)


class CythonSpecial(metaclass=_CythonSpecialMeta):
    def setup(self, name, args, N, api):
        self.py_func = getattr(cython_special, '_bench_{}_py'.format(name))
        self.cy_func = getattr(cython_special, '_bench_{}_cy'.format(name))
        m = re.match('^(.*)_[dDl]+$', name)
        self.np_func = getattr(special, m.group(1))

        self.obj = []
        for arg in args:
            self.obj.append(arg*np.ones(N))
        self.obj = tuple(self.obj)
