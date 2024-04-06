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
    'hypergeometric': _KlassMap('hypergeom', ('r', 'n', 'N')),
    'non_central_f': _KlassMap('ncf', ('dfn', 'dfd', 'nc')),
    'non_central_chi_squared': _KlassMap('ncx2', ('df', 'nc')),
    'non_central_t': _KlassMap('nct', ('df', 'nc')),
    'skew_normal': _KlassMap('skewnorm', ('loc', 'scale', 'a',)),
    'inverse_gaussian': _KlassMap('invgauss', ('mu', 'mean')),
}

# functions that take ctor params and parameter "x"
_x_funcs = ('pdf', 'cdf', 'sf', 'ppf', 'isf')

# functions that take only ctor params
_no_x_funcs = ('mean', 'variance', 'skewness', 'kurtosis_excess')
