import numpy as np
from mpmath import mp

mp.dps = 100

class Distribution:
    @staticmethod
    def _make_mpf_array(x):
        shape = np.shape(x) or (1,)
        x = np.asarray(x, dtype=np.float64).ravel()
        return np.asarray([mp.mpf(xi) for xi in x]).reshape(shape)

    def __init__(self, **kwargs):
        _shape_info = self._shape_info()

        invalid_params = set(kwargs) - set(_shape_info)
        if invalid_params:
            raise ValueError(f"Invalid parameters {invalid_params}.")

        missing_params = set(_shape_info) - set(kwargs)
        if missing_params:
            raise ValueError(f"Parameters {missing_params} not specified.")

        shapes = []
        for key, val in kwargs.items():
            shapes.append(np.shape(val))
        self._param_shape = np.broadcast_shapes(*shapes)

        for key, val in kwargs.items():
            val = np.broadcast_to(self._make_mpf_array(val),
                                  self._param_shape or (1,))
            setattr(self, key, val)

        self._check_domain()

    def _check_domain(self):
        # add check for relationship between shapes
        _shape_info = self._shape_info()
        conditions = []
        for param_name in _shape_info:
            param_val = getattr(self, param_name)
            a, b = _shape_info[param_name]['domain']
            a_inclusive, b_inclusive = _shape_info[param_name]['inclusive']
            a_compare = np.less_equal if a_inclusive else np.less
            b_compare = np.less_equal if b_inclusive else np.less
            conditions.append(a_compare(a, param_val)
                              & b_compare(param_val, b))
        self._in_domain = np.all(conditions, axis=0)

    def _result_shape(self, x=None):
        return (np.broadcast_shapes(np.shape(x), self._param_shape)
                if x is not None else self._param_shape)

    def support(self):
        a, b = self._support()
        return (np.broadcast_to(a, self._param_shape),
                np.broadcast_to(b, self._param_shape))

    def _support(self):
        return -mp.inf, mp.inf

    def supported(self, x):
        # TODO: add inclusive support
        a, b = self.support()
        return (a < x) & (x < b)


class Normal(Distribution):

    def _shape_info(self):
        return {}

    def _pdf(self, x):
        return mp.npdf(x)

class SkewNormal(Distribution):

    def _shape_info(self):
        return {'a': dict(domain=(-mp.inf, mp.inf), inclusive=(False, False))}

    def _pdf(self, x):
        return 2 * mp.npdf(x) * mp.ncdf(self.a * x)

class Beta(Distribution):
    def _support(self):
        return 0, 1

    def _shape_info(self):
        return {'a': dict(domain=(0, 1), inclusive=(False, False)),
                'b': dict(domain=(0, 1), inclusive=(False, False))}

dist = Normal()
# SkewNormal(a=0.5)
