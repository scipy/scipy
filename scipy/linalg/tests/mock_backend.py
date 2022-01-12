import numpy as np


class _MockFunction:
    def __init__(self, return_value=None):
        self.number_calls = 0
        self.return_value = return_value
        self.last_args = ([], {})

    def __call__(self, *args, **kwargs):
        self.number_calls += 1
        self.last_args = (args, kwargs)
        return self.return_value


method_names = [
    # solvers
    'solve_sylvester',
    'solve_continuous_lyapunov', 'solve_discrete_lyapunov',
    'solve_lyapunov',
    'solve_continuous_are', 'solve_discrete_are'
]

for name in method_names:
    globals()[name] = _MockFunction(np.array([[0, 0], [1, 1]]))


__ua_domain__ = "numpy.scipy.linalg"


def __ua_function__(method, args, kwargs):
    fn = globals().get(method.__name__)
    return (fn(*args, **kwargs) if fn is not None else NotImplemented)
