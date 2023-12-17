import numpy as np

from .common import run_monitored, set_mem_rlimit, Benchmark, safe_import

with safe_import():
    from scipy.stats import spearmanr

with safe_import():
    import scipy.interpolate as interpolate


class Leaks(Benchmark):
    unit = "relative increase with repeats"

    def track_leaks(self):
        set_mem_rlimit()

        # Setup temp file, make it fit in memory
        repeats = [2, 5, 10, 50, 200]
        peak_mems = []

        for repeat in repeats:
            code = """
            import numpy as np
            from scipy.interpolate import griddata

            def func(x, y):
                return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

            grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
            points = np.random.rand(1000, 2)
            values = func(points[:,0], points[:,1])

            for t in range(%(repeat)d):
                for method in ['nearest', 'linear', 'cubic']:
                    griddata(points, values, (grid_x, grid_y), method=method)

            """ % dict(repeat=repeat)

            _, peak_mem = run_monitored(code)
            peak_mems.append(peak_mem)

        corr, p = spearmanr(repeats, peak_mems)
        if p < 0.05:
            print("*"*79)
            print("PROBABLE MEMORY LEAK")
            print("*"*79)
        else:
            print("PROBABLY NO MEMORY LEAK")

        return max(peak_mems) / min(peak_mems)


class BenchPPoly(Benchmark):

    def setup(self):
        rng = np.random.default_rng(1234)
        m, k = 55, 3
        x = np.sort(rng.random(m+1))
        c = rng.random((k, m))
        self.pp = interpolate.PPoly(c, x)

        npts = 100
        self.xp = np.linspace(0, 1, npts)

    def time_evaluation(self):
        self.pp(self.xp)


class GridData(Benchmark):
    param_names = ['n_grids', 'method']
    params = [
        [10j, 100j, 1000j],
        ['nearest', 'linear', 'cubic']
    ]

    def setup(self, n_grids, method):
        self.func = lambda x, y: x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
        self.grid_x, self.grid_y = np.mgrid[0:1:n_grids, 0:1:n_grids]
        self.points = np.random.rand(1000, 2)
        self.values = self.func(self.points[:, 0], self.points[:, 1])

    def time_evaluation(self, n_grids, method):
        interpolate.griddata(self.points, self.values, (self.grid_x, self.grid_y),
                             method=method)


class Interpolate1d(Benchmark):
    param_names = ['n_samples', 'method']
    params = [
        [10, 50, 100, 1000, 10000],
        ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic'],
    ]

    def setup(self, n_samples, method):
        self.x = np.arange(n_samples)
        self.y = np.exp(-self.x/3.0)
        self.interpolator = interpolate.interp1d(self.x, self.y, kind=method)
        self.xp = np.linspace(self.x[0], self.x[-1], 4*n_samples)

    def time_interpolate(self, n_samples, method):
        """Time the construction overhead."""
        interpolate.interp1d(self.x, self.y, kind=method)

    def time_interpolate_eval(self, n_samples, method):
        """Time the evaluation."""
        self.interpolator(self.xp)


class Interpolate2d(Benchmark):
    param_names = ['n_samples', 'method']
    params = [
        [10, 50, 100],
        ['linear', 'cubic', 'quintic'],
    ]

    def setup(self, n_samples, method):
        r_samples = n_samples / 2.
        self.x = np.arange(-r_samples, r_samples, 0.25)
        self.y = np.arange(-r_samples, r_samples, 0.25)
        self.xx, self.yy = np.meshgrid(self.x, self.y)
        self.z = np.sin(self.xx**2+self.yy**2)

    def time_interpolate(self, n_samples, method):
        interpolate.interp2d(self.x, self.y, self.z, kind=method)


class Rbf(Benchmark):
    param_names = ['n_samples', 'function']
    params = [
        [10, 50, 100],
        ['multiquadric', 'inverse', 'gaussian', 'linear',
         'cubic', 'quintic', 'thin_plate']
    ]

    def setup(self, n_samples, function):
        self.x = np.arange(n_samples)
        self.y = np.sin(self.x)
        r_samples = n_samples / 2.
        self.X = np.arange(-r_samples, r_samples, 0.25)
        self.Y = np.arange(-r_samples, r_samples, 0.25)
        self.z = np.exp(-self.X**2-self.Y**2)

    def time_rbf_1d(self, n_samples, function):
        interpolate.Rbf(self.x, self.y, function=function)

    def time_rbf_2d(self, n_samples, function):
        interpolate.Rbf(self.X, self.Y, self.z, function=function)


class RBFInterpolator(Benchmark):
    param_names = ['neighbors', 'n_samples', 'kernel']
    params = [
        [None, 50],
        [10, 100, 1000],
        ['linear', 'thin_plate_spline', 'cubic', 'quintic', 'multiquadric',
         'inverse_multiquadric', 'inverse_quadratic', 'gaussian']
    ]

    def setup(self, neighbors, n_samples, kernel):
        rng = np.random.RandomState(0)
        self.y = rng.uniform(-1, 1, (n_samples, 2))
        self.x = rng.uniform(-1, 1, (n_samples, 2))
        self.d = np.sum(self.y, axis=1)*np.exp(-6*np.sum(self.y**2, axis=1))

    def time_rbf_interpolator(self, neighbors, n_samples, kernel):
        interp = interpolate.RBFInterpolator(
            self.y,
            self.d,
            neighbors=neighbors,
            epsilon=5.0,
            kernel=kernel
            )
        interp(self.x)


class UnivariateSpline(Benchmark):
    param_names = ['n_samples', 'degree']
    params = [
        [10, 50, 100],
        [3, 4, 5]
    ]

    def setup(self, n_samples, degree):
        r_samples = n_samples / 2.
        self.x = np.arange(-r_samples, r_samples, 0.25)
        self.y = np.exp(-self.x**2) + 0.1 * np.random.randn(*self.x.shape)

    def time_univariate_spline(self, n_samples, degree):
        interpolate.UnivariateSpline(self.x, self.y, k=degree)


class BivariateSpline(Benchmark):
    """
    Author: josef-pktd and scipy mailinglist example
    'http://scipy-user.10969.n7.nabble.com/BivariateSpline-examples\
    -and-my-crashing-python-td14801.html'
    """
    param_names = ['n_samples']
    params = [
        [10, 20, 30]
    ]

    def setup(self, n_samples):
        x = np.arange(0, n_samples, 0.5)
        y = np.arange(0, n_samples, 0.5)
        x, y = np.meshgrid(x, y)
        x = x.ravel()
        y = y.ravel()
        xmin = x.min()-1
        xmax = x.max()+1
        ymin = y.min()-1
        ymax = y.max()+1
        s = 1.1
        self.yknots = np.linspace(ymin+s, ymax-s, 10)
        self.xknots = np.linspace(xmin+s, xmax-s, 10)
        self.z = np.sin(x) + 0.1*np.random.normal(size=x.shape)
        self.x = x
        self.y = y

    def time_smooth_bivariate_spline(self, n_samples):
        interpolate.SmoothBivariateSpline(self.x, self.y, self.z)

    def time_lsq_bivariate_spline(self, n_samples):
        interpolate.LSQBivariateSpline(self.x, self.y, self.z,
                                       self.xknots.flat, self.yknots.flat)


class Interpolate(Benchmark):
    """
    Linear Interpolate in scipy and numpy
    """
    param_names = ['n_samples', 'module']
    params = [
        [10, 50, 100],
        ['numpy', 'scipy']
    ]

    def setup(self, n_samples, module):
        self.x = np.arange(n_samples)
        self.y = np.exp(-self.x/3.0)
        self.z = np.random.normal(size=self.x.shape)

    def time_interpolate(self, n_samples, module):
        if module == 'scipy':
            interpolate.interp1d(self.x, self.y, kind="linear")
        else:
            np.interp(self.z, self.x, self.y)


class RegularGridInterpolator(Benchmark):
    """
    Benchmark RegularGridInterpolator with method="linear".
    """
    param_names = ['ndim', 'max_coord_size', 'n_samples', 'flipped']
    params = [
        [2, 3, 4],
        [10, 40, 200],
        [10, 100, 1000, 10000],
        [1, -1]
    ]

    def setup(self, ndim, max_coord_size, n_samples, flipped):
        rng = np.random.default_rng(314159)

        # coordinates halve in size over the dimensions
        coord_sizes = [max_coord_size // 2**i for i in range(ndim)]
        self.points = [np.sort(rng.random(size=s))[::flipped]
                       for s in coord_sizes]
        self.values = rng.random(size=coord_sizes)

        # choose in-bounds sample points xi
        bounds = [(p.min(), p.max()) for p in self.points]
        xi = [rng.uniform(low, high, size=n_samples)
              for low, high in bounds]
        self.xi = np.array(xi).T

        self.interp = interpolate.RegularGridInterpolator(
            self.points,
            self.values,
        )

    def time_rgi_setup_interpolator(self, ndim, max_coord_size,
                                    n_samples, flipped):
        self.interp = interpolate.RegularGridInterpolator(
            self.points,
            self.values,
        )

    def time_rgi(self, ndim, max_coord_size, n_samples, flipped):
        self.interp(self.xi)


class RegularGridInterpolatorValues(interpolate.RegularGridInterpolator):
    def __init__(self, points, xi, **kwargs):
        # create fake values for initialization
        values = np.zeros(tuple([len(pt) for pt in points]))
        super().__init__(points, values, **kwargs)
        self._is_initialized = False
        # precompute values
        (self.xi, self.xi_shape, self.ndim,
         self.nans, self.out_of_bounds) = self._prepare_xi(xi)
        self.indices, self.norm_distances = self._find_indices(xi.T)
        self._is_initialized = True

    def _prepare_xi(self, xi):
        if not self._is_initialized:
            return super()._prepare_xi(xi)
        else:
            # just give back precomputed values
            return (self.xi, self.xi_shape, self.ndim,
                    self.nans, self.out_of_bounds)

    def _find_indices(self, xi):
        if not self._is_initialized:
            return super()._find_indices(xi)
        else:
            # just give back pre-computed values
            return self.indices, self.norm_distances

    def __call__(self, values, method=None):
        values = self._check_values(values)
        # check fillvalue
        self._check_fill_value(values, self.fill_value)
        # check dimensionality
        self._check_dimensionality(self.grid, values)
        # flip, if needed
        self.values = np.flip(values, axis=self._descending_dimensions)
        return super().__call__(self.xi, method=method)


class RegularGridInterpolatorSubclass(Benchmark):
    """
    Benchmark RegularGridInterpolator with method="linear".
    """
    param_names = ['ndim', 'max_coord_size', 'n_samples', 'flipped']
    params = [
        [2, 3, 4],
        [10, 40, 200],
        [10, 100, 1000, 10000],
        [1, -1]
    ]

    def setup(self, ndim, max_coord_size, n_samples, flipped):
        rng = np.random.default_rng(314159)

        # coordinates halve in size over the dimensions
        coord_sizes = [max_coord_size // 2**i for i in range(ndim)]
        self.points = [np.sort(rng.random(size=s))[::flipped]
                       for s in coord_sizes]
        self.values = rng.random(size=coord_sizes)

        # choose in-bounds sample points xi
        bounds = [(p.min(), p.max()) for p in self.points]
        xi = [rng.uniform(low, high, size=n_samples)
              for low, high in bounds]
        self.xi = np.array(xi).T

        self.interp = RegularGridInterpolatorValues(
            self.points,
            self.xi,
        )

    def time_rgi_setup_interpolator(self, ndim, max_coord_size,
                                    n_samples, flipped):
        self.interp = RegularGridInterpolatorValues(
            self.points,
            self.xi,
        )

    def time_rgi(self, ndim, max_coord_size, n_samples, flipped):
        self.interp(self.values)


class CloughTocherInterpolatorValues(interpolate.CloughTocher2DInterpolator):
    """Subclass of the CT2DInterpolator with optional `values`.

    This is mainly a demo of the functionality. See
    https://github.com/scipy/scipy/pull/18376 for discussion
    """
    def __init__(self, points, xi, tol=1e-6, maxiter=400, **kwargs):
        interpolate.CloughTocher2DInterpolator.__init__(self, points, None,
                                                        tol=tol, maxiter=maxiter)
        self.xi = None
        self._preprocess_xi(*xi)
        self.simplices, self.c = (
            interpolate.CloughTocher2DInterpolator._find_simplicies(self, self.xi)
        )

    def _preprocess_xi(self, *args):
        if self.xi is None:
            self.xi, self.interpolation_points_shape = (
                interpolate.CloughTocher2DInterpolator._preprocess_xi(self, *args)
            )
        return self.xi, self.interpolation_points_shape
    
    def _find_simplicies(self, xi):
        return self.simplices, self.c

    def __call__(self, values):
        self._set_values(values)
        return super().__call__(self.xi)


class CloughTocherInterpolatorSubclass(Benchmark):
    """
    Benchmark CloughTocherInterpolatorValues.

    Derived from the docstring example,
    https://docs.scipy.org/doc/scipy-1.11.2/reference/generated/scipy.interpolate.CloughTocher2DInterpolator.html
    """
    param_names = ['n_samples']
    params = [10, 50, 100]

    def setup(self, n_samples):
        rng = np.random.default_rng(314159)

        x = rng.random(n_samples) - 0.5
        y = rng.random(n_samples) - 0.5


        self.z = np.hypot(x, y)
        X = np.linspace(min(x), max(x))
        Y = np.linspace(min(y), max(y))
        self.X, self.Y = np.meshgrid(X, Y)

        self.interp = CloughTocherInterpolatorValues(
            list(zip(x, y)), (self.X, self.Y)
        )

    def time_clough_tocher(self, n_samples):
            self.interp(self.z)

