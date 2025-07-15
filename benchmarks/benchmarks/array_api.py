from .common import XPBenchmark, safe_import

with safe_import():
    from scipy._lib.array_api_compat import array_namespace as compat_namespace
    from scipy._lib._array_api import array_namespace


class ArrayNamespace(XPBenchmark):
    def setup(self, backend):
        def f(x):
            _ = array_namespace(x)
            return x
    
        super().setup(backend, f)
        self.x = self.synchronize(self.xp.empty(0))
        # Populate @lru_cache and jax.jit. Note that this benefits all backends.
        self.func(self.x)

    def time_array_namespace(self, backend):
        """scipy wrapper around array_api_compat.array_namespace"""
        array_namespace(self.x)

    def time_compat_namespace(self, backend):
        """Bare array_api_compat.array_namespace"""
        compat_namespace(self.x)

    def time_trivial_func(self, backend):
        """Trivial function that internally calls `xp=array_namespace(*args)`"""
        self.func(self.x)
