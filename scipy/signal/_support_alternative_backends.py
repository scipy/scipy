
from scipy._lib._array_api import SCIPY_ARRAY_API


class JaxBackend:
    name = "jax"
    primary_types = ["~jax:Array"]  # allow subclasses otherwise need _jax.ArrayImpl?
    secondary_types = []
    requires_opt_in = False
    # Patched below while decorating we add the function
    # based on how the original is decorated.
    functions = {}


class CupyBackend:
    name = "cupy"
    primary_types = ["~cupy:ndarray"]
    secondary_types = []
    requires_opt_in = False
    # Patched below while decorating we add the function
    # based on how the original is decorated.
    functions = {}


if SCIPY_ARRAY_API:
    from scipy._lib.spatch.backend_system import BackendSystem

    _bs = BackendSystem(
        None,  # don't load entry-points for now (no 3rd party backends)
        "_SCIPY_INTERNAL_BACKENDS",  # spatch env-var prefix
        default_primary_types=["numpy:ndarray"],
        backends=[CupyBackend, JaxBackend],
    )

    def _dispatchable(*args, cupy=True, jax=False, **kwargs):
        # Cupy supports most (but some are incompatible) while JAX only
        # supports a select few.  As these are builtin, we build the list
        # of supported functions on the fly here (normally hardcoded).
        def decorator(func):
            name = func.__qualname__
            if cupy:
                cupy_name = cupy if isinstance(cupy, str) else name
                CupyBackend.functions[f"scipy.signal:{name}"] = {
                    "function": f"cupyx.scipy.signal:{cupy_name}",
                }
            if jax:
                jax_name = jax if isinstance(jax, str) else name
                JaxBackend.functions[f"scipy.signal:{name}"] = {
                    "function": f"jax.scipy.signal:{jax_name}",
                }

            # Use the normal spatch decoration now, but change module
            # to scipy.signal (to have a nicer name for mapping).
            return _bs.dispatchable(
                *args, module="scipy.signal", **kwargs)(func)

        return decorator

else:
    def _do_nothing(func):
        return func

    def _dispatchable(*args, **kwargs):
        return _do_nothing
