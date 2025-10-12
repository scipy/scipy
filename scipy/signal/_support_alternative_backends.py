
from scipy._lib._array_api import SCIPY_ARRAY_API


class JaxBackend:
    name = "jax"
    # A class, just for convenience (may need to change)
    primary_types = ["~jax:Array"]  # allow subclasses otherwise need _jax.ArrayImpl?
    secondary_types = []
    requires_opt_in = False
    # Patched below while decorating we add the function
    # (unless it is blocklisted)
    functions = {}


class CupyBackend:
    name = "cupy"
    # A class, just for convenience (may need to change)
    primary_types = ["cupy:ndarray"]
    secondary_types = []
    requires_opt_in = False
    # Patched below while decorating we add the function
    # (unless it is blocklisted)
    RENAMES = {'freqz_sos': 'sosfreqz'}
    functions = {}

if SCIPY_ARRAY_API:
    from scipy._lib.spatch.backend_system import BackendSystem

    _bs = BackendSystem(
        None,  # don't load entry-points for now (no 3rd party backends)
        "_SCIPY_INTERNAL_BACKENDS",  # spatch env-var prefix
        default_primary_types=["numpy:ndarray"],
        backends=[JaxBackend],
    )

    def _dispatchable(*args, cupy=True, jax=False, **kwargs):
        # Cupy supports most (but some are incompatible) while JAX only
        # supports a select few.  As these are builtin, we build the list
        # of supported functions on the fly here (normally hardcoded).
        def decorator(func):
            name = func.__qualname__
            if cupy:
                cupy_name = CupyBackend.RENAMES.get(name, name)
                CupyBackend.functions[name] = {
                    "function": f"cupyx.signal:{name}"
                }
            if jax:
                JaxBackend.functions[name] = {
                    "function": f"jax.scipy.signal:{cupy_name}"
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
