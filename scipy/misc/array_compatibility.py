import numpy as np


def get_namespace(*xs):
    # `xs` contains one or more arrays, or possibly Python scalars (accepting
    # those is a matter of taste, but doesn't seem unreasonable).
    namespaces = {
        x.__array_namespace__() if hasattr(x, '__array_namespace__') else None
        for x in xs if not isinstance(x, (bool, int, float, complex))
    }

    if not namespaces:
        # one could special-case np.ndarray above or use np.asarray here if
        # older numpy versions need to be supported.
        raise ValueError("Unrecognized array input")

    if len(namespaces) != 1:
        raise ValueError(f"Multiple namespaces for array inputs: {namespaces}")

    xp, = namespaces
    if xp is None:
        # Use numpy as default
        return np, False
    return xp, True
