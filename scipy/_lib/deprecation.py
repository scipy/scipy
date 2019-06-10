import warnings

__all__ = ["_deprecated"]


def _deprecated(msg):

    def wrap(fun):
        def call(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return fun(*args, **kwargs)
        call.__doc__ = msg
        return call

    return wrap
