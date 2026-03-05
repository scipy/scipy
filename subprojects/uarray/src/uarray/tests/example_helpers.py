import uarray as ua


class _TypedBackend:
    __ua_domain__ = "ua_examples"

    def __init__(self, *my_types):
        self.my_types = my_types

    def __ua_convert__(self, dispatchables, coerce):
        if not all(type(d.value) in self.my_types for d in dispatchables):
            return NotImplemented
        return tuple(d.value for d in dispatchables)

    def __ua_function__(self, func, args, kwargs):
        return self.my_types[0]()


class TypeA:
    @classmethod
    def __repr__(cls):
        return cls.__name__


class TypeB(TypeA):
    pass


class TypeC(TypeA):
    pass


BackendA = _TypedBackend(TypeA)
BackendB = _TypedBackend(TypeB)
BackendC = _TypedBackend(TypeC)
BackendAB = _TypedBackend(TypeA, TypeB)
BackendBC = _TypedBackend(TypeB, TypeC)

creation_multimethod = ua.generate_multimethod(
    lambda: (), lambda a, kw, d: (a, kw), "ua_examples"
)
call_multimethod = ua.generate_multimethod(
    lambda *a: tuple(ua.Dispatchable(x, "mark") for x in a),
    lambda a, kw, d: (a, kw),
    "ua_examples",
)
