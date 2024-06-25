from ._biasedurn cimport CFishersNCHypergeometric, StochasticLib3
cimport numpy as np
import numpy as np
from libc.stdint cimport uint32_t, uint64_t
from libcpp.memory cimport unique_ptr

np.import_array()

from cpython.pycapsule cimport PyCapsule_GetPointer, PyCapsule_IsValid

cdef extern from "distributions.h":
    ctypedef struct bitgen_t:
        void *state
        uint64_t (*next_uint64)(void *st) nogil
        uint32_t (*next_uint32)(void *st) nogil
        double (*next_double)(void *st) nogil
        uint64_t (*next_raw)(void *st) nogil

    double random_normal(bitgen_t*, double, double) nogil

cdef class _PyFishersNCHypergeometric:
    cdef unique_ptr[CFishersNCHypergeometric] c_fnch

    def __cinit__(self, int n, int m, int N, double odds, double accuracy):
        self.c_fnch = unique_ptr[CFishersNCHypergeometric](new CFishersNCHypergeometric(n, m, N, odds, accuracy))

    def mode(self):
        return self.c_fnch.get().mode()

    def mean(self):
        return self.c_fnch.get().mean()

    def variance(self):
        return self.c_fnch.get().variance()

    def probability(self, int x):
        return self.c_fnch.get().probability(x)

    def moments(self):
        cdef double mean, var
        self.c_fnch.get().moments(&mean, &var)
        return mean, var


cdef class _PyWalleniusNCHypergeometric:
    cdef unique_ptr[CWalleniusNCHypergeometric] c_wnch

    def __cinit__(self, int n, int m, int N, double odds, double accuracy):
        self.c_wnch = unique_ptr[CWalleniusNCHypergeometric](new CWalleniusNCHypergeometric(n, m, N, odds, accuracy))

    def mode(self):
        return self.c_wnch.get().mode()

    def mean(self):
        return self.c_wnch.get().mean()

    def variance(self):
        return self.c_wnch.get().variance()

    def probability(self, int x):
        return self.c_wnch.get().probability(x)

    def moments(self):
        cdef double mean, var
        self.c_wnch.get().moments(&mean, &var)
        return mean, var


cdef bitgen_t* _glob_rng
cdef double next_double() noexcept nogil:
    global _glob_rng
    return _glob_rng.next_double(_glob_rng.state)
cdef double next_normal(const double m, const double s) noexcept nogil:
    global _glob_rng
    return random_normal(_glob_rng, m, s)

cdef object make_rng(random_state=None):
    # get a bit_generator object
    if random_state is None or isinstance(random_state, int):
        bg = np.random.RandomState(random_state)._bit_generator
    elif isinstance(random_state, np.random.RandomState):
        bg = random_state._bit_generator
    elif isinstance(random_state, np.random.Generator):
        bg = random_state.bit_generator
    else:
        raise ValueError('random_state is not one of None, int, RandomState, Generator')
    capsule = bg.capsule
    return capsule


cdef class _PyStochasticLib3:
    cdef unique_ptr[StochasticLib3] c_sl3
    cdef object capsule
    cdef bitgen_t* bit_generator

    def __cinit__(self):
        self.c_sl3 = unique_ptr[StochasticLib3](new StochasticLib3(0))
        self.c_sl3.get().next_double = &next_double
        self.c_sl3.get().next_normal = &next_normal

    def Random(self):
        return self.c_sl3.get().Random()

    def SetAccuracy(self, double accur):
        return self.c_sl3.get().SetAccuracy(accur)

    cdef void HandleRng(self, random_state=None) noexcept:
        self.capsule = make_rng(random_state)

        # get the bitgen_t pointer
        cdef const char *capsule_name = "BitGenerator"
        if not PyCapsule_IsValid(self.capsule, capsule_name):
            raise ValueError("Invalid pointer to anon_func_state")
        self.bit_generator = <bitgen_t *> PyCapsule_GetPointer(self.capsule, capsule_name)
        global _glob_rng
        _glob_rng = self.bit_generator

    def rvs_fisher(self, int n, int m, int N, double odds, int size, random_state=None):
        # handle random state
        self.HandleRng(random_state)

        # call for each
        rvs = np.empty(size, dtype=np.float64)
        for ii in range(size):
            rvs[ii] = self.c_sl3.get().FishersNCHyp(n, m, N, odds)
        return rvs

    def rvs_wallenius(self, int n, int m, int N, double odds, int size, random_state=None):
        # handle random state
        self.HandleRng(random_state)

        # call for each
        rvs = np.empty(size, dtype=np.float64)
        for ii in range(size):
            rvs[ii] = self.c_sl3.get().WalleniusNCHyp(n, m, N, odds)
        return rvs

    def FishersNCHyp(self, int n, int m, int N, double odds):
        self.HandleRng(None)  # get default rng
        return self.c_sl3.get().FishersNCHyp(n, m, N, odds)

    def WalleniusNCHyp(self, int n, int m, int N, double odds):
        self.HandleRng(None)  # get default rng
        return self.c_sl3.get().WalleniusNCHyp(n, m, N, odds)
