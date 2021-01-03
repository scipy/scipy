from numpy.random cimport bitgen_t

# Declare the class with cdef
cdef extern from "biasedurn/stocc.h" nogil:
    cdef cppclass CFishersNCHypergeometric:
        CFishersNCHypergeometric() except +
        CFishersNCHypergeometric(int, int, int, double, double) except +
        int mode()
        double mean()
        double variance()
        double probability(int x)
        double moments(double * mean, double * var)

    cdef cppclass CWalleniusNCHypergeometric:
        CWalleniusNCHypergeometric() except +
        CWalleniusNCHypergeometric(int, int, int, double, double) except +
        int mode()
        double mean()
        double variance()
        double probability(int x)
        double moments(double * mean, double * var)

    cdef cppclass StochasticLib3:
        StochasticLib3() except +
        StochasticLib3(int seed) except +
        void SetBitGen(bitgen_t *that_bitgen_state)
        double Random() except +
        void SetAccuracy(double accur)
        int FishersNCHyp (int n, int m, int N, double odds) except +
        int WalleniusNCHyp (int n, int m, int N, double odds) except +
