# Declare the class with cdef
cdef extern from "biasedurn/stocc.h" nogil:
    cdef cppclass CFishersNCHypergeometric:
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

cdef extern from "biasedurn/stocc.h":
    cdef cppclass StochasticLib3:
        StochasticLib3(int seed) nogil except +
        double Random() nogil except +
        void SetAccuracy(double accur) nogil
        int FishersNCHyp (int n, int m, int N, double odds) nogil except +
        int WalleniusNCHyp (int n, int m, int N, double odds) nogil except +
        double(*next_double)() nogil
