# -*- cython -*-
"""
Many Scipy special functions originally cast silently double input
arguments to integers.

Here, we define such unsafe wrappers manually.

"""

cdef extern from "cephes.h":
    double bdtrc(int k, int n, double p) nogil
    double bdtr(int k, int n, double p) nogil 
    double bdtri(int k, int n, double y) nogil
    double expn(int n, double x) nogil
    double hyp2f0(double a, double b, double x, int type, double *err) nogil
    double nbdtrc(int k, int n, double p) nogil
    double nbdtr(int k, int n, double p) nogil 
    double nbdtri(int k, int n, double p) nogil
    double pdtrc(int k, double m) nogil
    double pdtr(int k, double m) nogil
    double pdtri(int k, double y) nogil
    double kn(int n, double x) nogil
    double yn(int n, double x) nogil
    double smirnov(int n, double e) nogil
    double smirnovi(int n, double p) nogil

cdef inline double bdtrc_unsafe(double k, double n, double p) nogil:
    return bdtrc(<int>k, <int>n, p)
cdef inline double bdtr_unsafe(double k, double n, double p) nogil:
    return bdtr(<int>k, <int>n, p)
cdef inline double bdtri_unsafe(double k, double n, double y) nogil:
    return bdtri(<int>k, <int>n, y)
cdef inline double expn_unsafe(double n, double x) nogil:
    return expn(<int>n, x)
cdef inline double hyp2f0_unsafe(double a, double b, double x, double type, double *err) nogil:
    return hyp2f0(a, b, x, <int>type, err)
cdef inline double nbdtrc_unsafe(double k, double n, double p) nogil:
    return nbdtrc(<int>k, <int>n, p)
cdef inline double nbdtr_unsafe(double k, double n, double p)  nogil:
    return nbdtr(<int>k, <int>n, p)
cdef inline double nbdtri_unsafe(double k, double n, double p) nogil:
    return nbdtri(<int>k, <int>n, p)
cdef inline double pdtrc_unsafe(double k, double m) nogil:
    return pdtrc(<int>k, m)
cdef inline double pdtr_unsafe(double k, double m) nogil:
    return pdtr(<int>k, m)
cdef inline double pdtri_unsafe(double k, double y) nogil:
    return pdtri(<int>k, y)
cdef inline double kn_unsafe(double n, double x) nogil:
    return kn(<int>n, x)
cdef inline double yn_unsafe(double n, double x) nogil:
    return yn(<int>n, x)
cdef inline double smirnov_unsafe(double n, double e) nogil:
    return smirnov(<int>n, e)
cdef inline double smirnovi_unsafe(double n, double p) nogil:
    return smirnovi(<int>n, p)
