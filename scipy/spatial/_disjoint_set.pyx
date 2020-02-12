import numpy as np
cimport numpy as np
cimport cython


__all__ = ['DisjointSet']


cdef class DisjointSet:
    cdef:
        readonly np.npy_intp[:] parents
        readonly np.npy_intp[:] sizes

    def __init__(DisjointSet self, np.intp_t n):
        self.sizes = np.ones(n, dtype=int)
        self.parents = np.arange(n)

    def find(DisjointSet self, np.intp_t index):
        cdef np.npy_intp parent
        parents = self.parents
        parent = parents[index]
        while parent != parents[parent]:
            parent = parents[parent]
        parents[index] = parent
        return parent

    def merge(DisjointSet self, np.intp_t a, np.intp_t b):
        a = self.find(a)
        b = self.find(b)
        if a == b:
            return False

        sizes = self.sizes
        parents = self.parents
        if sizes[a] < sizes[b]:
            parents[a] = b
            sizes[b] += sizes[a]
        else:
            parents[b] = a
            sizes[a] += sizes[b]
        return True
