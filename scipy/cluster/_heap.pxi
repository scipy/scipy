# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np


cdef class Heap:
    cdef int[:] index_by_key
    cdef int[:] key_by_index
    cdef double[:] values
    cdef int size

    def __init__(self, double[:] values):
        self.size = values.shape[0]
        self.index_by_key = np.arange(self.size, dtype=np.intc)
        self.key_by_index = np.arange(self.size, dtype=np.intc)
        self.values = values.copy()
        cdef int i
        for i in reversed(range(self.size / 2)):
            self.sift_down(i)

    cdef get_min(self):
        return self.key_by_index[0], self.values[0]

    cdef remove_min(self):
        self.swap(0, self.size - 1)
        self.size -= 1
        self.sift_down(0)

    cdef change_value(self, int key, double value):
        cdef int index = self.index_by_key[key]
        cdef double old_value = self.values[index]
        self.values[index] = value
        if value < old_value:
            self.sift_up(index)
        else:
            self.sift_down(index)

    cdef sift_up(self, int index):
        cdef int parent = Heap.parent(index)
        while index > 0 and self.values[parent] > self.values[index]:
            self.swap(index, parent)
            index = parent
            parent = Heap.parent(index)

    cdef sift_down(self, int index):
        cdef int child = Heap.left_child(index)
        while child < self.size:
            if (child + 1 < self.size and
                    self.values[child + 1] < self.values[child]):
                child += 1

            if self.values[index] > self.values[child]:
                self.swap(index, child)
                index = child
                child = Heap.left_child(index)
            else:
                break

    @staticmethod
    cdef left_child(int parent):
        return (parent << 1) + 1

    @staticmethod
    cdef parent(int child):
        return (child - 1) >> 1

    cdef swap(self, int i, int j):
        self.values[i], self.values[j] = self.values[j], self.values[i]
        cdef int key_i = self.key_by_index[i]
        cdef int key_j = self.key_by_index[j]
        self.key_by_index[i] = key_j
        self.key_by_index[j] = key_i
        self.index_by_key[key_i] = j
        self.index_by_key[key_j] = i
