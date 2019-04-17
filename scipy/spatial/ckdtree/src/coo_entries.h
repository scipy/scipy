#ifndef CKDTREE_COO_ENTRIES
#define CKDTREE_COO_ENTRIES

#include <Python.h>
#include "numpy/arrayobject.h"

struct coo_entry {
    npy_intp i;
    npy_intp j;
    npy_float64 v;
};

#endif

