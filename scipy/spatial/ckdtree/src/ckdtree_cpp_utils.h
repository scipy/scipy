#include <vector>
#include <cstring>
#include <Python.h>

#include "ckdtree_cpp_ordered_pair.h"
#include "ckdtree_cpp_decl.h"
#include "ckdtree_cpp_exc.h"


#if PY_MAJOR_VERSION < 3
    #define ckdtree_PyBytes_FromStringAndSize(v,len) PyString_FromStringAndSize(v,len)
    #define ckdtree_PyBytes_Size(o) PyString_Size(o)
    #define ckdtree_PyBytes_AsString(o) PyString_AsString(o)
#else
    #define ckdtree_PyBytes_FromStringAndSize(v,len) PyBytes_FromStringAndSize(v,len)
    #define ckdtree_PyBytes_Size(o) PyBytes_Size(o)
    #define ckdtree_PyBytes_AsString(o) PyBytes_AsString(o)
#endif


inline void*
tree_buffer_pointer(std::vector<ckdtreenode> *buf)
{
    std::vector<ckdtreenode> &tmp = *buf;
    return (void*)&tmp[0];
}


inline ckdtreenode*
tree_buffer_root(std::vector<ckdtreenode> *buf)
{
    std::vector<ckdtreenode> &tmp = *buf;
    return &tmp[0];
}

inline ordered_pair *
ordered_pair_vector_buf(std::vector<ordered_pair> *buf)
{
    std::vector<ordered_pair> &tmp = *buf;
    return &tmp[0];
}


static PyObject *
pickle_tree_buffer(std::vector<ckdtreenode> *buf)
{
    char *v = (char*) &(buf->front());
    Py_ssize_t len = buf->size() * sizeof(ckdtreenode);
    return ckdtree_PyBytes_FromStringAndSize(v,len);
}


static PyObject *
unpickle_tree_buffer(std::vector<ckdtreenode> *buf, PyObject *src)
{
    Py_ssize_t s, n;
    ckdtreenode *target, *cur;
    s = ckdtree_PyBytes_Size(src);
    if (PyErr_Occurred()) return NULL;
    n = s / sizeof(ckdtreenode);
    cur = (ckdtreenode *)ckdtree_PyBytes_AsString(src);   
    if (PyErr_Occurred()) return NULL;
    try {
        buf->resize(n);
        target = &(buf->front());
        std::memcpy((void*)target,(void*)cur,s);
    } catch (...) {
        translate_cpp_exception();
        return NULL;
    }
    Py_RETURN_NONE;
}
