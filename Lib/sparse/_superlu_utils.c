
#include "Python.h"
#include <setjmp.h>

jmp_buf _superlu_py_jmpbuf;
PyObject *_superlumodule_memory_dict=NULL;

/* Abort to be used inside the superlu module so that memory allocation 
   errors don't exit Python and memory allocated internal to SuperLU is freed.
   Calling program should deallocate (using SUPERLU_FREE) all memory that could have 
   been allocated.  (It's ok to FREE unallocated memory)---will be ignored.
*/

void superlu_python_module_abort(char *msg)
{
  PyErr_SetString(PyExc_RuntimeError, msg);
  longjmp(_superlu_py_jmpbuf, -1);
}

void *superlu_python_module_malloc(size_t size)
{
  PyObject *key=NULL;
  long keyval;
  void *mem_ptr; 

  if (_superlumodule_memory_dict == NULL) {
    _superlumodule_memory_dict = PyDict_New();
  }
  mem_ptr = malloc(size);
  if (mem_ptr == NULL) return NULL;
  keyval = (long) mem_ptr;
  key = PyInt_FromLong(keyval);
  if (key == NULL) goto fail;
  if (PyDict_SetItem(_superlumodule_memory_dict, key, Py_None)) goto fail;
  Py_DECREF(key);
  return mem_ptr;

 fail:
  Py_XDECREF(key);
  free(mem_ptr);
  superlu_python_module_abort("superlu_malloc: Cannot set dictionary key value in malloc.");
  return NULL;
 
}

void superlu_python_module_free(void *ptr)
{
  PyObject *key;
  long keyval;
  PyObject *ptype, *pvalue, *ptraceback;

  if (ptr == NULL) return;
  PyErr_Fetch(&ptype, &pvalue, &ptraceback);
  keyval = (long )ptr;
  key = PyInt_FromLong(keyval);
  /* This will only free the pointer if it could find it in the dictionary
     of already allocated pointers --- thus after abort, the module can free all
     the memory that "might" have been allocated to avoid memory leaks on abort 
     calls.
   */ 
  if (!(PyDict_DelItem(_superlumodule_memory_dict, key))) {
    free(ptr);
  }
  Py_DECREF(key);
  PyErr_Restore(ptype, pvalue, ptraceback);
  return; 
}

