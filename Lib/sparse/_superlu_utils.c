
#include "Python.h"
#include <setjmp.h>

jmp_buf _superlu_py_jmpbuf;
PyObject *_superlumodule_memory_dict=NULL;

void superlu_delete_allkeys() 
{
  PyObject *keys=NULL, *key=NULL;
  void *mem_ptr;
  int i;
 
  if (_superlumodule_memory_dict == NULL) {return;}
  keys = PyDict_Keys(_superlumodule_memory_dict);
  
  for (i = 0; i < PyList_Size(keys); i++) {
    key = PyList_GET_ITEM(keys, i);
    mem_ptr = (void *)PyInt_AS_LONG((PyIntObject *)key);
    free(mem_ptr);
    PyDict_DelItem(_superlumodule_memory_dict, key);
  }
  
  Py_XDECREF(keys);
}


/* Abort to be used inside the superlu module so that memory allocation 
   errors don't exit Python and memory allocated internal to SuperLU is freed.
   This should free all the pointers in the dictionary of allocated memory.
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

  mem_ptr = malloc(size);
  keyval = (long) mem_ptr;
  if (mem_ptr == NULL) return NULL;
  key = PyInt_FromLong(keyval);
  if (key == NULL) goto fail;
  if (PyDict_SetItem(_superlumodule_memory_dict, key, Py_None)) goto fail;
  Py_DECREF(key);
  return mem_ptr;

 fail:
  Py_XDECREF(key);
  free(mem_ptr);
  superlu_python_module_abort("superlu_malloc: Cannot set dictionary key value in malloc.");
 
}


void superlu_python_module_free(void *ptr)
{
  PyObject *key;
  long keyval;

  keyval = (long )ptr;
  key = PyInt_FromLong(keyval);
  PyDict_DelItem(_superlumodule_memory_dict, key);
  free(ptr);
  Py_DECREF(key);
  return; 
}

