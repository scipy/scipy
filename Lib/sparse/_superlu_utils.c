
#include "Python.h"
#include <setjmp.h>

jmp_buf _superlu_py_jmpbuf;
PyObject *_superlumodule_memory_dict=NULL;
PyObject *_superlumodule_newmemory_list=NULL;
int _superlumodule_flagnewmemory=0;

void superlu_flag_new_keys()
{
  _superlumodule_flagnewmemory = 1;
}

void superlu_end_new_keys()
{
  _superlumodule_flagnewmemory = 0;
  Py_XDECREF(_superlumodule_newmemory_list);
}

void superlu_delete_newkeys()
{
  PyObject *keys=NULL, *key=NULL;
  void *mem_ptr;
  int i, N;
 
  if (_superlumodule_memory_dict == NULL) {return;}
  if (_superlumodule_newmemory_list == NULL) {return;}
  keys = _superlumodule_newmemory_list;

  N = PyList_Size(keys);
  fprintf(stderr, "N = %d\n", N);
  fflush(stderr);
  for (i = 0; i < N; i++) {
    key = PyList_GET_ITEM(keys, i);
    mem_ptr = (void *)PyInt_AS_LONG((PyIntObject *)key);
    /* Only free if key not already removed from dictionary */
    if (!(PyDict_DelItem(_superlumodule_memory_dict, key))) {
        free(mem_ptr);
        fprintf(stderr, "Freeing %p\n", mem_ptr);
        fflush(stderr);
    }
  }
  
  fprintf(stderr, "Here..5\n");
  fflush(stderr);
  Py_XDECREF(keys);
  fprintf(stderr,"Here..7\n");
  fflush(stderr);
 
  _superlumodule_flagnewmemory=0;
}


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
   Calling program should free all the pointers in the dictionary of allocated memory.
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
  if ((_superlumodule_newmemory_list == NULL) & (_superlumodule_flagnewmemory)) {
    _superlumodule_newmemory_list = PyList_New(0);
  }
  fprintf(stderr, "Deep Inside One..\n");
  fflush(stderr);
  mem_ptr = malloc(size);
  keyval = (long) mem_ptr;
  if (mem_ptr == NULL) return NULL;
  key = PyInt_FromLong(keyval);
  if (key == NULL) goto fail;
  if (PyDict_SetItem(_superlumodule_memory_dict, key, Py_None)) goto fail;
  if (_superlumodule_flagnewmemory) {
    if (PyList_Append(_superlumodule_newmemory_list, key)) goto fail;
    fprintf(stderr, "Putting in %p\n", mem_ptr);
    fflush(stderr);
  }
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

  keyval = (long )ptr;
  key = PyInt_FromLong(keyval);
  if (PyDict_DelItem(_superlumodule_memory_dict, key)) {fprintf(stderr,"Problem... "); fflush(stderr);}
  fprintf(stderr, "Graceful Freeing %p\n", ptr);
  fflush(stderr);
  free(ptr);
  Py_DECREF(key);
  return; 
}

