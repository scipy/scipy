// below are the functions in Python's import.c file that 
// need changing to remove the global module lock.
// pretty much lock_import and unlock_import are the only
// places changes are needed.
// These changes were to Python 1.5.2.  I'm not sure if they'll
// work on Python 2.1

static PyObject *import_locks = NULL;
typedef struct {
 PyThread_type_lock lock;
 long thread;
 int level;
} mod_lock;

static void
lock_import(char* mod_name)
{
 PyObject* py_lock = NULL;
 mod_lock *lock = NULL;
 long current_thread = PyThread_get_thread_ident();

 if (current_thread == -1)
  return; /* Too bad */

 if (import_locks == NULL) {
     import_locks = PyDict_New();
  if( !import_locks)
         Py_FatalError("lock_import: failed to create lock list");
 }
    py_lock = PyDict_GetItemString(import_locks,mod_name);
    if (!py_lock)

        lock = malloc(sizeof(mod_lock));
        lock->lock = PyThread_allocate_lock();
        lock->thread = -1;
        lock->level = 0;
  py_lock = PyInt_FromLong((long)lock);
  if(!py_lock)
   Py_FatalError("lock_import: failed to create lock");
  PyDict_SetItemString(import_locks,mod_name,py_lock);
 }
 else
     lock = (mod_lock*)PyInt_AsLong(py_lock);
 if( lock->thread == current_thread ) {
     lock->level++;
     return; /* ??what do we need to do to clean up REFs */
 }
 if (lock->thread != 1 || !PyThread_acquire_lock(lock->lock,0)) {
    PyThreadState *tstate = PyEval_SaveThread();
  PyThread_acquire_lock(lock->lock, 1);
  PyEval_RestoreThread(tstate);
 }
    lock->thread = current_thread;
    lock->level = 1;
}

static void
unlock_import(char* mod_name)
{
 PyObject* py_lock = NULL;
 mod_lock *lock = NULL;
 long current_thread = PyThread_get_thread_ident();

 if (current_thread == -1)
  return; /* Too bad */
    if (!import_locks)
        Py_FatalError("unlock_import: no import lock list");

    py_lock = PyDict_GetItemString(import_locks,mod_name);
    if (!py_lock)
        Py_FatalError("unlock_import: module was never locked");

 lock = (mod_lock*)PyInt_AsLong(py_lock);
 if (lock->thread != current_thread)
  Py_FatalError("unlock_import: not holding the import lock");
 lock->level--;
 if (lock->level == 0) {
  lock->thread = -1; /* not necessary */
  PyThread_release_lock(lock->lock);
  PyDict_DelItemString(import_locks,mod_name);
  free(lock);
  if(!PyDict_Size(import_locks))
   py_DECREF(import_locks)
 }
}

PyObject *
PyImport_ImportModuleEx(name, globals, locals, fromlist)
 char *name;
 PyObject *globals;
 PyObject *locals;
 PyObject *fromlist;
{
 PyObject *result;
 lock_import(name);
 result = import_module_ex(name, globals, locals, fromlist);
 unlock_import(name);
 return result;
}
