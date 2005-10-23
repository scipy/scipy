module_name = 'spmatrix'
c_api_name = 'SpMatrix_API'
c_api_types = ['LLMatType', 'CSRMatType', 'SSSMatType']
c_api_functions = [('SpMatrix_ParseVecOpArgs', 'int', '(PyObject *args, double **x_data, double **y_data, int n)'),
                   ('SpMatrix_GetShape', 'int', '(PyObject *, int[])'),
                   ('SpMatrix_GetOrder', 'int', '(PyObject *, int*)'),
                   ('SpMatrix_GetItem', 'double', '(PyObject *, int, int)'),
                   ('SpMatrix_Matvec', 'int', '(PyObject *matrix, int nx, double *x, int ny, double *y)'),
                   ('SpMatrix_Precon', 'int', '(PyObject *repc, int n, double *x, double *y)'),
                   ('SpMatrix_NewLLMatObject', 'PyObject *', '(int dim[], int sym, int sizeHint)'),
                   ('SpMatrix_LLMatGetItem', 'double', '(LLMatObject *a, int i, int j)'),
                   ('SpMatrix_LLMatSetItem', 'int', '(LLMatObject *a, int i, int j, double x)'),
                   ('SpMatrix_LLMatUpdateItemAdd', 'int', '(LLMatObject *a, int i, int j, double x)'),
                   ('SpMatrix_LLMatBuildColIndex', 'int', '(struct llColIndex **idx, LLMatObject *self, int includeDiagonal)'),
                   ('SpMatrix_LLMatDestroyColIndex', 'void', '(struct llColIndex **idx)'),
                   ('ItSolvers_Solve', 'int', """(PyObject *linsolver, PyObject *A, int n, \\
                double *b, double *x, double tol, int itmax, PyObject *K, \\
                int *info, int *iter, double *relres)""")
                   ]

nof_types = len(c_api_types)
nof_functions = len(c_api_functions)

print """#ifndef %s_H
#define %s_H

/*
 * C API
 */

/* Type definitions */
""" % (c_api_name.upper(), c_api_name.upper())

for i in xrange(nof_types):
    print '#define %s_NUM %d' % (c_api_types[i], i)

print """
/* Function definitions */
"""

for i in xrange(nof_functions):
    print '#define %s_NUM %d' % (c_api_functions[i][0], nof_types + i)
    print '#define %s_RET %s' % (c_api_functions[i][0], c_api_functions[i][1])
    print '#define %s_PROTO %s' % (c_api_functions[i][0], c_api_functions[i][2])
    print

print """
/* Total number of C API pointers */"""

print "#define %s_pointers %d " % (c_api_name, nof_types + nof_functions)
print

print "#ifdef %s_MODULE" % (module_name.upper())
print

for f in c_api_functions:
    print 'static %s_RET %s %s_PROTO;' % (f[0], f[0], f[0])
print

print """#define init_c_api(MODULE_DICT) \\
{ \\
  static void *%s[%s_pointers]; \\
  PyObject *c_api_object; \\
\\
  /* initialize C API pointer array */ \\""" % (c_api_name, c_api_name)
for t in c_api_types:
    print '  %s[%s_NUM] = (void *) &%s; \\' % (c_api_name, t, t)
for f in c_api_functions:
    print '  %s[%s_NUM] = (void *) %s; \\' % (c_api_name, f[0], f[0])
print """\\
  /* Create a CObject containing the API pointer array s address */ \\
  c_api_object = PyCObject_FromVoidPtr((void *)%s, NULL); \\
\\
  /* create a name for this object in the module's namespace */ \\
  if (c_api_object != NULL) { \\
    PyDict_SetItemString((MODULE_DICT), "_C_API", c_api_object); \\
    Py_DECREF(c_api_object); \\
  } \\
}""" % (c_api_name, )

print """
#else

#ifdef %s_UNIQUE_SYMBOL
#define %s %s_UNIQUE_SYMBOL
#endif

/* C API address pointer */
#ifdef NO_IMPORT_%s
extern void **%s;
#else
#ifdef %s_UNIQUE_SYMBOL
void **%s;
#else
static void **%s;
#endif
#endif
"""% (module_name.upper(), c_api_name, module_name.upper(), module_name.upper(), c_api_name,
      module_name.upper(), c_api_name, c_api_name)

for t in c_api_types:
    print '#define %s *(PyTypeObject *)%s[%s_NUM]' % (t, c_api_name, t)
print

for f in c_api_functions:
    print """#define %s \\
  (*(%s_RET (*)%s_PROTO) \\
  %s[%s_NUM])
  """ % (f[0], f[0], f[0], c_api_name, f[0])

print """#define import_%s() \\
{ \\
  PyObject *%s = PyImport_ImportModule("%s"); \\
  if (%s != NULL) { \\
    PyObject *module_dict = PyModule_GetDict(%s); \\
    PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API"); \\
    if (PyCObject_Check(c_api_object)) { \\
      %s = (void **)PyCObject_AsVoidPtr(c_api_object); \\
    } \\
  } \\
  assert(%s != NULL); \\
}

#endif

#endif
""" % (module_name, module_name, module_name, module_name, module_name, c_api_name, c_api_name)
