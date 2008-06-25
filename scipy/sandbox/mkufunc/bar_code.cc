
import_ufunc();

return_val = PyUFunc_FromFuncAndData(
    foo_functions,
    foo_data,
    foo_signatures,
    2,         /* ntypes */
    1,            /* nin */
    1,                  /* nout */
    PyUFunc_None,       /* identity */
    "foo",              /* name */
    "",                 /* doc */
    0);
