#include "c_lib.h"
#include "c_tdates.h"
#include "c_tseries.h"


static PyMethodDef cseries_methods[] = {

    {"MA_mov_sum", (PyCFunction)MaskedArray_mov_sum,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_median", (PyCFunction)MaskedArray_mov_median,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_average", (PyCFunction)MaskedArray_mov_average,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"MA_mov_stddev", (PyCFunction)MaskedArray_mov_stddev,
     METH_VARARGS | METH_KEYWORDS, ""},

    {"TS_convert", (PyCFunction)TimeSeries_convert,
     METH_VARARGS, ""},

    {"DA_asfreq", (PyCFunction)DateArray_asfreq,
     METH_VARARGS, ""},
    {"DA_getDateInfo", (PyCFunction)DateArray_getDateInfo,
     METH_VARARGS, ""},

    {"thisday", (PyCFunction)c_tdates_thisday,
     METH_VARARGS, c_tdates_thisday_doc},
    {"check_freq", (PyCFunction)c_tdates_check_freq,
     METH_VARARGS, c_tdates_check_freq_doc},
    {"check_freq_str", (PyCFunction)c_tdates_check_freq_str,
     METH_VARARGS, c_tdates_check_freq_str_doc},
    {"get_freq_group", (PyCFunction)c_tdates_get_freq_group,
     METH_VARARGS, c_tdates_get_freq_group_doc},

    {"set_callback_DateFromString", (PyCFunction)set_callback_DateFromString,
     METH_VARARGS, ""},
    {"set_callback_DateTimeFromString", (PyCFunction)set_callback_DateTimeFromString,
     METH_VARARGS, ""},

    {NULL, NULL}
};

PyMODINIT_FUNC
initcseries(void)
{
    PyObject *m;

    m = Py_InitModule("cseries", cseries_methods);
    if (m == NULL)
      return;

    import_c_lib(m);
    import_c_tdates(m);
    import_c_tseries(m);
}
