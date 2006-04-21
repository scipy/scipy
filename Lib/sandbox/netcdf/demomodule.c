/*
 * Demonstration of using the netCDF module from another C module.
 *
 * Written by Konrad Hinsen
 * last revision: 2001-1-3
 */

#include "Python.h"
#include "numpy/arrayobject.h"
#include "netcdf/netcdfmodule.h"


/*
 * Create a file with two dimensions and one variable.
 *
 * Attention: there is no error checking in this code! All return
 * values ought to be tested.
 */
static PyObject *
create_file(self, args)
     PyObject *self; /* Not used */
     PyObject *args; /* Not used either */
{
  /* Pointer to file object */
  PyNetCDFFileObject *file;
  /* Dimension names for variable foo */
  char *dimensions[] = {"n", "xyz", "string_length"};
  /* Pointer to variable object */
  PyNetCDFVariableObject *foo, *bar;
  /* Pointer to indices */
  PyNetCDFIndex *indices;

  /* Open file */
  file = PyNetCDFFile_Open("demo.nc", "w");
  /* Add file attribute */
  PyNetCDFFile_SetAttributeString(file, "title", "useless data");
  /* Add history line */
  PyNetCDFFile_AddHistoryLine(file, "Created some day");
  /* Create two dimensions */
  PyNetCDFFile_CreateDimension(file, dimensions[0], 10);
  PyNetCDFFile_CreateDimension(file, dimensions[1], 3);
  PyNetCDFFile_CreateDimension(file, dimensions[2], 100);
  /* Create variable */
  foo = PyNetCDFFile_CreateVariable(file, "foo", 'l', dimensions, 2);
  /* Add variable attribute */
  PyNetCDFVariable_SetAttributeString(foo, "units", "arbitrary");
  /* Create index array */
  indices = PyNetCDFVariable_Indices(foo);
  /* Write zeros everywhere */
  PyNetCDFVariable_WriteArray(foo, indices, PyInt_FromLong(0));
  /* Create variable */
  bar = PyNetCDFFile_CreateVariable(file, "bar", 'c', dimensions+2, 1);
  /* Write string */
  PyNetCDFVariable_WriteString(bar,
                (PyStringObject *)PyString_FromString("nothing important"));
  /* Close file */
  PyNetCDFFile_Close(file);

  /* Return None */
  Py_INCREF(Py_None);
  return Py_None;
}

/* Table of functions defined in the module */

static PyMethodDef demo_methods[] = {
  {"createDemoFile",	create_file, 1},
  {NULL,		NULL}		/* sentinel */
};


/* Module initialization */

DL_EXPORT(void) initdemo(void)
{
  PyObject *module;
  PyObject *netcdf, *netcdf_dict;
  PyObject *c_api_object;

  /* Create the module and add the functions */
  module = Py_InitModule("demo", demo_methods);

  /* Import the array module */
  import_array();

  /* Import netcdf and retrieve its C API address array */
  netcdf = PyImport_ImportModule("Scientific.IO.NetCDF");
  if (netcdf != NULL) {
    netcdf_dict = PyModule_GetDict(netcdf);
    c_api_object = PyDict_GetItemString(netcdf_dict, "_C_API");
    if (PyCObject_Check(c_api_object)) {
      PyNetCDF_API = (void **)PyCObject_AsVoidPtr(c_api_object);
    }
  }

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module demo");
}
