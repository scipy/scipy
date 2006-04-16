#ifndef Py_NETCDFMODULE_H
#define Py_NETCDFMODULE_H
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Include file for netCDF files and variables.
 *
 * Written by Konrad Hinsen
 * last revision: 2001-5-4
 */


#include <stdio.h>

/* NetCDFFile object */

typedef struct {
  PyObject_HEAD
  PyObject *dimensions;   /* dictionary */
  PyObject *variables;    /* dictionary */
  PyObject *attributes;   /* dictionary */
  PyObject *name;         /* string */
  PyObject *mode;         /* string */
  int id;
  char open;
  char define;
  char write;
  int recdim;
} PyNetCDFFileObject;


/* NetCDFVariable object */

typedef struct {
  PyObject_HEAD
  PyNetCDFFileObject *file;
  PyObject *attributes;   /* dictionary */
  char *name;
  int *dimids;
  size_t *dimensions;
  int type;               /* same as array types */
  int nd;
  int id;
  char unlimited;
} PyNetCDFVariableObject;


/* Variable index structure */

typedef struct {
  int start;
  int stop;
  int stride;
  int item;
} PyNetCDFIndex;

/*
 * C API functions
 */

/* Type definitions */
#define PyNetCDFFile_Type_NUM 0
#define PyNetCDFVariable_Type_NUM 1

/* Open a netCDF file (i.e. create a new file object) */
#define PyNetCDFFile_Open_RET PyNetCDFFileObject *
#define PyNetCDFFile_Open_PROTO Py_PROTO((char *filename, char *mode))
#define PyNetCDFFile_Open_NUM 2

/* Close a netCDF file. Returns -1 if there was an error. */
#define PyNetCDFFile_Close_RET int
#define PyNetCDFFile_Close_PROTO Py_PROTO((PyNetCDFFileObject *file))
#define PyNetCDFFile_Close_NUM 3

/* Ensure that all data is written to the disk file.
   Returns 0 if there was an error. */
#define PyNetCDFFile_Sync_RET int
#define PyNetCDFFile_Sync_PROTO Py_PROTO((PyNetCDFFileObject *file))
#define PyNetCDFFile_Sync_NUM 4

/* Create a new dimension. Returns -1 if there was an error. */
#define PyNetCDFFile_CreateDimension_RET int
#define PyNetCDFFile_CreateDimension_PROTO \
        Py_PROTO((PyNetCDFFileObject *file, char *name, long size))
#define PyNetCDFFile_CreateDimension_NUM 5

/* Create a netCDF variable and return the variable object */
#define PyNetCDFFile_CreateVariable_RET PyNetCDFVariableObject *
#define PyNetCDFFile_CreateVariable_PROTO \
      Py_PROTO((PyNetCDFFileObject *file, char *name, int typecode, \
                char **dimension_names, int ndim))
#define PyNetCDFFile_CreateVariable_NUM 6

/* Return an object referring to an existing variable */
#define PyNetCDFFile_GetVariable_RET PyNetCDFVariableObject *
#define PyNetCDFFile_GetVariable_PROTO \
	  Py_PROTO((PyNetCDFFileObject *file, char *name))
#define PyNetCDFFile_GetVariable_NUM 7

/* Get variable rank */
#define PyNetCDFVariable_GetRank_RET int
#define PyNetCDFVariable_GetRank_PROTO Py_PROTO((PyNetCDFVariableObject *var))
#define PyNetCDFVariable_GetRank_NUM 8

/* Get variable shape */
#define PyNetCDFVariable_GetShape_RET size_t *
#define PyNetCDFVariable_GetShape_PROTO Py_PROTO((PyNetCDFVariableObject *var))
#define PyNetCDFVariable_GetShape_NUM 9

/* Allocate and initialize index structures for reading/writing data */
#define PyNetCDFVariable_Indices_RET PyNetCDFIndex *
#define PyNetCDFVariable_Indices_PROTO Py_PROTO((PyNetCDFVariableObject *var))
#define PyNetCDFVariable_Indices_NUM 10

/* Read data and return an array object */
#define PyNetCDFVariable_ReadAsArray_RET PyArrayObject *
#define PyNetCDFVariable_ReadAsArray_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var, PyNetCDFIndex *indices))
#define PyNetCDFVariable_ReadAsArray_NUM 11

/* Write array. Returns -1 if there was an error.  */
#define PyNetCDFVariable_WriteArray_RET int
#define PyNetCDFVariable_WriteArray_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var, PyNetCDFIndex *indices, \
		    PyObject *array))
#define PyNetCDFVariable_WriteArray_NUM 12

/* Get file attribute */
#define PyNetCDFFile_GetAttribute_RET PyObject *
#define PyNetCDFFile_GetAttribute_PROTO \
	  Py_PROTO((PyNetCDFFileObject *var, char *name))
#define PyNetCDFFile_GetAttribute_NUM 13

/* Set file attribute */
#define PyNetCDFFile_SetAttribute_RET int
#define PyNetCDFFile_SetAttribute_PROTO \
	  Py_PROTO((PyNetCDFFileObject *var, char *name, PyObject *value))
#define PyNetCDFFile_SetAttribute_NUM 14

/* Set file attribute to string value */
#define PyNetCDFFile_SetAttributeString_RET int
#define PyNetCDFFile_SetAttributeString_PROTO \
	  Py_PROTO((PyNetCDFFileObject *var, char *name, char *value))
#define PyNetCDFFile_SetAttributeString_NUM 15

/* Get variable attribute */
#define PyNetCDFVariable_GetAttribute_RET PyObject *
#define PyNetCDFVariable_GetAttribute_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var, char *name))
#define PyNetCDFVariable_GetAttribute_NUM 16

/* Set variable attribute */
#define PyNetCDFVariable_SetAttribute_RET int
#define PyNetCDFVariable_SetAttribute_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var, char *name, PyObject *value))
#define PyNetCDFVariable_SetAttribute_NUM 17

/* Set variable attribute to string value */
#define PyNetCDFVariable_SetAttributeString_RET int
#define PyNetCDFVariable_SetAttributeString_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var, char *name, char *value))
#define PyNetCDFVariable_SetAttributeString_NUM 18

/* Add entry to the history */
#define PyNetCDFFile_AddHistoryLine_RET int
#define PyNetCDFFile_AddHistoryLine_PROTO \
	  Py_PROTO((PyNetCDFFileObject *self, char *text))
#define PyNetCDFFile_AddHistoryLine_NUM 19

/* Write string. Returns -1 if there was an error.  */
#define PyNetCDFVariable_WriteString_RET int
#define PyNetCDFVariable_WriteString_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var, PyStringObject *value))
#define PyNetCDFVariable_WriteString_NUM 20

/* Read string  */
#define PyNetCDFVariable_ReadAsString_RET PyStringObject *
#define PyNetCDFVariable_ReadAsString_PROTO \
	  Py_PROTO((PyNetCDFVariableObject *var))
#define PyNetCDFVariable_ReadAsString_NUM 21

/* Total number of C API pointers */
#define PyNetCDF_API_pointers 22



#ifdef _NETCDF_MODULE

/* Type object declarations */
staticforward PyTypeObject PyNetCDFFile_Type;
staticforward PyTypeObject PyNetCDFVariable_Type;

/* Type check macros */
#define PyNetCDFFile_Check(op) ((op)->ob_type == &PyNetCDFFile_Type)
#define PyNetCDFVariable_Check(op) ((op)->ob_type == &PyNetCDFVariable_Type)

/* C API function declarations */
static PyNetCDFFile_Open_RET PyNetCDFFile_Open PyNetCDFFile_Open_PROTO;
static PyNetCDFFile_Close_RET PyNetCDFFile_Close PyNetCDFFile_Close_PROTO;
static PyNetCDFFile_Sync_RET PyNetCDFFile_Sync PyNetCDFFile_Sync_PROTO;
static PyNetCDFFile_CreateDimension_RET PyNetCDFFile_CreateDimension \
  PyNetCDFFile_CreateDimension_PROTO;
static PyNetCDFFile_CreateVariable_RET PyNetCDFFile_CreateVariable \
  PyNetCDFFile_CreateVariable_PROTO;
static PyNetCDFFile_GetVariable_RET PyNetCDFFile_GetVariable \
  PyNetCDFFile_GetVariable_PROTO;
static PyNetCDFVariable_GetRank_RET PyNetCDFVariable_GetRank \
  PyNetCDFVariable_GetRank_PROTO;
static PyNetCDFVariable_GetShape_RET PyNetCDFVariable_GetShape \
  PyNetCDFVariable_GetShape_PROTO;
static PyNetCDFVariable_Indices_RET PyNetCDFVariable_Indices \
  PyNetCDFVariable_Indices_PROTO;
static PyNetCDFVariable_ReadAsArray_RET PyNetCDFVariable_ReadAsArray \
  PyNetCDFVariable_ReadAsArray_PROTO;
static PyNetCDFVariable_ReadAsString_RET PyNetCDFVariable_ReadAsString \
  PyNetCDFVariable_ReadAsString_PROTO;
static PyNetCDFVariable_WriteArray_RET PyNetCDFVariable_WriteArray \
  PyNetCDFVariable_WriteArray_PROTO;
static PyNetCDFVariable_WriteString_RET PyNetCDFVariable_WriteString \
  PyNetCDFVariable_WriteString_PROTO;
static PyNetCDFFile_GetAttribute_RET PyNetCDFFile_GetAttribute \
  PyNetCDFFile_GetAttribute_PROTO;
static PyNetCDFFile_SetAttribute_RET PyNetCDFFile_SetAttribute \
  PyNetCDFFile_SetAttribute_PROTO;
static PyNetCDFFile_SetAttributeString_RET PyNetCDFFile_SetAttributeString \
  PyNetCDFFile_SetAttributeString_PROTO;
static PyNetCDFVariable_GetAttribute_RET PyNetCDFVariable_GetAttribute \
  PyNetCDFVariable_GetAttribute_PROTO;
static PyNetCDFVariable_SetAttribute_RET PyNetCDFVariable_SetAttribute \
  PyNetCDFVariable_SetAttribute_PROTO;
static PyNetCDFVariable_SetAttributeString_RET \
  PyNetCDFVariable_SetAttributeString \
  PyNetCDFVariable_SetAttributeString_PROTO;
static PyNetCDFFile_AddHistoryLine_RET PyNetCDFFile_AddHistoryLine \
  PyNetCDFFile_AddHistoryLine_PROTO;

#else

/* C API address pointer */ 
static void **PyNetCDF_API;

/* Type check macros */
#define PyNetCDFFile_Check(op) \
   ((op)->ob_type == (PyTypeObject *)PyNetCDF_API[PyNetCDFFile_Type_NUM])
#define PyNetCDFVariable_Check(op) \
   ((op)->ob_type == (PyTypeObject *)PyNetCDF_API[PyNetCDFVariable_Type_NUM])

/* C API function declarations */
#define PyNetCDFFile_Open \
  (*(PyNetCDFFile_Open_RET (*)PyNetCDFFile_Open_PROTO) \
   PyNetCDF_API[PyNetCDFFile_Open_NUM])
#define PyNetCDFFile_Close \
  (*(PyNetCDFFile_Close_RET (*)PyNetCDFFile_Close_PROTO) \
   PyNetCDF_API[PyNetCDFFile_Close_NUM])
#define PyNetCDFFile_Sync \
  (*(PyNetCDFFile_Sync_RET (*)PyNetCDFFile_Sync_PROTO) \
   PyNetCDF_API[PyNetCDFFile_Sync_NUM])
#define PyNetCDFFile_CreateDimension \
  (*(PyNetCDFFile_CreateDimension_RET (*)PyNetCDFFile_CreateDimension_PROTO) \
   PyNetCDF_API[PyNetCDFFile_CreateDimension_NUM])
#define PyNetCDFFile_CreateVariable \
  (*(PyNetCDFFile_CreateVariable_RET (*)PyNetCDFFile_CreateVariable_PROTO) \
   PyNetCDF_API[PyNetCDFFile_CreateVariable_NUM])
#define PyNetCDFFile_GetVariable \
  (*(PyNetCDFFile_GetVariable_RET (*)PyNetCDFFile_GetVariable_PROTO) \
   PyNetCDF_API[PyNetCDFFile_GetVariable_NUM])
#define PyNetCDFVariable_GetRank \
  (*(PyNetCDFVariable_GetRank_RET (*)PyNetCDFVariable_GetRank_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_GetRank_NUM])
#define PyNetCDFVariable_GetShape \
  (*(PyNetCDFVariable_GetShape_RET (*)PyNetCDFVariable_GetShape_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_GetShape_NUM])
#define PyNetCDFVariable_Indices \
  (*(PyNetCDFVariable_Indices_RET (*)PyNetCDFVariable_Indices_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_Indices_NUM])
#define PyNetCDFVariable_ReadAsArray \
  (*(PyNetCDFVariable_ReadAsArray_RET (*)PyNetCDFVariable_ReadAsArray_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_ReadAsArray_NUM])
#define PyNetCDFVariable_ReadAsString \
  (*(PyNetCDFVariable_ReadAsString_RET (*)PyNetCDFVariable_ReadAsString_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_ReadAsString_NUM])
#define PyNetCDFVariable_WriteArray \
  (*(PyNetCDFVariable_WriteArray_RET (*)PyNetCDFVariable_WriteArray_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_WriteArray_NUM])
#define PyNetCDFVariable_WriteString \
  (*(PyNetCDFVariable_WriteString_RET (*)PyNetCDFVariable_WriteString_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_WriteString_NUM])
#define PyNetCDFFile_GetAttribute \
  (*(PyNetCDFFile_GetAttribute_RET (*)PyNetCDFFile_GetAttribute_PROTO) \
   PyNetCDF_API[PyNetCDFFile_GetAttribute_NUM])
#define PyNetCDFFile_SetAttribute \
  (*(PyNetCDFFile_SetAttribute_RET (*)PyNetCDFFile_SetAttribute_PROTO) \
   PyNetCDF_API[PyNetCDFFile_SetAttribute_NUM])
#define PyNetCDFFile_SetAttributeString \
  (*(PyNetCDFFile_SetAttributeString_RET \
     (*)PyNetCDFFile_SetAttributeString_PROTO) \
   PyNetCDF_API[PyNetCDFFile_SetAttributeString_NUM])
#define PyNetCDFVariable_GetAttribute \
  (*(PyNetCDFVariable_GetAttribute_RET (*)PyNetCDFVariable_GetAttribute_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_GetAttribute_NUM])
#define PyNetCDFVariable_SetAttribute \
  (*(PyNetCDFVariable_SetAttribute_RET (*)PyNetCDFVariable_SetAttribute_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_SetAttribute_NUM])
#define PyNetCDFVariable_SetAttributeString \
  (*(PyNetCDFVariable_SetAttributeString_RET \
     (*)PyNetCDFVariable_SetAttributeString_PROTO) \
   PyNetCDF_API[PyNetCDFVariable_SetAttributeString_NUM])
#define PyNetCDFFile_AddHistoryLine \
  (*(PyNetCDFFile_AddHistoryLine_RET \
     (*)PyNetCDFFile_AddHistoryLine_PROTO) \
   PyNetCDF_API[PyNetCDFFile_AddHistoryLine_NUM])

#define import_netcdf() \
{ \
  PyObject *module = PyImport_ImportModule("Scientific.IO.NetCDF"); \
  if (module != NULL) { \
    PyObject *module_dict = PyModule_GetDict(module); \
    PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API"); \
    if (PyCObject_Check(c_api_object)) { \
      PyNetCDF_API = (void **)PyCObject_AsVoidPtr(c_api_object); \
    } \
  } \
}

#endif



#ifdef __cplusplus
}
#endif
#endif /* Py_NETCDFMODULE_H */
