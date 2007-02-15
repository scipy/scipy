#include <Python.h>
#include <structmember.h>
#include <arrayobject.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hli.h>
#include <chlilib.c>

//Constants
#define MAXOBJNAME 64
#define MAXNLLENGTH 1000

#define CALLFAME(cmd) Py_BEGIN_ALLOW_THREADS; cmd; Py_END_ALLOW_THREADS; if (checkError(status)) return NULL

#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/**********************************************************************/

static float  nmistt[3];        //Numeric
static double pmistt[3];        //Precision
static int    bmistt[3];        //Boolean
static int    dmistt[3];        //Date

//Numeric
static float N_ND = (float)1.701419e+038;
static float N_NC = (float)1.701418e+038;
static float N_NA = (float)1.701417e+038;

// Precision
static double P_ND = 1.70141507979e+038;
static double P_NC = 1.70141507978e+038;
static double P_NA = 1.70141507977e+038;
// Boolean
static int B_ND = 127;
static int B_NC = 126;
static int B_NA = 125;
// Date
static int D_ND = -1;
static int D_NC = -2;
static int D_NA = -3;

static char cfame_doc[] = "Module providing access to FAME functionality.";

//if there was an error, we need to set the Python error status before returning
static checkError(int status)
{
    if (status != HSUCC && status != HTRUNC)
    {
        char message[1000];
        PyErr_SetString(PyExc_RuntimeError, getsta(status, message));
        cfmfin(&status);
        return 1;
    }
    return 0;
}

static int makeTranslationTables(void)
{
    //Set up translation tables for ND, NC, NA mappings
    int status;

    cfmspm(&status, P_NC, P_ND, P_NA, pmistt);  //Precision
    cfmsnm(&status, N_NC, N_ND, N_NA, nmistt);  //Numeric
    cfmsbm(&status, B_NC, B_ND, B_NA, bmistt);  //Boolean
    cfmsdm(&status, D_NC, D_ND, D_NA, dmistt);  //Date
    return 0;
}

static char cfame_set_option_doc[] = "wrapper to cfmsopt";
static PyObject *
cfame_set_option(PyObject *self, PyObject *args)
{
    int status;
    char *name, *val;
    if (!PyArg_ParseTuple(args, "ss:set_option", &name, &val)) return NULL;

    CALLFAME(cfmsopt(&status, name, val));

    Py_RETURN_NONE;
}



static char cfame_open_doc[] = "C level open method. This is called from the __init__ method of FameDb in fame.py";
static PyObject *
cfame_open(PyObject *self, PyObject *args)
{
    int status;
    int dbkey, access;
    char *dbname;
    if (!PyArg_ParseTuple(args, "si:open", &dbname, &access)) return NULL;

    CALLFAME(cfmopdb (&status, &dbkey, dbname, access));

    return PyInt_FromLong(dbkey);
}

static char cfame_close_doc[] = "C level portion of the close method.";
static PyObject *
cfame_close(PyObject *self, PyObject *args)
{
    int status;
    int dbkey;
    if (!PyArg_ParseTuple(args, "i:close", &dbkey)) return NULL;

    cfmcldb (&status, dbkey);
    if (checkError(status)) return NULL;

    return PyInt_FromLong(0);
}

static char cfame_wildlist_doc[] = "C level portion of the wildlist method.";
static PyObject *
cfame_wildlist(PyObject *self, PyObject *args)
{
    int status;
    int dbkey;
    char *expression;
    int class, type, freq;
    char objnam[MAXOBJNAME+1];
    PyObject *result = PyList_New(0);

    if (!PyArg_ParseTuple(args, "is:wildlist(dbkey, expression)", &dbkey, &expression)) return NULL;

    // initialize wildlist
    CALLFAME(cfminwc (&status, dbkey, expression));

    // get first data object matching wildlist expression
    cfmnxwc (&status, dbkey, objnam, &class, &type, &freq);

    if (status == HNOOBJ)
        // no matching objects, return empty list
        return result;
    else
        if (checkError(status)) return NULL;

    while (status != HNOOBJ)
    {
        // append objnam to list
        if (PyList_Append(result, PyString_FromString(objnam))) return NULL;

        // get next item
        cfmnxwc (&status, dbkey, objnam, &class, &type, &freq);
        if (status != HNOOBJ && status != HSUCC)
            if (checkError(status)) return NULL;
    }

    return result;
}

// Make appropriate boolean mask for data (based on special constants)
static PyObject *make_mask(void *data, int arraylen, int type) {

    PyArrayObject *mask;
    int i;
    int *mask_raw;

    if ((mask_raw = malloc(arraylen * sizeof(int))) == NULL) return PyErr_NoMemory();

    switch(type)
    {

        case HNUMRC: { // numeric
            float *castdata = (float*)data;
            float val;

            for (i = 0; i < arraylen; i++) {
                val = castdata[i];
                if (val == N_ND || val == N_NC || val == N_NA) {
                    mask_raw[i] = 1;
                } else {
                    mask_raw[i] = 0;
                }
            }
        } break;
        case HBOOLN: { // boolean
            int *castdata = (int*)data;
            int val;

            for (i = 0; i < arraylen; i++) {
                val = castdata[i];
                if (val == B_ND || val == B_NC || val == B_NA) {
                    mask_raw[i] = 1;
                } else { mask_raw[i] = 0;}
            }
        } break;
        case HSTRNG: { // string
            char **castdata = (char**)data;
            char *val;
            for (i = 0; i < arraylen; i++) {
                val = castdata[i];
                if (val == "") {
                    mask_raw[i] = 1;
                } else { mask_raw[i] = 0; }
            }
        } break;
        case HPRECN: { // precision
            double *castdata = (double*)data;
            double val;
            for (i = 0; i < arraylen; i++) {
                val = castdata[i];
                if (val == P_ND || val == P_NC || val == P_NA) {
                    mask_raw[i] = 1;
                } else { mask_raw[i] = 0; }
            }
        } break;
        default:
            if (type >= 8) {
                int *castdata = (int*)data;
                int val;
                for (i = 0; i < arraylen; i++) {
                    val = castdata[i];
                    if (val == D_ND || val == D_NC || val == D_NA) {
                        mask_raw[i] = 1;
                    } else { mask_raw[i] = 0; }
                }
            } else {
                PyErr_SetString(PyExc_ValueError, "unsupported datatype");
                return NULL;
            }
    }

    mask = (PyArrayObject*)PyArray_SimpleNewFromData(1, &arraylen, PyArray_INT32, mask_raw);
    mask->flags = (mask->flags) | NPY_OWNDATA;

    return (PyObject*)mask;
}


static char cfame_read_doc[] = "C level portion of read method.";
static PyObject *
cfame_read(PyObject *self, PyObject *args)
{
    int status, dbkey, i;
    int dataFlag;

    char *object_name;

    int first_point, last_point; //this defines the custom range to read (-1 for both means read all)

    int max_string_len;

    int class, type, freq, start_year, start_period, end_year, end_period;  // data fields returned by cfmosiz
    int basis, observed, created_year, created_month, created_day, mod_year, mod_month, mod_day;  //additional fields for cfmwhat
    char desc[1], doc[1];
    PyObject * returnVal = NULL;
    PyObject * values = NULL;
    int numobjs, typeNum;
    int range[3];
    void* dbValues;

    desc[0] = 0x0;
    doc[0]  = 0x0;

    if (!PyArg_ParseTuple(args, "isiii:read",
                                &dbkey,
                                &object_name,
                                &first_point,
                                &last_point,
                                &max_string_len)) return NULL;   //get params

    CALLFAME(cfmwhat(&status, dbkey, object_name, &class, &type, &freq, &basis, &observed, &start_year, &start_period, &end_year, &end_period, &created_year, &created_month, &created_day, &mod_year, &mod_month, &mod_day, desc, doc));

    returnVal = PyDict_New();

    PyDict_SetItemString(returnVal, "type", PyInt_FromLong(type));
    PyDict_SetItemString(returnVal, "freq", PyInt_FromLong(freq));
    PyDict_SetItemString(returnVal, "class", PyInt_FromLong(class));
    PyDict_SetItemString(returnVal, "mod_year", PyInt_FromLong(mod_year));
    PyDict_SetItemString(returnVal, "mod_month", PyInt_FromLong(mod_month));
    PyDict_SetItemString(returnVal, "mod_day", PyInt_FromLong(mod_day));
    PyDict_SetItemString(returnVal, "observed", PyInt_FromLong(observed));
    PyDict_SetItemString(returnVal, "basis", PyInt_FromLong(basis));

    if (type == HNAMEL)     //namelists
    {
        int length;
        char names[MAXOBJNAME*MAXNLLENGTH+1];

        CALLFAME(cfmgtnl(&status, dbkey, object_name, HNLALL, names, MAXOBJNAME*MAXNLLENGTH, &length));
        PyDict_SetItemString(returnVal, "data", PyString_FromStringAndSize(names, length)); //just return the namelist as a comma delimited string
    }
    else
    {
        dataFlag = 1;

        switch (class)
        {
            case HSERIE:
                //initialize custom range

                //if we are dealing with a date we need to convert
                //'begin' and 'end' dates to year/period format

                if (first_point != -1) {
                    if (freq == HCASEX) {
                        start_year = 0;
                        start_period = first_point;
                    } else {
                        CALLFAME(cfmdatp(&status, freq, first_point, &start_year, &start_period));
                    }
                } else {
                    if (freq == HCASEX) {
                        /* for case series, if first_point not explicitly defined, always
                        read starting at index 1 (not the first data point like with
                        time series */
                        start_year = 0;
                        start_period = 1;
                    }
                }

                if (last_point != -1) {
                    if (freq == HCASEX) {
                        end_year = 0;
                        end_period = last_point;
                    } else {
                        CALLFAME(cfmdatp(&status, freq, last_point, &end_year, &end_period));
                    }
                }

                if (end_year < start_year ||
                    (start_year == end_year && end_period < start_period) ||
                    (start_period == -1)) {
                    dataFlag = 0;
                    break;
                }

                numobjs = -1;
                CALLFAME(cfmsrng(&status, freq, &start_year, &start_period, &end_year, &end_period, range, &numobjs)); //set the range of data to get
                break;

            case HSCALA:
                numobjs = 1;
                break;
            default:  //This should never happen
                PyErr_SetString(PyExc_RuntimeError, "Critical internal error #0 in CFAMEMODULE");
                return NULL;
        }

        if (dataFlag)
        {
            switch (type)   //initialize an array of the correct type to get the data from Fame
            {
                case HNUMRC:
                    if ((dbValues = malloc(numobjs * sizeof(float))) == NULL) return PyErr_NoMemory();
                    typeNum = PyArray_FLOAT;
                    break;
                case HPRECN:
                    if ((dbValues = malloc(numobjs * sizeof(double))) == NULL) return PyErr_NoMemory();
                    typeNum = PyArray_DOUBLE;
                    break;
                case HSTRNG:
                    typeNum = PyArray_OBJECT;
                    break;
                default:
                    if ((dbValues = malloc(numobjs * sizeof(int))) == NULL) return PyErr_NoMemory();
                    typeNum = PyArray_INT;
                    break;
            }
            if (type == HSTRNG)     //additional initilization for getting strings
            {

                if (class == HSERIE)
                {
                    PyObject** temp;

                    //string series
                    int* missing;
                    int* outlen;
                    int inlen[1];
                    int *mask_raw;

                    if ( ((dbValues = malloc(numobjs * sizeof(char*))) == NULL) ||
                         ((temp = malloc(numobjs * sizeof(PyObject*))) == NULL) ||
                         ((mask_raw = malloc(numobjs * sizeof(int))) == NULL) ||
                         ((missing = malloc(numobjs * sizeof(int))) == NULL) ||
                         ((outlen = malloc(numobjs * sizeof(int))) == NULL) ) {
                        return PyErr_NoMemory();
                    }

                    for (i = 0; i < numobjs; i++) {
                        if ((((char**)dbValues)[i] = malloc((max_string_len+1) * sizeof(char))) == NULL) {
                            return PyErr_NoMemory();
                        }
                    }

                    inlen[0] = -max_string_len;

                    //we need to know how big each string will be so that we can set up room for it
                    CALLFAME(cfmgtsts(&status, dbkey, object_name, range, dbValues, missing, inlen, outlen));
                    for (i = 0; i < numobjs; i++) {
                        if (outlen[i] > max_string_len) {
                            PyErr_SetString(PyExc_RuntimeError, "FAME returned a string longer than the max_string_len. Adjust max_string_len parameter.");
                            return NULL;
                        } else {

                            if (missing[i] != HNMVAL) {
                                if ((temp[i] = PyString_FromString("")) == NULL) {
                                    PyErr_SetString(PyExc_RuntimeError, "Failed to initialize missing string element.");
                                    return NULL;
                                }
                                mask_raw[i] = 1;
                            } else {
                                if ((temp[i] = PyString_FromStringAndSize(((char**)dbValues)[i], outlen[i])) == NULL) {
                                    return PyErr_NoMemory();
                                }
                                mask_raw[i] = 0;
                            }

                            free(((char**)dbValues)[i]);
                        }
                    }

                    free(dbValues);
                    dbValues = temp;

                    {
                        PyArrayObject* data = (PyArrayObject *)PyArray_SimpleNewFromData(1, &numobjs, typeNum, dbValues);
                        PyArrayObject* mask = (PyArrayObject*)PyArray_SimpleNewFromData(1, &numobjs, PyArray_INT32, mask_raw);
                        PyObject* startindex = PyInt_FromLong((long)range[1]);

                        // transfer ownership of dbValues to the array
                        data->flags = (data->flags) | NPY_OWNDATA;
                        mask->flags = (mask->flags) | NPY_OWNDATA;

                        PyDict_SetItemString(returnVal, "data", (PyObject*)data);
                        PyDict_SetItemString(returnVal, "mask", (PyObject*)mask);
                        PyDict_SetItemString(returnVal, "startindex", startindex);

                        Py_DECREF(data);
                        Py_DECREF(mask);
                        Py_DECREF(startindex);
                    }

                    free(missing);
                    free(outlen);
                }
                else
                {
                    //get one string
                    int missing;
                    int length;

                    if ((dbValues = malloc((max_string_len+1) * sizeof(char))) == NULL) return PyErr_NoMemory();

                    CALLFAME(cfmgtstr(&status, dbkey, object_name, NULL, dbValues, &missing, max_string_len, &length));

                    if (length > max_string_len) {
                        PyErr_SetString(PyExc_RuntimeError, "FAME returned a string longer than the maxlength. Use extra long string parameter");
                        return NULL;
                    }

                    {
                        PyObject* data = PyString_FromString((char*)dbValues);
                        PyObject* mask;
                        PyObject* startindex = PyInt_FromLong(-1);

                        if (missing != HNMVAL) { mask = PyBool_FromLong(1); }
                        else {                   mask = PyBool_FromLong(0); }

                        PyDict_SetItemString(returnVal, "data", data);
                        PyDict_SetItemString(returnVal, "mask", mask);
                        PyDict_SetItemString(returnVal, "startindex", startindex);

                        Py_DECREF(data);
                        Py_DECREF(mask);
                        Py_DECREF(startindex);
                    }

                }
            } else {
                switch(type)
                {

                    case HNUMRC:
                        CALLFAME(cfmrrng(&status, dbkey, object_name, range, dbValues, HTMIS, nmistt));
                        break;
                    case HBOOLN:
                        CALLFAME(cfmrrng(&status, dbkey, object_name, range, dbValues, HTMIS, bmistt));
                        break;
                    case HPRECN:
                        CALLFAME(cfmrrng(&status, dbkey, object_name, range, dbValues, HTMIS, pmistt));
                        break;
                    default:
                        if (type >= 8) {
                            CALLFAME(cfmrrng(&status, dbkey, object_name, range, dbValues, HTMIS, dmistt));
                        } else {
                            PyErr_SetString(PyExc_ValueError, "unsupported datatype");
                            return NULL;
                        }
                }

                {
                    PyArrayObject* data = (PyArrayObject *)PyArray_SimpleNewFromData(1, &numobjs, typeNum, dbValues);
                    PyObject* mask = make_mask(dbValues, numobjs, type);
                    PyObject* startindex = PyInt_FromLong((long)range[1]);

                    // transfer ownership of dbValues to the array
                    data->flags = (data->flags) | NPY_OWNDATA;

                    PyDict_SetItemString(returnVal, "data", (PyObject*)data);
                    PyDict_SetItemString(returnVal, "mask", mask);
                    PyDict_SetItemString(returnVal, "startindex", startindex);

                    Py_DECREF(data);
                    Py_DECREF(mask);
                    Py_DECREF(startindex);
                }
            }



        } // (dataFlag)


        //if dataFlag was set return an object with no data
        if (!dataFlag) {
            return returnVal;
        }
    } // (type == HNAMEL)     //namelists
    return returnVal;
}


// replace masked values with the ND constant
static PyArrayObject *replace_mask(PyObject *orig_data, PyObject *orig_mask, int type) {

    PyArrayObject *data_copy, *data, *mask;
    PyObject *valMask, *fillVal;
    int i;

    data_copy = (PyArrayObject *)PyArray_Copy((PyArrayObject *)orig_data);
    data = PyArray_GETCONTIGUOUS(data_copy);

    // don't care if mask is contiguous or not
    mask = (PyArrayObject *)orig_mask;

    switch(type)
    {

        case HNUMRC:
            fillVal = PyFloat_FromDouble(N_ND);
            break;
        case HBOOLN:
            fillVal = PyInt_FromLong(B_ND);
            break;
        case HSTRNG:
            fillVal = PyString_FromString("");
            break;
        case HPRECN:
            fillVal = PyFloat_FromDouble(P_ND);
            break;
        default:
            if (type >= 8) {
                fillVal = PyInt_FromLong(D_ND);
            } else {
                PyErr_SetString(PyExc_ValueError, "unsupported datatype");
                return NULL;
            }
    }

    for (i = 0; i < data->dimensions[0]; i++) {
        valMask = PyArray_GETITEM(mask, PyArray_GetPtr(mask, &i));
        if (PyInt_AsLong(valMask)) {
            PyArray_SETITEM(data, PyArray_GetPtr(data, &i), fillVal);
        }
        Py_DECREF(valMask);
    }
    return data;
}

static char cfame_write_series_doc[] =
    "C level portion of code for write_tser and write_cser method.";
static PyObject *
cfame_write_series(PyObject *self, PyObject *args)
{
    int status, dbkey;
    PyObject *dataArrayTemp, *maskArrayTemp;
    PyArrayObject *dataArray, *maskArray;
    char* name;
    char errMsg[500];
    int class, start_index, end_index, numobjs, source_type, type,
        source_freq, freq, start_year, start_period, end_year, end_period;
    PyObject * returnVal = NULL;
    int range[3];

    if (!PyArg_ParseTuple(args, "isOOiiii:write_series",
                                    &dbkey, &name,
                                    &dataArrayTemp, &maskArrayTemp,
                                    &start_index, &end_index,
                                    &source_type, &source_freq)) return NULL;

    CALLFAME(cfmosiz(&status, dbkey, name,
                     &class, &type, &freq, &start_year,
                     &start_period, &end_year, &end_period));

    if (source_type != type) {
        PyErr_SetString(PyExc_RuntimeError,
                        "received a non-matching type, cannot write");
        return NULL;
    }

    if (source_freq != freq) {
        PyErr_SetString(PyExc_RuntimeError,
                        "received a non-matching frequency, cannot write");
        return NULL;
    }

    numobjs = -1;
    if (freq == HCASEX) {
        start_year = 0;
        end_year = 0;
        start_period = start_index;
        end_period = end_index;
    } else if (freq >= 226) {  // HOURLY, MINUTELY, or SECONDLY
        CALLFAME(timeper(&status, freq, start_index, &start_year, &start_period));
        CALLFAME(timeper(&status, freq, end_index, &end_year, &end_period));
    } else {   //convert int dates to fame period dates
        CALLFAME(cfmdatp(&status, freq, start_index, &start_year, &start_period));
        CALLFAME(cfmdatp(&status, freq, end_index, &end_year, &end_period));
    }
    //set the range that we will be writing to
    CALLFAME(cfmsrng(&status, freq, &start_year, &start_period, &end_year, &end_period, range, &numobjs));
    if (!PyArray_Check(dataArrayTemp)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "write_series was passed something other than an ndarray");
        return NULL;
    }

    if (type == HSTRNG) {

        //setting strings requires a different function call
        int* missing;
        int* lengths;
        char** values;
        int i;

        PyObject *str, *mask;

        if (((missing = malloc(numobjs * sizeof(int))) == NULL) ||
            ((lengths = malloc(numobjs * sizeof(int))) == NULL) ||
            ((values = malloc(numobjs * sizeof(char*))) == NULL)) {
            return PyErr_NoMemory();
        }

        dataArray = (PyArrayObject*)dataArrayTemp;
        maskArray = (PyArrayObject*)maskArrayTemp;

        for (i = 0; i < numobjs; i++) {
            //extract a string and add it to the array to be written

            str = PyArray_GETITEM(dataArray, PyArray_GetPtr(dataArray, &i));
            mask = PyArray_GETITEM(maskArray, PyArray_GetPtr(maskArray, &i));

            lengths[i] = PyString_Size(str);
            if ((values[i] = malloc((lengths[i]+1) * sizeof(char))) == NULL) return PyErr_NoMemory();
            values[i] = PyString_AsString(str);

            if (PyInt_AsLong(mask)) {
                missing[i] = HNDVAL;
            } else {
                missing[i] = HNMVAL;
            }
        }
        //write all the strings to Fame
        CALLFAME(cfmwsts(&status, dbkey, name, range, values, missing, lengths));

        //clear the extra memory that the strings are using
        for (i = 0; i < numobjs; i++) {
            free(values[i]);
        }

        free(missing);
        free(lengths);
        free(values);

    } else {

        // replace masked values with the ND constant
        dataArray = replace_mask(dataArrayTemp, maskArrayTemp, type);

        switch (type) {
            case HNUMRC: { // numeric
                    CALLFAME(cfmwrng(&status, dbkey, name, range, dataArray->data, HTMIS, nmistt));
                    break;
                }
            case HPRECN: { // precision
                    CALLFAME(cfmwrng(&status, dbkey, name, range, dataArray->data, HTMIS, pmistt));
                    break;
                }
            case HBOOLN: { // boolean
                    CALLFAME(cfmwrng(&status, dbkey, name, range, dataArray->data, HTMIS, bmistt));
                    break;
                }
            default:
                if(type >= 8) { // date type
                    CALLFAME(cfmwrng(&status, dbkey, name, range, dataArray->data, HTMIS, dmistt));
                } else {
                    Py_DECREF(dataArray);
                    sprintf(errMsg, "unsupported data type: %i", type);
                    PyErr_SetString(PyExc_RuntimeError, errMsg);
                    return NULL;
                }

        }
        Py_DECREF(dataArray);
    }

    Py_RETURN_NONE;
}


static char cfame_write_scalar_doc[] = "C level portion of write_scalar method.";
static PyObject *
cfame_write_scalar(PyObject *self, PyObject *args)
{
    int status, dbkey;
    PyObject* object;
    char* name;
    int class, freq, start_year, start_period, end_year, end_period;  // data fields returned by cfmosiz
    PyObject * returnVal = NULL;
    int source_type, type;
    int range[3];

    if (!PyArg_ParseTuple(args, "isOi:write_scalar", &dbkey, &name, &object, &source_type)) return NULL;   //get params
    CALLFAME(cfmosiz(&status, dbkey, name, &class, &type, &freq, &start_year, &start_period, &end_year, &end_period));   //get object info

    if (source_type != type) {
        PyErr_SetString(PyExc_RuntimeError, "received a non-matching type, cannot write");
        return NULL;
    }

    switch (type) {
        case HSTRNG: {
                char* value;
                int length;

                length  = PyString_Size(object);
                value = malloc((length + 1) * sizeof(char));
                value = PyString_AsString(object);

                CALLFAME(cfmwstr(&status, dbkey, name, range, value, HNMVAL, length));
                free(value);
            } break;
        case HNUMRC: {
                float values[1];
                values[0] = (float)PyFloat_AsDouble(object);
                CALLFAME(cfmwrng(&status, dbkey, name, range, values, HNTMIS, NULL));
            } break;
        case HPRECN: {
                double values[1];
                values[0] = PyFloat_AsDouble(object);
                CALLFAME(cfmwrng(&status, dbkey, name, range, values, HNTMIS, NULL));
            } break;
        case HBOOLN: {
                int values[1];
                values[0] = (int)PyInt_AsLong(object);
                CALLFAME(cfmwrng(&status, dbkey, name, range, values, HNTMIS, NULL));
            } break;
        default:
            if (type >= 8) {
                // date data type
                int values[1];
                values[0] = (int)PyInt_AsLong(object);
                CALLFAME(cfmwrng(&status, dbkey, name, range, values, HNTMIS, NULL));
            } else {
                PyErr_SetString(PyExc_ValueError, "Unrecognized type, cannot write");
                return NULL;
            }
    }
    Py_RETURN_NONE;
}


static char cfame_write_namelist_doc[] = "C level portion of code for writing namelists.";
static PyObject *
cfame_write_namelist(PyObject *self, PyObject *args)
{
    int status, dbkey;
    char *name, *namelist;

    if (!PyArg_ParseTuple(args, "iss:writeNamelist", &dbkey, &name, &namelist)) return NULL;

    CALLFAME(cfmwtnl(&status, dbkey, name, HNLALL, namelist));

    Py_RETURN_NONE;
}

static char cfame_create_doc[] = "C level portion of code for creating objects.";
static PyObject *
cfame_create(PyObject *self, PyObject *args)
{
    int status, dbkey;
    char* object_name;
    int class_arg, freq_arg, type_arg, basis_arg, observed_arg;

    if (!PyArg_ParseTuple(args, "isiiiii:create",
                                &dbkey, &object_name, &class_arg,
                                &freq_arg, &type_arg, &basis_arg,
                                &observed_arg)) return NULL;

    CALLFAME(cfmnwob(&status, dbkey, object_name, class_arg, freq_arg, type_arg, basis_arg, observed_arg));

    Py_RETURN_NONE;
}

static char cfame_remove_doc[] = "C level portion of code for deleting objects.";
static PyObject*
cfame_remove(PyObject* self, PyObject* args)
{
    int status, dbkey;
    char* object_name;

    if (!PyArg_ParseTuple(args, "is:remove", &dbkey, &object_name)) return NULL;
    CALLFAME(cfmdlob(&status, dbkey, object_name));

    Py_RETURN_NONE;
}

static char cfame_exists_doc[] = "C level portion of code for checking existence of object.";
static PyObject*
cfame_exists(PyObject* self, PyObject* args)
{
    int status, dbkey;
    char* object_name;
    int deslen, doclen;

    if (!PyArg_ParseTuple(args, "is:exists", &dbkey, &object_name)) return NULL;

    cfmdlen (&status, dbkey, object_name, &deslen, &doclen);
    if (status == HNOOBJ)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

static char cfame_updated_doc[] = "C level portion of updated function.";
static PyObject*
cfame_updated(PyObject* self, PyObject* args)
{
    int status, dbkey;
    char *object_name;
    int class, type, freq, start_year, start_period, end_year, end_period;
    int basis, observ, created_year, created_month, created_day, mod_year, mod_month, mod_day;
    char desc[1], doc[1];
    PyObject * returnVal = NULL;

    desc[0] = 0x0;
    doc[0]  = 0x0;

    if (!PyArg_ParseTuple(args, "is:updated", &dbkey, &object_name)) return NULL;

    CALLFAME(cfmwhat(&status, dbkey, object_name, &class, &type, &freq, &basis, &observ,
                     &start_year, &start_period, &end_year, &end_period,
                     &created_year, &created_month, &created_day,
                     &mod_year, &mod_month, &mod_day,
                     desc, doc));

    returnVal = PyDict_New();

    PyDict_SetItemString(returnVal, "mod_year", PyInt_FromLong(mod_year));
    PyDict_SetItemString(returnVal, "mod_month", PyInt_FromLong(mod_month));
    PyDict_SetItemString(returnVal, "mod_day", PyInt_FromLong(mod_day));

    return returnVal;
}

static char cfame_whats_doc[] = "C level portion of whats function.";
static PyObject *
cfame_whats(PyObject *self, PyObject *args)
{
    //arguments
    char *object_name;
    int dbkey;

    //return val
    PyObject *returnVal = NULL;

    int deslen, doclen; /* Length of desc and doc string return from cfmdlen */

    int status;
    void *tempdesc, *tempdoc;
    char *desc, *doc;
    int class, type, freq, start_year, start_period, end_year, end_period;  // data fields returned by cfmosiz
    int basis, observ, created_year, created_month, created_day, mod_year, mod_month, mod_day;  //additional fields for cfmwhat

    PyObject *py_class, *py_type, *py_freq, *py_basis, *py_observ,
             *py_start_year, *py_start_period, *py_end_year, *py_end_period,
             *py_mod_year, *py_mod_month, *py_mod_day,
             *py_desc, *py_doc;

    if (!PyArg_ParseTuple(args, "is:whats",
                                &dbkey,
                                &object_name)) return NULL;

    /* Get the length of the desc and doc strings */
    CALLFAME(cfmdlen(&status, dbkey, object_name, &deslen, &doclen));

    /* allocate the memory needed for the desc/doc strings */
    if ((tempdesc = malloc((deslen + 1) * sizeof(char))) == NULL) return PyErr_NoMemory();
    if ((tempdoc = malloc((doclen + 1) * sizeof(char))) == NULL) return PyErr_NoMemory();

    /* set the memory to non-null chars to tell fame the length of string we will accept */
    memset(tempdesc, 'A', deslen);
    memset(tempdoc, 'A', doclen);

    /* cast to char array */
    desc = (char*)tempdesc;
    doc = (char*)tempdoc;

    /* terminate the string with a null */
    desc[deslen] = 0x0;
    doc[doclen] = 0x0;

    CALLFAME(cfmwhat(&status, dbkey, object_name, &class, &type, &freq, &basis, &observ,
                     &start_year, &start_period, &end_year, &end_period,
                     &created_year, &created_month, &created_day,
                     &mod_year, &mod_month, &mod_day,
                     desc, doc));

    py_class = PyInt_FromLong(class);
    py_type = PyInt_FromLong(type);
    py_freq = PyInt_FromLong(freq);
    py_basis = PyInt_FromLong(basis);
    py_observ = PyInt_FromLong(observ);
    py_start_year = PyInt_FromLong(start_year);
    py_start_period = PyInt_FromLong(start_period);
    py_end_year = PyInt_FromLong(end_year);
    py_end_period = PyInt_FromLong(end_period);
    py_mod_year = PyInt_FromLong(mod_year);
    py_mod_month = PyInt_FromLong(mod_month);
    py_mod_day = PyInt_FromLong(mod_day);

    py_desc = PyString_FromString(desc);
    py_doc = PyString_FromString(doc);

    returnVal = PyDict_New();

    PyDict_SetItemString(returnVal, "class", py_class);
    PyDict_SetItemString(returnVal, "type", py_type);
    PyDict_SetItemString(returnVal, "freq", py_freq);
    PyDict_SetItemString(returnVal, "basis", py_basis);
    PyDict_SetItemString(returnVal, "observ", py_observ);
    PyDict_SetItemString(returnVal, "start_year", py_start_year);
    PyDict_SetItemString(returnVal, "start_period", py_start_period);
    PyDict_SetItemString(returnVal, "end_year", py_end_year);
    PyDict_SetItemString(returnVal, "end_period", py_end_period);
    PyDict_SetItemString(returnVal, "mod_year", py_mod_year);
    PyDict_SetItemString(returnVal, "mod_month", py_mod_month);
    PyDict_SetItemString(returnVal, "mod_day", py_mod_day);

    PyDict_SetItemString(returnVal, "desc", py_desc);
    PyDict_SetItemString(returnVal, "doc", py_doc);

    free((void*)desc);
    free((void*)doc);

    Py_DECREF(py_class);
    Py_DECREF(py_type);
    Py_DECREF(py_freq);
    Py_DECREF(py_basis);
    Py_DECREF(py_observ);
    Py_DECREF(py_start_year);
    Py_DECREF(py_start_period);
    Py_DECREF(py_end_year);
    Py_DECREF(py_end_period);
    Py_DECREF(py_mod_year);
    Py_DECREF(py_mod_month);
    Py_DECREF(py_mod_day);

    Py_DECREF(py_desc);
    Py_DECREF(py_doc);

    return returnVal;
}

static char cfame_size_doc[] = "C level portion of size method.";
static PyObject *
cfame_size(PyObject *self, PyObject *args)
{
    //arguments
    char *object_name;
    int dbkey;

    //return val
    PyObject *returnVal = NULL;

    PyObject *py_class, *py_type, *py_freq, *py_start_year, *py_start_period,
             *py_end_year, *py_end_period;

    int status;
    int class, type, freq, start_year, start_period, end_year, end_period;  // data fields returned by cfmosiz

    if (!PyArg_ParseTuple(args, "is:size",
                                &dbkey,
                                &object_name)) return NULL;

    CALLFAME(cfmosiz(&status, dbkey, object_name, &class, &type, &freq, &start_year, &start_period, &end_year, &end_period));

    py_class = PyInt_FromLong(class);
    py_type = PyInt_FromLong(type);
    py_freq = PyInt_FromLong(freq);
    py_start_year = PyInt_FromLong(start_year);
    py_start_period = PyInt_FromLong(start_period);
    py_end_year = PyInt_FromLong(end_year);
    py_end_period = PyInt_FromLong(end_period);

    returnVal = PyDict_New();

    PyDict_SetItemString(returnVal, "class", py_class);
    PyDict_SetItemString(returnVal, "type", py_type);
    PyDict_SetItemString(returnVal, "freq", py_freq);
    PyDict_SetItemString(returnVal, "start_year", py_start_year);
    PyDict_SetItemString(returnVal, "start_period", py_start_period);
    PyDict_SetItemString(returnVal, "end_year", py_end_year);
    PyDict_SetItemString(returnVal, "end_period", py_end_period);

    Py_DECREF(py_class);
    Py_DECREF(py_type);
    Py_DECREF(py_freq);
    Py_DECREF(py_start_year);
    Py_DECREF(py_start_period);
    Py_DECREF(py_end_year);
    Py_DECREF(py_end_period);

    return returnVal;
}

static char cfame_restore_doc[] = "C level portion of restore method.";
static PyObject *
cfame_restore(PyObject *self, PyObject *args)
{
    int status, dbkey;

    if (!PyArg_ParseTuple(args, "i:restore", &dbkey)) return NULL;

    CALLFAME(cfmrsdb(&status, dbkey));
    Py_RETURN_NONE;
}

///////////////////////////////////////////////////////////////////////

static PyMethodDef cfame_methods[] = {
    {"open", cfame_open, METH_VARARGS, cfame_open_doc},
    {"close", cfame_close, METH_VARARGS, cfame_close_doc},
    {"wildlist", cfame_wildlist, METH_VARARGS, cfame_wildlist_doc},
    {"read", cfame_read, METH_VARARGS, cfame_read_doc},
    {"whats", cfame_whats, METH_VARARGS, cfame_whats_doc},
    {"size", cfame_size, METH_VARARGS, cfame_size_doc},
    {"set_option", cfame_set_option, METH_VARARGS, cfame_set_option_doc},
    {"write_scalar", cfame_write_scalar, METH_VARARGS, cfame_write_scalar_doc},
    {"write_series", cfame_write_series, METH_VARARGS, cfame_write_series_doc},
    {"create", cfame_create, METH_VARARGS, cfame_create_doc},
    {"remove", cfame_remove, METH_VARARGS, cfame_remove_doc},
    {"exists", cfame_exists, METH_VARARGS, cfame_exists_doc},
    {"updated", cfame_updated, METH_VARARGS, cfame_updated_doc},
    {"write_namelist", cfame_write_namelist, METH_VARARGS, cfame_write_namelist_doc},
    {"restore", cfame_restore, METH_VARARGS, cfame_restore_doc},
    {NULL, NULL}
};

PyMODINIT_FUNC
initcfame(void)
{
    int status;
    cfmini(&status);
    Py_InitModule3("cfame", cfame_methods, cfame_doc);
    import_array();

    makeTranslationTables();
}