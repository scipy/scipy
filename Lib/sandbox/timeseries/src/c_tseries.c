#include "c_tdates.h"
#include "c_tseries.h"

/* Helper function for TimeSeries_convert:
    determine the size of the second dimension for the resulting
    converted array */
static long get_height(int fromFreq, int toFreq) {

    int maxBusDaysPerYear, maxBusDaysPerQuarter, maxBusDaysPerMonth;
    int maxDaysPerYear, maxDaysPerQuarter, maxDaysPerMonth;

    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    if (fromGroup == FR_UND) { fromGroup = FR_DAY; }

    maxBusDaysPerYear = 262;
    maxBusDaysPerQuarter = 66;
    maxBusDaysPerMonth = 23;

    maxDaysPerYear = 366;
    maxDaysPerQuarter = 92;
    maxDaysPerMonth = 31;

    switch(fromGroup)
    {
        case FR_ANN: return 1;
        case FR_QTR:
            switch(toGroup)
            {
                case FR_ANN: return 4;
                default: return 1;
            }
        case FR_MTH: //monthly
            switch(toGroup)
            {
                case FR_ANN: return 12;
                case FR_QTR: return 3;
                default: return 1;
            }
        case FR_WK: //weekly
            switch(toGroup)
            {
                case FR_ANN: return 53;
                case FR_QTR: return 13;
                case FR_MTH: return 4;
                default: return 1;
            }
        case FR_BUS: //business
            switch(toGroup)
            {
                case FR_ANN: return maxBusDaysPerYear;;
                case FR_QTR: return maxBusDaysPerQuarter;
                case FR_MTH: return maxBusDaysPerMonth;
                case FR_WK: return 5;
                default: return 1;
            }
        case FR_DAY: //daily
            switch(toGroup)
            {
                case FR_ANN: return maxDaysPerYear;;
                case FR_QTR: return maxDaysPerQuarter;
                case FR_MTH: return maxDaysPerMonth;
                case FR_WK: return 7;
                default: return 1;
            }
        case FR_HR: //hourly
            switch(toGroup)
            {
                case FR_ANN: return 24 * maxDaysPerYear;;
                case FR_QTR: return 24 * maxDaysPerQuarter;
                case FR_MTH: return 24 * maxDaysPerMonth;
                case FR_WK: return 24 * 7;
                case FR_DAY: return 24;
                case FR_BUS: return 24;
                default: return 1;
            }
        case FR_MIN: //minutely
            switch(toGroup)
            {
                case FR_ANN: return 24 * 60 * maxDaysPerYear;;
                case FR_QTR: return 24 * 60 * maxDaysPerQuarter;
                case FR_MTH: return 24 * 60 * maxDaysPerMonth;
                case FR_WK: return 24 * 60 * 7;
                case FR_DAY: return 24 * 60;
                case FR_BUS: return 24 * 60;
                case FR_HR: return 60;
                default: return 1;
            }
        case FR_SEC: //minutely
            switch(toGroup)
            {
                case FR_ANN: return 24 * 60 * 60 * maxDaysPerYear;;
                case FR_QTR: return 24 * 60 * 60 * maxDaysPerQuarter;
                case FR_MTH: return 24 * 60 * 60 * maxDaysPerMonth;
                case FR_WK: return 24 * 60 * 60 * 7;
                case FR_DAY: return 24 * 60 * 60;
                case FR_BUS: return 24 * 60 * 60;
                case FR_HR: return 60 * 60;
                case FR_MIN: return 60;
                default: return 1;
            }
        default: return 1;
    }
}

PyObject *
TimeSeries_convert(PyObject *self, PyObject *args)
{
    PyObject *arrayTest;
    PyArrayObject *array, *newArray;
    PyArrayObject *mask, *newMask;

    PyObject *returnVal = NULL;
    PyObject *start_index_retval;

    long startIndex;
    long newStart, newStartTemp;
    long newEnd, newEndTemp;
    long newLen, newHeight;
    int i;
    long currIndex, prevIndex;
    long nd;
    npy_intp *dim, *newIdx;
    long currPerLen;
    char *position;
    PyObject *fromFreq_arg, *toFreq_arg;
    int fromFreq, toFreq;
    char relation;
    asfreq_info af_info;

    PyObject *val, *valMask;

    long (*asfreq_main)(long, char, asfreq_info*) = NULL;
    long (*asfreq_endpoints)(long, char, asfreq_info*) = NULL;
    long (*asfreq_reverse)(long, char, asfreq_info*) = NULL;

    returnVal = PyDict_New();

    if (!PyArg_ParseTuple(args,
        "OOOslO:convert(array, fromfreq, tofreq, position, startIndex, mask)",
        &array, &fromFreq_arg, &toFreq_arg,
        &position, &startIndex, &mask)) return NULL;

    if((fromFreq = check_freq(fromFreq_arg)) == INT_ERR_CODE) return NULL;
    if((toFreq = check_freq(toFreq_arg)) == INT_ERR_CODE) return NULL;

    if (toFreq == fromFreq)
    {
        PyObject *sidx;
        newArray = (PyArrayObject *)PyArray_Copy(array);
        newMask = (PyArrayObject *)PyArray_Copy(mask);
        sidx = PyInt_FromLong(startIndex);

        PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
        PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);
        PyDict_SetItemString(returnVal, "startindex", sidx);

        Py_DECREF(newArray);
        Py_DECREF(newMask);
        Py_DECREF(sidx);

        return returnVal;
    }

    switch(position[0])
    {
        case 'S':
            // start -> before
            relation = 'B';
            break;
        case 'E':
            // end -> after
            relation = 'A';
            break;
        default:
            return NULL;
            break;
    }

    get_asfreq_info(fromFreq, toFreq, &af_info);

    asfreq_main = get_asfreq_func(fromFreq, toFreq, 1);
    asfreq_endpoints = get_asfreq_func(fromFreq, toFreq, 0);

    //convert start index to new frequency
    CHECK_ASFREQ(newStartTemp = asfreq_main(startIndex, 'B', &af_info));
    if (newStartTemp < 1) {
        CHECK_ASFREQ(newStart = asfreq_endpoints(startIndex, 'A', &af_info));
    }
    else { newStart = newStartTemp; }

    //convert end index to new frequency
    CHECK_ASFREQ(newEndTemp = asfreq_main(startIndex+array->dimensions[0]-1, 'A', &af_info));
    if (newEndTemp < 1) {
        CHECK_ASFREQ(newEnd = asfreq_endpoints(startIndex+array->dimensions[0]-1, 'B', &af_info));
    }
    else { newEnd = newEndTemp; }

    if (newStart < 1) {
        PyErr_SetString(PyExc_ValueError, "start_date outside allowable range for destination frequency");
        return NULL;
    }

    newLen = newEnd - newStart + 1;
    newHeight = get_height(fromFreq, toFreq);

    if (newHeight > 1) {
        long tempval;
        asfreq_info af_info_rev;

        get_asfreq_info(toFreq, fromFreq, &af_info_rev);
        asfreq_reverse = get_asfreq_func(toFreq, fromFreq, 0);

        CHECK_ASFREQ(tempval = asfreq_reverse(newStart, 'B', &af_info_rev));
        currPerLen = startIndex - tempval;

        nd = 2;
        dim = PyDimMem_NEW(nd);
        dim[0] = (npy_intp)newLen;
        dim[1] = (npy_intp)newHeight;
    } else {
        nd = 1;
        dim = PyDimMem_NEW(nd);
        dim[0] = (npy_intp)newLen;
    }

    newIdx = PyDimMem_NEW(nd);
    arrayTest = PyArray_SimpleNew(nd, dim, array->descr->type_num);
    if (arrayTest == NULL) { return NULL; }
    newArray = (PyArrayObject*)arrayTest;
    newMask  = (PyArrayObject*)PyArray_SimpleNew(nd, dim, mask->descr->type_num);

    PyDimMem_FREE(dim);

    PyArray_FILLWBYTE(newArray,0);
    PyArray_FILLWBYTE(newMask,1);

    prevIndex = newStart;

    //set values in the new array
    for (i = 0; i < array->dimensions[0]; i++) {

        val = PyArray_GETITEM(array, PyArray_GetPtr(array, &i));
        valMask = PyArray_GETITEM(mask, PyArray_GetPtr(mask, &i));

        CHECK_ASFREQ(currIndex = asfreq_main(startIndex + i, relation, &af_info));

        newIdx[0] = currIndex-newStart;

        if (newHeight > 1) {

                if (currIndex != prevIndex)
                {
                    //reset period length
                    currPerLen = 0;
                    prevIndex = currIndex;
                }

                newIdx[1] = currPerLen;
                currPerLen++;
        }

        if (newIdx[0] > -1) {
            PyArray_SETITEM(newArray, PyArray_GetPtr(newArray, newIdx), val);
            PyArray_SETITEM(newMask, PyArray_GetPtr(newMask, newIdx), valMask);
        }

        Py_DECREF(val);
        Py_DECREF(valMask);

    }

    PyDimMem_FREE(newIdx);

    start_index_retval = (PyObject*)PyInt_FromLong(newStart);

    PyDict_SetItemString(returnVal, "values", (PyObject*)newArray);
    PyDict_SetItemString(returnVal, "mask", (PyObject*)newMask);
    PyDict_SetItemString(returnVal, "startindex", start_index_retval);

    Py_DECREF(newArray);
    Py_DECREF(newMask);
    Py_DECREF(start_index_retval);

    return returnVal;
}


/* This function is directly copied from direct copy of function in  */
/* Return typenumber from dtype2 unless it is NULL, then return
   NPY_DOUBLE if dtype1->type_num is integer or bool
   and dtype1->type_num otherwise.
*/
static int
_get_type_num_double(PyArray_Descr *dtype1, PyArray_Descr *dtype2)
{
        if (dtype2 != NULL)
                return dtype2->type_num;

        /* For integer or bool data-types */
        if (dtype1->type_num < NPY_FLOAT) {
                return NPY_DOUBLE;
        }
        else {
                return dtype1->type_num;
        }
}

#define _CHKTYPENUM(typ) ((typ) ? (typ)->type_num : PyArray_NOTYPE)

/* validates the standard arguments to moving functions and set the original
   mask, original ndarray, and mask for the result */
static PyObject *
check_mov_args(PyObject *orig_arrayobj, int span, int min_win_size,
               PyObject **orig_ndarray, PyObject **result_mask) {

    PyObject *orig_mask=NULL;
    PyArrayObject **orig_ndarray_tmp, **result_mask_tmp;
    int *raw_result_mask;

    if (!PyArray_Check(orig_arrayobj)) {
        PyErr_SetString(PyExc_ValueError, "array must be a valid subtype of ndarray");
        return NULL;
    }

    // check if array has a mask, and if that mask is an array
    if (PyObject_HasAttrString(orig_arrayobj, "_mask")) {
        PyObject *tempMask = PyObject_GetAttrString(orig_arrayobj, "_mask");
        if (PyArray_Check(tempMask)) {
            orig_mask = PyArray_EnsureArray(tempMask);
        } else {
            Py_DECREF(tempMask);
        }
    }

    *orig_ndarray = PyArray_EnsureArray(orig_arrayobj);
    orig_ndarray_tmp = (PyArrayObject**)orig_ndarray;

    if ((*orig_ndarray_tmp)->nd != 1) {
        PyErr_SetString(PyExc_ValueError, "array must be 1 dimensional");
        return NULL;
    }

    if (span < min_win_size) {
        char *error_str;
        error_str = malloc(60 * sizeof(char));
        MEM_CHECK(error_str)
        sprintf(error_str,
                "span must be greater than or equal to %i",
                min_win_size);
        PyErr_SetString(PyExc_ValueError, error_str);
        free(error_str);
        return NULL;
    }

    raw_result_mask = malloc((*orig_ndarray_tmp)->dimensions[0] * sizeof(int));
    MEM_CHECK(raw_result_mask)

    {
        PyArrayObject *orig_mask_tmp;
        int i, valid_points=0, is_masked;

        orig_mask_tmp = (PyArrayObject*)orig_mask;

        for (i=0; i<((*orig_ndarray_tmp)->dimensions[0]); i++) {

            is_masked=0;

            if (orig_mask != NULL) {
                PyObject *valMask;
                valMask = PyArray_GETITEM(orig_mask_tmp,
                                          PyArray_GetPtr(orig_mask_tmp, &i));
                is_masked = (int)PyInt_AsLong(valMask);
                Py_DECREF(valMask);
            }

            if (is_masked) {
                valid_points=0;
            } else {
                if (valid_points < span) { valid_points += 1; }
                if (valid_points < span) { is_masked = 1; }
            }

            raw_result_mask[i] = is_masked;
        }
    }

    *result_mask = PyArray_SimpleNewFromData(
                             1, (*orig_ndarray_tmp)->dimensions,
                             PyArray_INT32, raw_result_mask);
    MEM_CHECK(*result_mask)
    result_mask_tmp = (PyArrayObject**)result_mask;
    (*result_mask_tmp)->flags = ((*result_mask_tmp)->flags) | NPY_OWNDATA;
    return 0;
}

/* computation portion of moving sum. Appropriate mask is overlayed on top
   afterwards */
static PyObject*
calc_mov_sum(PyArrayObject *orig_ndarray, int span, int rtype)
{
    PyArrayObject *result_ndarray=NULL;
    int i;

    result_ndarray = (PyArrayObject*)PyArray_ZEROS(
                                       orig_ndarray->nd,
                                       orig_ndarray->dimensions,
                                       rtype, 0);
    ERR_CHECK(result_ndarray)

    for (i=0; i<orig_ndarray->dimensions[0]; i++) {

        PyObject *val=NULL, *mov_sum_val=NULL;

        val = PyArray_GETITEM(orig_ndarray, PyArray_GetPtr(orig_ndarray, &i));

        if (i == 0) {
            mov_sum_val = val;
        } else {
            int prev_idx = i-1;
            PyObject *mov_sum_prevval;
            mov_sum_prevval= PyArray_GETITEM(result_ndarray,
                                   PyArray_GetPtr(result_ndarray, &prev_idx));
            mov_sum_val = np_add(val, mov_sum_prevval);
            Py_DECREF(mov_sum_prevval);
            ERR_CHECK(mov_sum_val)

            if (i >= span) {
                PyObject *temp_val, *rem_val;
                int rem_idx = i-span;
                temp_val = mov_sum_val;
                rem_val = PyArray_GETITEM(orig_ndarray,
                                   PyArray_GetPtr(orig_ndarray, &rem_idx));

                mov_sum_val = np_subtract(temp_val, rem_val);
                ERR_CHECK(mov_sum_val)

                Py_DECREF(temp_val);
                Py_DECREF(rem_val);
            }
        }

        PyArray_SETITEM(result_ndarray,
                        PyArray_GetPtr(result_ndarray, &i),
                        mov_sum_val);

        if (mov_sum_val != val) { Py_DECREF(val); }

        Py_DECREF(mov_sum_val);
    }

    return (PyObject*)result_ndarray;

}

PyObject *
MaskedArray_mov_sum(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL,
             *result_ndarray=NULL, *result_mask=NULL,
             *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_sum(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;

    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &result_mask);

    rtype = _CHKTYPENUM(dtype);

    result_ndarray = calc_mov_sum((PyArrayObject*)orig_ndarray,
                                  span, rtype);
    ERR_CHECK(result_ndarray)

    result_dict = PyDict_New();
    MEM_CHECK(result_dict)
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

PyObject *
MaskedArray_mov_average(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL,
             *result_ndarray=NULL, *result_mask=NULL,
             *result_dict=NULL,
             *mov_sum=NULL, *denom=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_average(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;


    check_mov_args(orig_arrayobj, span, 2,
                   &orig_ndarray, &result_mask);

    rtype = _get_type_num_double(((PyArrayObject*)orig_ndarray)->descr, dtype);

    mov_sum = calc_mov_sum((PyArrayObject*)orig_ndarray, span, rtype);
    ERR_CHECK(mov_sum)

    denom = PyFloat_FromDouble(1.0/(double)(span));

    result_ndarray = np_multiply(mov_sum, denom);
    ERR_CHECK(result_ndarray)

    Py_DECREF(mov_sum);
    Py_DECREF(denom);

    result_dict = PyDict_New();
    MEM_CHECK(result_dict)
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}


/* computation portion of moving median. Appropriate mask is overlayed on top
   afterwards.

   The algorithm used here is based on the code found at:
    http://cran.r-project.org/src/contrib/Devel/runStat_1.1.tar.gz

   This code was originally released under the GPL, but the author
   (David Brahm) has granted me (and scipy) permission to use it under the BSD
   license. */
PyObject*
calc_mov_median(PyArrayObject *orig_ndarray, int span, int rtype)
{
    PyArrayObject *result_ndarray=NULL;
    PyObject **result_array, **ref_array, **even_array=NULL;
    PyObject *new_val, *old_val;
    PyObject *temp_add, *one_half;
    int a, i, k, R, arr_size, z;
    int *r;

    arr_size = orig_ndarray->dimensions[0];

    result_ndarray = (PyArrayObject*)PyArray_ZEROS(
                                       orig_ndarray->nd,
                                       orig_ndarray->dimensions,
                                       rtype, 0);
    ERR_CHECK(result_ndarray)

    if (arr_size >= span) {
        result_array = calloc(arr_size, sizeof(PyObject*));
        MEM_CHECK(result_array)

        /* this array will be used for quick access to the data in the original
           array (so PyArray_GETITEM doesn't have to be used over and over in the
           main loop) */
        ref_array = malloc(arr_size * sizeof(PyObject*));
        MEM_CHECK(ref_array)

        for (i=0; i<arr_size; i++) {
            ref_array[i] = PyArray_GETITEM(orig_ndarray, PyArray_GetPtr(orig_ndarray, &i));
        }

        /* this array wll be used for keeping track of the "ranks" of the values
           in the current window */
        r = malloc(span * sizeof(int));
        MEM_CHECK(r)

        for (i=0; i < span; i++) {
            r[i] = 1;
        }

        if ((span % 2) == 0) {
            // array to store two median values when span is an even #
            even_array = calloc(2, sizeof(PyObject*));
            MEM_CHECK(even_array)
        }

        R = (span + 1)/2;
        one_half = PyFloat_FromDouble(0.5);

        z = arr_size - span;

        /* Calculate initial ranks "r" */
        for (i=0; i < span; i++) {

            for (k=0;   k < i;  k++) {
                if (np_greater_equal(ref_array[z+i], ref_array[z+k])) {
                    r[i]++;
                }
            }
            for (k=i+1; k < span; k++) {
                if (np_greater(ref_array[z+i], ref_array[z+k])) {
                    r[i]++;
                }
            }

            /* If rank=R, this is the median */
            if (even_array != NULL) {
                if (r[i]==R) {
                    even_array[0] = ref_array[z+i];
                } else if (r[i] == (R+1)) {
                    even_array[1] = ref_array[z+i];
                }
            } else {
                if (r[i]==R) {
                    result_array[arr_size-1] = ref_array[z+i];
                }
            }
        }

        if (even_array != NULL) {
            temp_add = np_add(even_array[0], even_array[1]);
            result_array[arr_size-1] = np_multiply(temp_add, one_half);
            Py_DECREF(temp_add);
        }

        for (i=arr_size-2; i >= span-1; i--) {
            a = span;
            z = i - span + 1;
            old_val = ref_array[i+1];
            new_val = ref_array[i-span+1];

            for (k=span-1; k > 0; k--) {
                r[k] = r[k-1]; /* Shift previous iteration's ranks */
                if (np_greater_equal(ref_array[z+k], new_val)) {r[k]++; a--;}
                if (np_greater(ref_array[z+k], old_val)) {r[k]--;}

                if (r[k]==R) {
                    result_array[i] = ref_array[z+k];
                }

                if (even_array != NULL) {
                    if (r[k]==R) {
                        even_array[0] = ref_array[z+k];
                    } else if (r[k] == (R+1)) {
                        even_array[1] = ref_array[z+k];
                    }
                } else {
                    if (r[k]==R) {
                        result_array[i] = ref_array[z+k];
                    }
                }

            }

            r[0] = a;

            if (even_array != NULL) {
                if (a==R) {
                    even_array[0] = new_val;
                } else if (a == (R+1)) {
                    even_array[1] = new_val;
                }

                temp_add = np_add(even_array[0], even_array[1]);
                result_array[i] = np_multiply(temp_add, one_half);;
                Py_DECREF(temp_add);

            } else {
                if (a==R) {
                    result_array[i] = new_val;
                }
            }

        }

        Py_DECREF(one_half);

        for (i=span-1; i<arr_size; i++) {
            PyArray_SETITEM(result_ndarray,
                            PyArray_GetPtr(result_ndarray, &i),
                            result_array[i]);
        }

        for (i=0; i<arr_size; i++) {
            Py_DECREF(ref_array[i]);
        }

        if (even_array != NULL) {
            for (i=span-1; i<arr_size; i++) {
                Py_DECREF(result_array[i]);
            }
        }

        free(ref_array);
        free(result_array);
    }

    return (PyObject*)result_ndarray;

}

PyObject *
MaskedArray_mov_median(PyObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *orig_arrayobj=NULL, *orig_ndarray=NULL,
             *result_ndarray=NULL, *result_mask=NULL, *result_dict=NULL;
    PyArray_Descr *dtype=NULL;

    int rtype, span;

    static char *kwlist[] = {"array", "span", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oi|O&:mov_median(array, span, dtype)", kwlist,
                &orig_arrayobj, &span,
                PyArray_DescrConverter2, &dtype)) return NULL;

    check_mov_args(orig_arrayobj, span, 1,
                   &orig_ndarray, &result_mask);

    rtype = _CHKTYPENUM(dtype);

    result_ndarray = calc_mov_median((PyArrayObject*)orig_ndarray,
                                     span, rtype);
    ERR_CHECK(result_ndarray)

    result_dict = PyDict_New();
    MEM_CHECK(result_dict)
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}


PyObject *
MaskedArray_mov_stddev(PyObject *self, PyObject *args, PyObject *kwds)
{

    PyObject *orig_ndarray=NULL, *orig_arrayobj=NULL,
             *result_ndarray=NULL, *result_mask=NULL,
             *result_dict=NULL,
             *result_temp1=NULL, *result_temp2=NULL, *result_temp3=NULL,
             *mov_sum=NULL, *mov_sum_sq=NULL,
             *denom1=NULL, *denom2=NULL;

    PyArray_Descr *dtype=NULL;

    int rtype, span, is_variance, bias;

    static char *kwlist[] = {"array", "span", "is_variance", "bias", "dtype", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                "Oiii|O&:mov_stddev(array, span, is_variance, bias, dtype)",
                kwlist, &orig_arrayobj, &span, &is_variance, &bias,
                PyArray_DescrConverter2, &dtype)) return NULL;


    check_mov_args(orig_arrayobj, span, 2,
                   &orig_ndarray, &result_mask);

    rtype = _get_type_num_double(((PyArrayObject*)orig_ndarray)->descr, dtype);

    mov_sum = calc_mov_sum((PyArrayObject*)orig_ndarray, span, rtype);
    ERR_CHECK(mov_sum)

    result_temp1 = np_multiply(orig_ndarray, orig_ndarray);
    ERR_CHECK(result_temp1)

    mov_sum_sq = calc_mov_sum((PyArrayObject*)result_temp1, span, rtype);
    Py_DECREF(result_temp1);
    ERR_CHECK(mov_sum_sq)


    /*
      formulas from:
      http://en.wikipedia.org/wiki/Standard_deviation#Rapid_calculation_methods
     */
    if (bias == 0) {
        denom1 = PyFloat_FromDouble(1.0/(double)(span-1));
        denom2 = PyFloat_FromDouble(1.0/(double)(span*(span-1)));
    } else {
        denom1 = PyFloat_FromDouble(1.0/(double)span);
        denom2 = PyFloat_FromDouble(1.0/(double)(span*span));
    }

    result_temp1 = np_multiply(mov_sum_sq, denom1);
    ERR_CHECK(result_temp1)
    Py_DECREF(mov_sum_sq);
    Py_DECREF(denom1);

    result_temp3 = np_multiply(mov_sum, mov_sum);
    ERR_CHECK(result_temp3)
    Py_DECREF(mov_sum);

    result_temp2 = np_multiply(result_temp3, denom2);
    ERR_CHECK(result_temp2)
    Py_DECREF(result_temp3);
    Py_DECREF(denom2);

    result_temp3 = np_subtract(result_temp1, result_temp2);
    ERR_CHECK(result_temp3)
    Py_DECREF(result_temp1);
    Py_DECREF(result_temp2);

    if (is_variance) {
        result_ndarray = result_temp3;
    } else {
        result_temp1 = np_sqrt(result_temp3);
        ERR_CHECK(result_temp1)
        Py_DECREF(result_temp3);
        result_ndarray = result_temp1;
    }

    result_dict = PyDict_New();
    MEM_CHECK(result_dict)
    PyDict_SetItemString(result_dict, "array", result_ndarray);
    PyDict_SetItemString(result_dict, "mask", result_mask);

    Py_DECREF(result_ndarray);
    Py_DECREF(result_mask);
    return result_dict;
}

void import_c_tseries(PyObject *m) { import_array(); }
