#include "Python.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>

#include "numpy/libnumarray.h"

#define MAX_ARRAYS 1024

static PyObject *_Error;

typedef Float64 (*combiner)(int, int, int, Float64 temp[MAX_ARRAYS]);


static int
_mask_and_sort(int ninputs, int index, Float64 **inputs, UInt8 **masks,
	       Float64 temp[MAX_ARRAYS])
{
	int i, j, goodpix;
	if (masks) {
		for (i=j=0; i<ninputs; i++) {
			if (masks[i][index] == 0) 
				temp[j++] = inputs[i][index];
		}
	} else {
		for (i=j=0; i<ninputs; i++) {
			temp[j++] = inputs[i][index];
		}
	}
	goodpix = j;
	for(i=0; i<goodpix; i++) {
		for (j=i+1; j<goodpix; j++) {
			if (temp[j] < temp[i]) {
				Float64 temp2 = temp[j];
				temp[j] = temp[i];
				temp[i] = temp2;
			}
		}
	}
	return goodpix;
}

static Float64
_inner_median(int goodpix, int nlow, int nhigh, Float64 temp[MAX_ARRAYS])
{
	Float64 median;
	int midpoint, medianpix = goodpix-nhigh-nlow;
	if (medianpix <= 0) {
		median = 0;
	} else {
		midpoint = medianpix / 2;
		if (medianpix % 2) /* odd */ {
			median = temp[ midpoint + nlow ];
		} else {
			median = (temp[ midpoint + nlow ] + 
				  temp[ midpoint + nlow - 1 ]) / 2.0;
		}	
	}
	return median;
}

static Float64
_inner_average(int goodpix, int nlow, int nhigh, Float64 temp[MAX_ARRAYS])
{
	Float64 average;
	int i, averagepix = goodpix-nhigh-nlow;

	if (averagepix <= 0) {
		average = 0;
	} else {
		for(i=nlow, average=0; i<averagepix+nlow;  i++)
			average += temp[i];
		average /= averagepix;
	}
	return average;
}

static Float64
_inner_minimum(int goodpix, int nlow, int nhigh, Float64 temp[MAX_ARRAYS])
{
	int minimumpix = goodpix-nhigh-nlow;
	if (minimumpix <= 0) {
		return 0;
	} else {
	       return temp[nlow];
	}
}

static int
_combine(combiner f, int dim, int maxdim, int ninputs, int nlow, int nhigh,
	PyArrayObject *inputs[], PyArrayObject *masks[], PyArrayObject *output)
{
	int i, j;

	if (dim == maxdim-1) {
		Float64 sorted[MAX_ARRAYS];
		Float64 *tinputs[MAX_ARRAYS];
		UInt8    *tmasks[MAX_ARRAYS];
		Float64 *toutput;
		int cols = inputs[0]->dimensions[dim];

		/* Allocate and convert 1 temporary row at a time */
		for(i=0; i<ninputs; i++)
			tinputs[i] = (Float64 *) inputs[i]->data;
		if (masks) {
			for(i=0; i<ninputs; i++)
				tmasks[i] = (UInt8 *) masks[i]->data;
		}
		toutput = (Float64 *) output->data;
		
		for(j=0; j<cols; j++) {
			int goodpix = _mask_and_sort(
				ninputs, j, tinputs, masks ? tmasks : NULL, sorted);
			toutput[j] = f(goodpix, nlow, nhigh, sorted);
		}
	} else {
		for (i=0; i<inputs[0]->dimensions[dim]; i++) {
			for(j=0; j<ninputs; j++) {
				inputs[j]->data += inputs[j]->strides[dim]*i;
				if (masks) {
					masks[j]->data += masks[j]->strides[dim]*i;
				}
			}
			output->data += output->strides[dim]*i;
			_combine(f, dim+1, maxdim, ninputs, nlow, nhigh, 
				inputs, masks, output);
			for(j=0; j<ninputs; j++) {
				inputs[j]->data -= inputs[j]->strides[dim]*i;
				if (masks) {
					masks[j]->data -= masks[j]->strides[dim]*i;
				}
			}
			output->data -= output->strides[dim]*i;
		}
	}
	return 0;
}

typedef struct
{
	char *name;
	combiner fptr;
} fmapping;

static fmapping functions[] = {
	{"median", _inner_median},
	{"average", _inner_average},
	{"minimum", _inner_minimum},
};


static PyObject *
_Py_combine(PyObject *obj, PyObject *args, PyObject *kw)
{
	PyObject   *arrays, *output;
	int        nlow=0, nhigh=0, narrays;
	PyObject   *badmasks=Py_None;
	char       *keywds[] = { "arrays", "output", "nlow", "nhigh", 
				 "badmasks", "kind", NULL };
	char *kind;
	combiner f;
	PyArrayObject  *arr[MAX_ARRAYS], *bmk[MAX_ARRAYS], *toutput;
	int i;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "OO|iiOs:combine", keywds, 
			 &arrays, &output, &nlow, &nhigh, &badmasks, &kind))
		return NULL;

	narrays = PySequence_Length(arrays);
	if (narrays < 0) 
		return PyErr_Format(
			PyExc_TypeError, "combine: arrays is not a sequence");
	if (narrays > MAX_ARRAYS)
		return PyErr_Format(
			PyExc_TypeError, "combine: too many arrays.");

	for(i=0; i<narrays; i++) {
		PyObject *a = PySequence_GetItem(arrays, i);
		if (!a) return NULL;
		arr[i] = NA_InputArray(a, tFloat64, C_ARRAY);
		if (!arr[i]) return NULL;
		Py_DECREF(a);
		if (badmasks != Py_None) {
			a =  PySequence_GetItem(badmasks, i);
			if (!a) return NULL;	
			bmk[i] = NA_InputArray(a, tUInt8, C_ARRAY);
			if (!bmk[i]) return NULL;
			Py_DECREF(a);
		}
	}

	toutput = NA_OutputArray(output, tFloat64, C_ARRAY);
	if (!toutput) return NULL;
	
	for (i=0,f=0; i<ELEM(functions); i++)
		if  (!strcmp(kind, functions[i].name)) {
			f = functions[i].fptr;
			break;
		}
	if (!f)	return PyErr_Format(
		PyExc_ValueError, "Invalid comination function.");

	if (_combine( f, 0, arr[0]->nd, narrays, nlow, nhigh, 
		     arr, (badmasks != Py_None ? bmk : NULL), 
		     toutput) < 0)
		return NULL;

	for(i=0; i<narrays; i++) {
		Py_DECREF(arr[i]);
		if (badmasks != Py_None) {
			Py_DECREF(bmk[i]);
		}
	}
	Py_DECREF(toutput);

	Py_INCREF(Py_None);
	return Py_None;
}

static PyMethodDef _combineMethods[] = {
    {"combine", (PyCFunction) _Py_combine, METH_VARARGS | METH_KEYWORDS}, 
    {NULL, NULL} /* Sentinel */
};

PyMODINIT_FUNC init_combine(void)
{
	PyObject *m, *d;
	m = Py_InitModule("_combine", _combineMethods);
	d = PyModule_GetDict(m);
	_Error = PyErr_NewException("_combine.error", NULL, NULL);
	PyDict_SetItemString(d, "error", _Error);
	import_libnumarray();
}

/*
 * Local Variables:
 * mode: C
 * c-file-style: "python"
 * End:
 */
