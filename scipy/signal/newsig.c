#include <Python.h>

static int 
RawFilter2(const PyArrayObject *b, const PyArrayObject *a,
	   const PyArrayObject *x, const PyArrayObject *zi, 
	   const PyArrayObject *zf, PyArrayObject *y, int axis,
	   BasicFilterFunction *filter_func);

/*
 * XXX: Error checking not done yet
 */
static PyObject *
sigtools_linear_filter2(PyObject * dummy, PyObject * args)
{
	PyObject       *b = NULL, *a = NULL, *X = NULL, *Vi = NULL;
	PyArrayObject  *arY = NULL, *arb = NULL, *ara = NULL, *arX = NULL,
	               *arVi = NULL, *arVf = NULL;
	int             axis = -1, typenum, theaxis;
	char           *ara_ptr, input_flag = 0;
	intp na, nb;
	BasicFilterFunction *basic_filter;

	if (!PyArg_ParseTuple(args, "OOO|iO", &b, &a, &X, &axis, &Vi)) {
		return NULL;
	}

	typenum = PyArray_ObjectType(b, 0);
	typenum = PyArray_ObjectType(a, typenum);
	typenum = PyArray_ObjectType(X, typenum);
	if (Vi != NULL) {
		typenum = PyArray_ObjectType(Vi, typenum);
	}

	arY = NULL;
	arVf = NULL;
	ara = NULL;
	arb = NULL;
	arX = NULL;
	arVi = NULL;
	ara = (PyArrayObject *) PyArray_ContiguousFromObject(a, typenum, 1, 1);
	arb = (PyArrayObject *) PyArray_ContiguousFromObject(b, typenum, 1, 1);
	arX = (PyArrayObject *) PyArray_FromObject(X, typenum, 0, 0);
	if (ara == NULL || arb == NULL || arX == NULL) {
		goto fail;
	}

	if (axis < -arX->nd || axis > arX->nd - 1) {
		PyErr_SetString(PyExc_ValueError,
				"selected axis is out of range");
		goto fail;
	}
	if (axis < 0) {
		theaxis = arX->nd + axis;
	} else {
		theaxis = axis;
	}

	if (Vi != NULL) {
		arVi = (PyArrayObject *) PyArray_FromObject(Vi, typenum,
							    arX->nd, arX->nd);
		if (arVi == NULL)
			goto fail;
		input_flag = (PyArray_Size((PyObject *) arVi) > 0);
	}

	arY = (PyArrayObject *) PyArray_SimpleNew(arX->nd, 
						  arX->dimensions, typenum);
	if (arY == NULL) {
		goto fail;
	}

	if (input_flag) {
		arVf = (PyArrayObject *) PyArray_SimpleNew(arVi->nd, 
							   arVi->dimensions, 
						           typenum);
	}

	basic_filter = BasicFilterFunctions[(int) (arX->descr->type_num)];
	if (basic_filter == NULL) {
		PyErr_SetString(PyExc_ValueError,
				"linear_filter not available for this type");
		goto fail;
	}

	/* Skip over leading zeros in vector representing denominator (a) */
	// XXX: TODO
#if 0
	ara_ptr = ara->data;
	while (memcmp(ara_ptr, Va.zero, Va.elsize) == 0) {
		ara_ptr += Va.elsize;
		Va.data = ara_ptr;
		Va.numels--;
	}
#endif

	na = PyArray_SIZE(ara);
	nb = PyArray_SIZE(arb);
	if (input_flag) {
		if (arVi->dimensions[theaxis] != (na > nb ? na : nb) - 1) {
			PyErr_SetString(PyExc_ValueError,
					"The number of initial conditions must be max([len(a),len(b)]) - 1");
			goto fail;
		}
	}

	fprintf(stderr, "%s\n", __func__);
	RawFilter2(arb, ara, arX, arVi, arVf, arY, theaxis, basic_filter);

	Py_XDECREF(ara);
	Py_XDECREF(arb);
	Py_XDECREF(arX);
	Py_XDECREF(arVi);

	if (!input_flag) {
		return PyArray_Return(arY);
	} else {
		return Py_BuildValue("(NN)", arY, arVf);
	}


fail:
	Py_XDECREF(ara);
	Py_XDECREF(arb);
	Py_XDECREF(arX);
	Py_XDECREF(arVi);
	Py_XDECREF(arVf);
	Py_XDECREF(arY);
	return NULL;
}

static int
zfill(const PyArrayObject *x, intp nx, char* xzfilled, intp nxzfilled)
{
	char *xzero;
	intp i, nxl;

	nxl = PyArray_ITEMSIZE(x);

	xzero = PyArray_Zero(x);

	if (nx > 0) {
		memcpy(xzfilled, x->data, nx * nxl);
	}
	for(i = nx; i < nxzfilled; ++i) {
		memcpy(xzfilled + i * nxl, xzero, nxl);
	}

	PyDataMem_FREE(xzero);

	return 0;
}

/*
 * a and b assumed to be contiguous
 */
static int
RawFilter2(const PyArrayObject *b, const PyArrayObject *a,
	   const PyArrayObject *x, const PyArrayObject *zi, 
	   const PyArrayObject *zf, PyArrayObject *y, int axis,
	   BasicFilterFunction *filter_func)
{
	PyArrayIterObject *itx, *ity;
	intp nitx, i, nxl;
	intp na, nb, nal, nbl;
	intp nfilt;
	char *azfilled, *bzfilled, *zfzfilled;

	itx = (PyArrayIterObject *)PyArray_IterAllButAxis(
		(PyObject *)x, &axis);
	if (itx == NULL) {
		fprintf(stderr, "FAIL\n");
	}
	nitx = itx->size;

	ity = (PyArrayIterObject *)PyArray_IterAllButAxis(
		(PyObject *)y, &axis);
	if (ity == NULL) {
		fprintf(stderr, "FAIL\n");
	}

	na = PyArray_SIZE(a);
	nal = PyArray_ITEMSIZE(a);
	nb = PyArray_SIZE(b);
	nbl = PyArray_ITEMSIZE(b);

	nfilt = na > nb ? na : nb;

	azfilled = malloc(nal * nfilt);
	bzfilled = malloc(nbl * nfilt);

	nxl = PyArray_ITEMSIZE(x);
	zfzfilled = malloc(nxl * (nfilt-1) );

	zfill(a, na, azfilled, nfilt);
	zfill(b, nb, bzfilled, nfilt);

	if (zi != NULL) {
		fprintf(stderr, "%s: FAILS\n", __func__);
		return -1;
	} else {
		zfill(x, 0, zfzfilled, nfilt-1);
	}


#if 0
	fprintf(stderr, "%s: a and b are %f and %f\n", __func__,
((double*)azfilled)[0], ((double*)bzfilled)[0]);
	//fprintf(stderr, "%s: itx->size is %d\n", __func__, xsize);
#endif
	for(i = 0; i < nitx-1; ++i) {
#if 0
		fprintf(stderr, "item %d is %f, next is %d bytes away, "\
				"filter %d items\n", 
			i, ((double*)itx->dataptr)[0], itx->strides[axis],
			PyArray_DIM(x, axis));
#endif
		filter_func(bzfilled, azfilled, 
			    itx->dataptr, ity->dataptr, zfzfilled, 
			    nfilt, PyArray_DIM(x, axis), itx->strides[axis],
			    ity->strides[axis]);
		PyArray_ITER_NEXT(itx);
		PyArray_ITER_NEXT(ity);

		if (zi != NULL) {
			fprintf(stderr, "%s: FAIL\n", __func__);
			return -1;
		} else {
			/* XXX: inefficient because of the malloc in there */
			zfill(x, 0, zfzfilled, nfilt-1);
		}

	}
	//RawFilter(Vb, Va, x, y, vi, vf, basic_filter, theaxis);
	/* fprintf(stderr, "Now, Here.\n"); */

	Py_DECREF(ity);
	Py_DECREF(itx);

	return 0;
}
