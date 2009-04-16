#include <Python.h>

static int
RawFilter(const PyArrayObject *b, const PyArrayObject *a,
	   const PyArrayObject *x, const PyArrayObject *zi,
	   const PyArrayObject *zf, PyArrayObject *y, int axis,
	   BasicFilterFunction *filter_func);

static char doc_linear_filter[] = "(y,Vf) = _linear_filter(b,a,X,Dim=-1,Vi=None)  implemented using Direct Form II transposed flow diagram. If Vi is not given, Vf is not returned.";
 
/*
 * XXX: Error checking not done yet
 */
static PyObject *
sigtools_linear_filter(PyObject * dummy, PyObject * args)
{
	PyObject *b, *a, *X, *Vi;
	PyArrayObject *arY, *arb, *ara, *arX, *arVi, *arVf;
	int axis, typenum, theaxis, st;
	char *ara_ptr, input_flag = 0, *azero;
	intp na, nb, nal;
	BasicFilterFunction *basic_filter;

        axis = -1;
        Vi = NULL;
	if (!PyArg_ParseTuple(args, "OOO|iO", &b, &a, &X, &axis, &Vi)) {
		return NULL;
	}

	typenum = PyArray_ObjectType(b, 0);
	typenum = PyArray_ObjectType(a, typenum);
	typenum = PyArray_ObjectType(X, typenum);
	if (Vi != NULL) {
		typenum = PyArray_ObjectType(Vi, typenum);
	}

	arY = arVf = arVi = NULL;
	ara = (PyArrayObject *) PyArray_ContiguousFromObject(a, typenum, 1, 1);
	arb = (PyArrayObject *) PyArray_ContiguousFromObject(b, typenum, 1, 1);
	arX = (PyArrayObject *) PyArray_FromObject(X, typenum, 0, 0);
        /* XXX: fix failure handling here */
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
                Py_ssize_t nvi;
		arVi = (PyArrayObject *) PyArray_FromObject(Vi, typenum,
							    arX->nd, arX->nd);
		if (arVi == NULL)
			goto fail;

                nvi = PyArray_Size((PyObject *) arVi);
                if (nvi > 0) {
                        input_flag = 1;
                } else {
                        input_flag = 0;
                        Py_DECREF(arVi);
                        arVi = NULL;
                }
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
	/* XXX: handle this correctly */
	azero = PyArray_Zero(ara);
	ara_ptr = ara->data;
        nal = PyArray_ITEMSIZE(ara);
	if (memcmp(ara_ptr, azero, nal) == 0) {
		PyErr_SetString(PyExc_ValueError,
				"BUG: filter coefficient a[0] == 0 not supported yet");
                goto fail;
	}
	PyDataMem_FREE(azero);

	na = PyArray_SIZE(ara);
	nb = PyArray_SIZE(arb);
	if (input_flag) {
		if (arVi->dimensions[theaxis] != (na > nb ? na : nb) - 1) {
			PyErr_SetString(PyExc_ValueError,
					"The number of initial conditions must be max([len(a),len(b)]) - 1");
			goto fail;
		}
	}

	st = RawFilter(arb, ara, arX, arVi, arVf, arY, theaxis, basic_filter);
	if (st) {
		goto fail;
	}

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

	/* PyArray_Zero does not take const pointer, hence the cast */
	xzero = PyArray_Zero((PyArrayObject*)x);

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
 *
 * XXX: this code is very conservative, and could be considerably sped up for
 * the usual cases (like contiguity).
 *
 * XXX: the code should be refactored (at least with/without initial
 * condition), some code is wasteful here
 */
static int
RawFilter(const PyArrayObject *b, const PyArrayObject *a,
	   const PyArrayObject *x, const PyArrayObject *zi,
	   const PyArrayObject *zf, PyArrayObject *y, int axis,
	   BasicFilterFunction *filter_func)
{
	PyArrayIterObject *itx, *ity, *itzi, *itzf;
	intp nitx, i, nxl, nzfl, j;
	intp na, nb, nal, nbl;
	intp nfilt;
	char *azfilled, *bzfilled, *zfzfilled, *yoyo;
	PyArray_CopySwapFunc *copyswap = x->descr->f->copyswap;

	itx = (PyArrayIterObject *)PyArray_IterAllButAxis(
		(PyObject *)x, &axis);
	if (itx == NULL) {
		PyErr_SetString(PyExc_MemoryError,
				"Could not create itx");
		goto fail;
	}
	nitx = itx->size;

	ity = (PyArrayIterObject *)PyArray_IterAllButAxis(
		(PyObject *)y, &axis);
	if (ity == NULL) {
		PyErr_SetString(PyExc_MemoryError,
				"Could not create ity");
		goto clean_itx;
	}

        if (zi != NULL) {
                itzi = (PyArrayIterObject *)PyArray_IterAllButAxis(
                        (PyObject *)zi, &axis);
                if (itzi == NULL) {
			PyErr_SetString(PyExc_MemoryError,
					"Could not create itzi");
			goto clean_ity;
                }

                itzf = (PyArrayIterObject *)PyArray_IterAllButAxis(
                        (PyObject *)zf, &axis);
                if (itzf == NULL) {
			PyErr_SetString(PyExc_MemoryError,
					"Could not create itzf");
			goto clean_itzi;
                }
        }

	na = PyArray_SIZE(a);
	nal = PyArray_ITEMSIZE(a);
	nb = PyArray_SIZE(b);
	nbl = PyArray_ITEMSIZE(b);

	nfilt = na > nb ? na : nb;

	azfilled = malloc(nal * nfilt);
	if (azfilled == NULL) {
		PyErr_SetString(PyExc_MemoryError,
				"Could not create azfilled");
		goto clean_itzf;
	}
	bzfilled = malloc(nbl * nfilt);
	if (bzfilled == NULL) {
		PyErr_SetString(PyExc_MemoryError,
				"Could not create bzfilled");
		goto clean_azfilled;
	}

	nxl = PyArray_ITEMSIZE(x);
	zfzfilled = malloc(nxl * (nfilt-1) );
	if (zfzfilled == NULL) {
		PyErr_SetString(PyExc_MemoryError,
				"Could not create zfzfilled");
		goto clean_bzfilled;
	}
	/* Initialize zfzilled to 0, so that we can use Py_XINCREF/Py_XDECREF
	 * on it for object arrays (necessary for copyswap to work correctly).
	 * Stricly speaking, it is not needed for fundamental types (as values
	 * are copied instead of pointers, without refcounts), but oh well...
	 */
	memset(zfzfilled, 0, nxl * (nfilt-1));

	zfill(a, na, azfilled, nfilt);
	zfill(b, nb, bzfilled, nfilt);

        /* XXX: Check that zf and zi have same type ? */
        if (zf != NULL) {
                nzfl = PyArray_ITEMSIZE(zf);
        } else {
                nzfl = 0;
        }

        /* Iterate over the input array */
        for(i = 0; i < nitx; ++i) {
                if (zi != NULL) {
                        yoyo = itzi->dataptr;
                        /* Copy initial conditions zi in zfzfilled buffer */
                        for(j = 0; j < nfilt - 1; ++j) {
                                copyswap(zfzfilled + j * nzfl, yoyo, 0, NULL);
                                yoyo += itzi->strides[axis];
                        }
                        PyArray_ITER_NEXT(itzi);
                } else {
                        zfill(x, 0, zfzfilled, nfilt-1);
                }

                filter_func(bzfilled, azfilled,
                            itx->dataptr, ity->dataptr, zfzfilled,
                            nfilt, PyArray_DIM(x, axis), itx->strides[axis],
                            ity->strides[axis]);
                PyArray_ITER_NEXT(itx);
                PyArray_ITER_NEXT(ity);

                /* Copy tmp buffer fo final values back into zf output array */
                if (zi != NULL) {
                        yoyo = itzf->dataptr;
                        for(j = 0; j < nfilt - 1; ++j) {
                                copyswap(yoyo, zfzfilled + j * nzfl, 0, NULL);
                                yoyo += itzf->strides[axis];
                        }
                        PyArray_ITER_NEXT(itzf);
                }
	}

	/* Free up allocated memory */
	free(zfzfilled);
	free(bzfilled);
	free(azfilled);

	if (zi != NULL) {
                Py_DECREF(itzf);
                Py_DECREF(itzi);
        }
	Py_DECREF(ity);
	Py_DECREF(itx);

	return 0;

clean_bzfilled:
	free(bzfilled);
clean_azfilled:
	free(azfilled);
clean_itzf:
	if (zf != NULL) {
		Py_DECREF(itzf);
	}
clean_itzi:
	if (zi != NULL) {
		Py_DECREF(itzi);
	}
clean_ity:
	Py_DECREF(ity);
clean_itx:
	Py_DECREF(itx);
fail:
	return -1;
}
