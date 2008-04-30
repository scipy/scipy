/* Python extension interface code */

#include "Python.h"
#include "numpy/arrayobject.h"

static PyObject *Register_Histogram(PyObject *self, PyObject *args)
{
     /*
       joint histogram memory is created in python to avoid memory leak problem 
     */

    int num;
    int numM;
    int nd;
    int type;
    int itype;
    int nd_histo;
    int nd_rotmatrix;
    int nd_S;
    npy_intp *dimsF;
    npy_intp *dimsG;
    npy_intp *dims_histo;
    npy_intp *dims_rotmatrix;
    npy_intp *dims_S;
    unsigned char *imageG;
    unsigned char *imageF;
    double        *pHisto;
    double        *M;
    int           *S;
    PyObject *imgArray1 = NULL;
    PyObject *imgArray2 = NULL;
    PyObject *rotArray  = NULL;
    PyObject *SArray    = NULL;
    PyObject *hArray    = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOOO", &imgArray1, &imgArray2, &rotArray, &SArray, &hArray))
	goto exit;

    /* check in the Python code that F and G are the same dims, type */
    imageF = (unsigned char *)PyArray_DATA(imgArray1);
    imageG = (unsigned char *)PyArray_DATA(imgArray2);
    nd     = PyArray_NDIM(imgArray1);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    dimsF  = PyArray_DIMS(imgArray1);
    dimsG  = PyArray_DIMS(imgArray2);
    type   = PyArray_TYPE(imgArray1);
    num    = PyArray_SIZE(imgArray1);

    M = (double *)PyArray_DATA(rotArray);
    nd_rotmatrix   = PyArray_NDIM(rotArray);
    dims_rotmatrix = PyArray_DIMS(rotArray);
    numM           = PyArray_SIZE(rotArray);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    pHisto = (double *)PyArray_DATA(hArray);
    nd_histo   = PyArray_NDIM(hArray);
    dims_histo = PyArray_DIMS(hArray);
    /* check to make sure this is 256x256  */

    if(!NI_Histogram2D((int)dimsF[0], (int)dimsF[1], (int)dimsF[2], 
                       (int)dimsG[0], (int)dimsG[1], (int)dimsG[2], 
		       S, M, imageG, imageF, pHisto))
	    goto exit;

exit:

    /* return the 2D histogram */
    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}


static PyObject *Register_HistogramLite(PyObject *self, PyObject *args)
{
     /*
       joint histogram memory is created in python to avoid memory leak problem 
     */

    int num;
    int numG;
    int nd;
    int type;
    int itype;
    int nd_histo;
    int nd_rotmatrix;
    int nd_S;
    npy_intp *dimsF;
    npy_intp *dimsG;
    npy_intp *dims_histo;
    npy_intp *dims_rotmatrix;
    npy_intp *dims_S;
    unsigned char *imageG;
    unsigned char *imageF;
    double        *pHisto;
    double        *M;
    int           *S;
    PyObject *imgArray1 = NULL;
    PyObject *imgArray2 = NULL;
    PyObject *rotArray  = NULL;
    PyObject *SArray    = NULL;
    PyObject *hArray    = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOOO", &imgArray1, &imgArray2, &rotArray, &SArray, &hArray))
	goto exit;

    /* check in the Python code that F and G are the same dims, type */
    imageF = (unsigned char *)PyArray_DATA(imgArray1);
    imageG = (unsigned char *)PyArray_DATA(imgArray2);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArray1);
    dimsF  = PyArray_DIMS(imgArray1);
    dimsG  = PyArray_DIMS(imgArray2);
    type   = PyArray_TYPE(imgArray1);
    num    = PyArray_SIZE(imgArray1);
    numG   = PyArray_SIZE(imgArray2);

    M = (double *)PyArray_DATA(rotArray);
    nd_rotmatrix   = PyArray_NDIM(rotArray);
    dims_rotmatrix = PyArray_DIMS(rotArray);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    pHisto = (double *)PyArray_DATA(hArray);
    nd_histo   = PyArray_NDIM(hArray);
    dims_histo = PyArray_DIMS(hArray);
    /* check to make sure this is 256x256  */

    if(!NI_Histogram2DLite((int)dimsF[0], (int)dimsF[1], (int)dimsF[2], 
                           (int)dimsG[0], (int)dimsG[1], (int)dimsG[2], 
		            S, M, imageG, imageF, pHisto))
	    goto exit;

exit:

    /* return the 2D histogram */
    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}

static PyObject *Register_VolumeResample(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int itype;
    int mode;
    int scale;
    npy_intp *dimsF;
    npy_intp *dimsG;
    unsigned char *imageG;
    unsigned char *imageF;
    double        *Z;
    PyObject *imgArray1 = NULL;
    PyObject *imgArray2 = NULL;
    PyObject *coordZoom = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOii", &imgArray1, &imgArray2, &coordZoom, &scale, &mode))
	goto exit;

    /* check in the Python code that F and G are the same dims, type */
    imageF = (unsigned char *)PyArray_DATA(imgArray1);
    imageG = (unsigned char *)PyArray_DATA(imgArray2);
    Z = (double *)PyArray_DATA(coordZoom);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArray1);
    dimsF  = PyArray_DIMS(imgArray1);
    dimsG  = PyArray_DIMS(imgArray2);
    type   = PyArray_TYPE(imgArray1);
    num    = PyArray_SIZE(imgArray1);

    if(!NI_VolumeResample((int)dimsF[0], (int)dimsF[1], (int)dimsF[2], 
                         (int)dimsG[0], (int)dimsG[1], (int)dimsG[2], 
		          scale, mode, imageG, imageF, Z))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}



static PyObject *Register_CubicResample(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int itype;
    int nd_rotmatrix;
    int nd_S;
    npy_intp *dimsF;
    npy_intp *dimsG;
    npy_intp *dims_rotmatrix;
    npy_intp *dims_S;
    unsigned char *imageG;
    unsigned char *imageF;
    double        *M;
    int           *S;
    PyObject *imgArray1 = NULL;
    PyObject *imgArray2 = NULL;
    PyObject *rotArray  = NULL;
    PyObject *SArray    = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOO", &imgArray1, &imgArray2, &rotArray, &SArray))
	goto exit;

    /* check in the Python code that F and G are the same dims, type */
    imageF = (unsigned char *)PyArray_DATA(imgArray1);
    imageG = (unsigned char *)PyArray_DATA(imgArray2);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArray1);
    dimsF  = PyArray_DIMS(imgArray1);
    dimsG  = PyArray_DIMS(imgArray2);
    type   = PyArray_TYPE(imgArray1);
    num    = PyArray_SIZE(imgArray1);

    M = (double *)PyArray_DATA(rotArray);
    nd_rotmatrix   = PyArray_NDIM(rotArray);
    dims_rotmatrix = PyArray_DIMS(rotArray);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    if(!NI_CubicResample((int)dimsF[0], (int)dimsF[1], (int)dimsF[2], 
                         (int)dimsG[0], (int)dimsG[1], (int)dimsG[2], 
		          S, M, imageG, imageF))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}


static PyObject *Register_LinearResample(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int itype;
    int nd_rotmatrix;
    int nd_S;
    npy_intp *dimsF;
    npy_intp *dimsG;
    npy_intp *dims_rotmatrix;
    npy_intp *dims_S;
    unsigned char *imageG;
    unsigned char *imageF;
    double        *M;
    int           *S;
    PyObject *imgArray1 = NULL;
    PyObject *imgArray2 = NULL;
    PyObject *rotArray  = NULL;
    PyObject *SArray    = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOO", &imgArray1, &imgArray2, &rotArray, &SArray))
	goto exit;

    /* check in the Python code that F and G are the same dims, type */
    imageF = (unsigned char *)PyArray_DATA(imgArray1);
    imageG = (unsigned char *)PyArray_DATA(imgArray2);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArray1);
    dimsF  = PyArray_DIMS(imgArray1);
    dimsG  = PyArray_DIMS(imgArray2);
    type   = PyArray_TYPE(imgArray1);
    num    = PyArray_SIZE(imgArray1);

    M = (double *)PyArray_DATA(rotArray);
    nd_rotmatrix   = PyArray_NDIM(rotArray);
    dims_rotmatrix = PyArray_DIMS(rotArray);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    if(!NI_LinearResample((int)dimsF[0], (int)dimsF[1], (int)dimsF[2], 
                          (int)dimsG[0], (int)dimsG[1], (int)dimsG[2], 
		           S, M, imageG, imageF))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}


static PyObject *Register_ImageThreshold(PyObject *self, PyObject *args)
{

    /* set threshold from the volume integrated histogram */
    int num;
    int nd;
    int type;
    int itype;
    int histogram_elements;
    int	tindex;
    npy_intp *dimsImage;
    npy_intp *dimsHistogram;
    unsigned short *image;
    double *H;
    double *IH;
    double threshold;
    PyObject *imgArray   = NULL;
    PyObject *histogram  = NULL;
    PyObject *ihistogram = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOd", &imgArray, &histogram, &ihistogram, &threshold))
	goto exit;

    image = (unsigned short *)PyArray_DATA(imgArray);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd        = PyArray_NDIM(imgArray);
    dimsImage = PyArray_DIMS(imgArray);
    type      = PyArray_TYPE(imgArray);
    num       = PyArray_SIZE(imgArray);

    H  = (double *)PyArray_DATA(histogram);
    IH = (double *)PyArray_DATA(ihistogram);
    histogram_elements = PyArray_SIZE(histogram);

    if(!NI_ImageThreshold((int)dimsImage[0], (int)dimsImage[1], (int)dimsImage[2], 
		               image, H, IH, histogram_elements, threshold, &tindex))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("i", tindex); 

}


static PyObject *Register_ResampleWithGradient(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int itype;
    int nd_rotmatrix;
    int nd_S;
    npy_intp *dimsScale;
    npy_intp *dimsOffset;
    npy_intp *dimsS;
    npy_intp *dimsD;
    npy_intp *dims_rotmatrix;
    npy_intp *dims_S;
    unsigned char *imageS;
    unsigned char *imageD;
    double        *M;
    int           *S;
    double        *scale;
    int           *offset;
    double        *gradientX;
    double        *gradientY;
    double        *gradientZ;
    PyObject *imgArrayS   = NULL;
    PyObject *imgArrayD   = NULL;
    PyObject *rotArray    = NULL;
    PyObject *SArray      = NULL;
    PyObject *scaleArray  = NULL;
    PyObject *offsetArray = NULL;
    PyObject *gradXArray  = NULL;
    PyObject *gradYArray  = NULL;
    PyObject *gradZArray  = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOOOOOOO", &imgArrayS, &imgArrayD, &rotArray, &SArray, &scaleArray,
			                    &offsetArray, &gradXArray, &gradYArray, &gradZArray))
	goto exit;

    /* check in the Python code that S and D are the same dims, type */
    imageS = (unsigned char *)PyArray_DATA(imgArrayS);
    imageD = (unsigned char *)PyArray_DATA(imgArrayD);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArrayS);
    dimsS  = PyArray_DIMS(imgArrayS);
    dimsD  = PyArray_DIMS(imgArrayD);
    type   = PyArray_TYPE(imgArrayS);
    num    = PyArray_SIZE(imgArrayS);

    M = (double *)PyArray_DATA(rotArray);
    nd_rotmatrix   = PyArray_NDIM(rotArray);
    dims_rotmatrix = PyArray_DIMS(rotArray);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    scale  = (double *)PyArray_DATA(scaleArray);
    offset = (int *)PyArray_DATA(offsetArray);
    dimsScale  = PyArray_DIMS(scaleArray);
    dimsOffset = PyArray_DIMS(offsetArray);

    gradientX = (double *)PyArray_DATA(gradXArray);
    gradientY = (double *)PyArray_DATA(gradYArray);
    gradientZ = (double *)PyArray_DATA(gradZArray);

    if(!NI_ResampleWithGradient((int)dimsS[0], (int)dimsS[1], (int)dimsS[2], 
                                (int)dimsD[0], (int)dimsD[1], (int)dimsD[2], 
		                 S, M, imageD, imageS, scale, offset, gradientX,
				 gradientY, gradientZ))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}


static PyObject *Register_Find_Mask(PyObject *self, PyObject *args)
{

    int i;
    int num;
    int length;
    double  *X;
    double  *Y;
    double  *Z;
    int     *xLims;
    int     *yLims;
    int     *zLims;
    int     *mask;
    PyObject *MArray   = NULL;
    PyObject *XArray   = NULL;
    PyObject *YArray   = NULL;
    PyObject *ZArray   = NULL;
    PyObject *XLimits  = NULL;
    PyObject *YLimits  = NULL;
    PyObject *ZLimits  = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOOOOO", &MArray, &XArray, &YArray, &ZArray, &XLimits, &YLimits, &ZLimits))
	goto exit;

    num   = PyArray_SIZE(XArray);
    X     = (double *)PyArray_DATA(XArray);
    Y     = (double *)PyArray_DATA(YArray);
    Z     = (double *)PyArray_DATA(ZArray);
    mask  = (int *)PyArray_DATA(MArray);
    xLims = (int *)PyArray_DATA(XLimits);
    yLims = (int *)PyArray_DATA(YLimits);
    zLims = (int *)PyArray_DATA(ZLimits);

    for(length = 0, i = 0; i < num; ++i){
	if( ((X[i] >= xLims[0]) && (X[i] <= xLims[1])) &&
	    ((Y[i] >= yLims[0]) && (Y[i] <= yLims[1])) &&
	    ((Z[i] >= zLims[0]) && (Z[i] <= zLims[1])) ){
	    mask[length++] = i;
	}
    } 


exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("i", length); 

}



static PyObject *Register_Resample_Gradient_Coords(PyObject *self, PyObject *args)
{

    int num;
    int size;
    int nd;
    int type;
    int itype;
    int nd_S;
    npy_intp *dimsScale;
    npy_intp *dimsOffset;
    npy_intp *dimsS;
    npy_intp *dimsD;
    npy_intp *dims_S;
    npy_intp *dims_Coords;
    unsigned char *imageS;
    unsigned char *imageD;
    double        *X;
    double        *Y;
    double        *Z;
    int           *S;
    double        *scale;
    int           *offset;
    double        *gradientX;
    double        *gradientY;
    double        *gradientZ;
    PyObject *imgArrayS    = NULL;
    PyObject *imgArrayD    = NULL;
    PyObject *SArray       = NULL;
    PyObject *scaleArray   = NULL;
    PyObject *offsetArray  = NULL;
    PyObject *gradXArray   = NULL;
    PyObject *gradYArray   = NULL;
    PyObject *gradZArray   = NULL;
    PyObject *coordXArray  = NULL;
    PyObject *coordYArray  = NULL;
    PyObject *coordZArray  = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOOOOOOOOO", &coordZArray, &coordYArray, &coordXArray,
                                &imgArrayS, &imgArrayD, &SArray, &scaleArray, &offsetArray,
			        &gradXArray, &gradYArray, &gradZArray))
	goto exit;

    /* check in the Python code that S and D are the same dims, type */
    imageS = (unsigned char *)PyArray_DATA(imgArrayS);
    imageD = (unsigned char *)PyArray_DATA(imgArrayD);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArrayS);
    dimsS  = PyArray_DIMS(imgArrayS);
    dimsD  = PyArray_DIMS(imgArrayD);
    type   = PyArray_TYPE(imgArrayS);
    num    = PyArray_SIZE(imgArrayS);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    scale  = (double *)PyArray_DATA(scaleArray);
    offset = (int *)PyArray_DATA(offsetArray);
    dimsScale  = PyArray_DIMS(scaleArray);
    dimsOffset = PyArray_DIMS(offsetArray);

    gradientX = (double *)PyArray_DATA(gradXArray);
    gradientY = (double *)PyArray_DATA(gradYArray);
    gradientZ = (double *)PyArray_DATA(gradZArray);

    X = (double *)PyArray_DATA(coordXArray);
    Y = (double *)PyArray_DATA(coordYArray);
    Z = (double *)PyArray_DATA(coordZArray);

    dims_Coords = PyArray_DIMS(coordXArray);
    size = PyArray_SIZE(coordXArray);

    if(!NI_Resample_Gradient_Coords(size, (int)dimsS[0], (int)dimsS[1], (int)dimsS[2], (int)dimsD[0], 
			           (int)dimsD[1], (int)dimsD[2], S, X, Y, Z, imageD, imageS, scale, 
				   offset, gradientX, gradientY, gradientZ))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}



static PyObject *Register_Resample_Coords(PyObject *self, PyObject *args)
{

    int num;
    int size;
    int nd;
    int type;
    int itype;
    int nd_S;
    npy_intp *dimsScale;
    npy_intp *dimsOffset;
    npy_intp *dimsS;
    npy_intp *dimsD;
    npy_intp *dims_S;
    npy_intp *dims_Coords;
    unsigned char *imageS;
    unsigned char *imageD;
    double        *X;
    double        *Y;
    double        *Z;
    int           *S;
    double        *scale;
    int           *offset;
    PyObject *imgArrayS    = NULL;
    PyObject *imgArrayD    = NULL;
    PyObject *SArray       = NULL;
    PyObject *scaleArray   = NULL;
    PyObject *offsetArray  = NULL;
    PyObject *coordXArray  = NULL;
    PyObject *coordYArray  = NULL;
    PyObject *coordZArray  = NULL;
	
    if(!PyArg_ParseTuple(args, "OOOOOOOO", &coordZArray, &coordYArray, &coordXArray,
                                &imgArrayS, &imgArrayD, &SArray, &scaleArray, &offsetArray))
	goto exit;

    /* check in the Python code that S and D are the same dims, type */
    imageS = (unsigned char *)PyArray_DATA(imgArrayS);
    imageD = (unsigned char *)PyArray_DATA(imgArrayD);
    /* reads dims as 0 = layers, 1 = rows, 2 = cols */
    nd     = PyArray_NDIM(imgArrayS);
    dimsS  = PyArray_DIMS(imgArrayS);
    dimsD  = PyArray_DIMS(imgArrayD);
    type   = PyArray_TYPE(imgArrayS);
    num    = PyArray_SIZE(imgArrayS);

    S = (int *)PyArray_DATA(SArray);
    nd_S   = PyArray_NDIM(SArray);
    dims_S = PyArray_DIMS(SArray);

    scale  = (double *)PyArray_DATA(scaleArray);
    offset = (int *)PyArray_DATA(offsetArray);
    dimsScale  = PyArray_DIMS(scaleArray);
    dimsOffset = PyArray_DIMS(offsetArray);

    X = (double *)PyArray_DATA(coordXArray);
    Y = (double *)PyArray_DATA(coordYArray);
    Z = (double *)PyArray_DATA(coordZArray);

    dims_Coords = PyArray_DIMS(coordXArray);
    size = PyArray_SIZE(coordXArray);

    if(!NI_Resample_Coords(size, (int)dimsS[0], (int)dimsS[1], (int)dimsS[2], (int)dimsD[0], 
		           (int)dimsD[1], (int)dimsD[2], S, X, Y, Z, imageD, imageS, scale, offset)) 
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue(""); 

}



static PyMethodDef RegisterMethods[] =
{
    { "register_find_mask",                     Register_Find_Mask,                METH_VARARGS, NULL },
    { "register_resample_coords",               Register_Resample_Coords,          METH_VARARGS, NULL },
    { "register_resample_gradient_coords",      Register_Resample_Gradient_Coords, METH_VARARGS, NULL },
    { "register_resample_w_gradient",           Register_ResampleWithGradient,     METH_VARARGS, NULL },
    { "register_histogram",                     Register_Histogram,                METH_VARARGS, NULL },
    { "register_histogram_lite",                Register_HistogramLite,            METH_VARARGS, NULL },
    { "register_linear_resample",               Register_LinearResample,           METH_VARARGS, NULL },
    { "register_cubic_resample",                Register_CubicResample,            METH_VARARGS, NULL },
    { "register_volume_resample",               Register_VolumeResample,           METH_VARARGS, NULL },
    { "register_image_threshold",               Register_ImageThreshold,           METH_VARARGS, NULL },
    {  NULL, NULL, 0, NULL},
};

PyMODINIT_FUNC init_register(void)
{
    Py_InitModule("_register", RegisterMethods);
    import_array();
}


