#include "ndImage_Segmenter_structs.h"
#include "Python.h"
#include "numpy/arrayobject.h"

static PyObject *NDI_Segmenter_CannyEdges(PyObject *self, PyObject *args)
{

    double sigma;
    double cannyLow;
    double cannyHigh;
    double BPHigh;
    int lowThreshold;
    int highThreshold;
    int apearture;
    int num;
    int nd;
    int type;
    int itype;
    int mode;
    int groups;
    npy_intp *dims;
    double *fP1;
    unsigned short *fP2;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    //
    // pass in 2D LPF coefficients
    if(!PyArg_Parse(args, "(dddiiidiO)", &sigma, &cannyLow, &cannyHigh, &mode, &lowThreshold, &highThreshold,
			                 &BPHigh, &apearture, &iArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    itype  = 4;
    eArray = (PyObject*)PyArray_SimpleNew(nd, dims, itype);
    fP2    = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    if(!NI_CannyEdges(num, (int)dims[0], (int)dims[1], sigma, cannyLow, cannyHigh, mode, lowThreshold,
		      highThreshold, BPHigh, apearture, fP1, fP2, &groups))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("Oi", eArray, groups);

}

static PyObject *NDI_Segmenter_SobelEdges(PyObject *self, PyObject *args)
{

    double sobelLow;
    double BPHigh;
    int lowThreshold;
    int highThreshold;
    int apearture;
    int num;
    int nd;
    int type;
    int itype;
    int groups;
    int mode;
    npy_intp *dims;
    double *fP1;
    unsigned short *fP2;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    //
    // pass in 2D LPF coefficients
    if(!PyArg_Parse(args, "(diiidiO)", &sobelLow, &mode, &lowThreshold, &highThreshold, &BPHigh, &apearture, &iArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    // this is int type and hard-wirred. pass this in from Python code
    itype  = 4; // unsigned short
    eArray = (PyObject*)PyArray_SimpleNew(nd, dims, itype);
    fP2    = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    
    if(!NI_SobelEdges(num, (int)dims[0], (int)dims[1], sobelLow, mode, lowThreshold, highThreshold, BPHigh, apearture,
		      fP1, fP2, &groups))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("Oi", eArray, groups-1);

}



static PyObject *NDI_Segmenter_ShenCastanEdges(PyObject *self, PyObject *args)
{
    int window;
    int lowThreshold;
    int highThreshold;
    double ShenCastanLow;
    double b;
    int num;
    int nd;
    int type;
    int itype;
    npy_intp *dims;
    double *fP1;
    unsigned short *fP2;
    int groups;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    if(!PyArg_Parse(args, "(ddiiiO)", &ShenCastanLow, &b, &window, &lowThreshold, &highThreshold, &iArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    // this is int type and hard-wirred. pass this in from Python code
    itype  = 4; // unsigned short
    eArray = (PyObject*)PyArray_SimpleNew(nd, dims, itype);
    fP2    = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    if(!NI_ShenCastanEdges(num, (int)dims[0], (int)dims[1], b, ShenCastanLow, window, lowThreshold, highThreshold, 
			   fP1, fP2, &groups))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("Oi", eArray, groups-1);

}

static PyObject *NDI_Segmenter_GetObjectStats(PyObject *self, PyObject *args)
{


    int num;
    int nd;
    int type;
    npy_intp *dims;
    npy_intp *objNumber;
    unsigned short *fP1;
    PyObject  *iArray = NULL;
    PyObject  *nArray = NULL;
    objStruct *myData;

    if(!PyArg_Parse(args, "(OO)", &iArray, &nArray))
	    goto exit;

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(nArray))
	    goto exit;

    	//
	//   PyArray_ContiguousFromObject or PyArray_ContiguousFromAny to be explored 
	//   for non-contiguous
	//

	
    // pointer to the edge-labeled image
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);
    fP1  = (unsigned short *)PyArray_DATA(iArray);

    // the object descriptor array that was allocated from numpy
    objNumber = PyArray_DIMS(nArray); // this is the number of labels in the edge image
    myData = (objStruct*)PyArray_DATA(nArray);

    if(!NI_GetObjectStats((int)dims[0], (int)dims[1], (int)objNumber[0], fP1, myData))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *NDI_Segmenter_MorphoThinFilt(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    npy_intp *objNumber;
    unsigned short *fP1;
    PyObject  *iArray = NULL;
    PyObject  *nArray = NULL;
    objStruct *ROIList;

    if(!PyArg_Parse(args, "(OO)", &iArray, &nArray))
	    goto exit;

    fP1  = (unsigned short *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    objNumber = PyArray_DIMS(nArray); // this is the number of labels in the edge image
    ROIList = (objStruct*)PyArray_DATA(nArray);

    if(!PyArray_ISCONTIGUOUS(iArray))
	    goto exit;

    if(!NI_ThinFilter(num, (int)dims[0], (int)dims[1], (int)objNumber[0], fP1, ROIList))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *NDI_Segmenter_BuildBoundary(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    npy_intp *objNumber;
    unsigned short *fP1;
    PyObject  *iArray = NULL;
    PyObject  *nArray = NULL;
    objStruct *ROIList;

    if(!PyArg_Parse(args, "(OO)", &iArray, &nArray))
	    goto exit;

    fP1  = (unsigned short *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);
    //
    // this is int type and hard-wirred. pass this in from Python code

    objNumber = PyArray_DIMS(nArray); // this is the number of labels in the edge image
    ROIList = (objStruct*)PyArray_DATA(nArray);

    if(!PyArray_ISCONTIGUOUS(iArray))
	    goto exit;

    //
    // pass in ROI list and labeled edges
    // return an augmented ROI list
    // replace the edgeImage with maskImage
    //
    if(!NI_BuildBoundary(num, (int)dims[0], (int)dims[1], (int)objNumber[0], fP1, ROIList))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}


static PyObject *NDI_Segmenter_VoxelMeasures(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    npy_intp *objNumber;
    double *fP1;
    unsigned short *fP2;
    PyObject  *iArray = NULL;
    PyObject  *nArray = NULL;
    PyObject  *eArray = NULL;
    objStruct *ROIList;

    if(!PyArg_Parse(args, "(OOO)", &iArray, &eArray, &nArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    // eArray and iArray are same dims
    fP2  = (unsigned short *)PyArray_DATA(eArray);

    objNumber = PyArray_DIMS(nArray); // this is the number of labels in the edge image
    ROIList = (objStruct*)PyArray_DATA(nArray);

    if(!PyArray_ISCONTIGUOUS(iArray))
	    goto exit;

    //
    // pass in ROI list and labeled edges
    // return an augmented ROI list
    // replace the edgeImage with maskImage
    //

    if(!NI_VoxelMeasures(num, (int)dims[0], (int)dims[1], (int)objNumber[0], fP1, fP2, ROIList))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *NDI_Segmenter_TextureMeasures(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    npy_intp *objNumber;
    double *fP1;
    unsigned short *fP2;
    PyObject  *iArray = NULL;
    PyObject  *nArray = NULL;
    PyObject  *eArray = NULL;
    objStruct *ROIList;

    if(!PyArg_Parse(args, "(OOO)", &iArray, &eArray, &nArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    // eArray and iArray are same dims
    fP2  = (unsigned short *)PyArray_DATA(eArray);

    objNumber = PyArray_DIMS(nArray); // this is the number of labels in the edge image
    ROIList = (objStruct*)PyArray_DATA(nArray);

    if(!PyArray_ISCONTIGUOUS(iArray))
	    goto exit;

    //
    // pass in ROI list and labeled edges
    // return an augmented ROI list
    // replace the edgeImage with maskImage
    //

    if(!NI_TextureMeasures(num, (int)dims[0], (int)dims[1], (int)objNumber[0], fP1, fP2, ROIList))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *NDI_Segmenter_RegionGrow(PyObject *self, PyObject *args)
{

    int lowThreshold;
    int highThreshold;
    int closeWindow;
    int openWindow;
    int num;
    int nd;
    int type;
    int itype;
    int groups;
    npy_intp *dims;
    double *fP1;
    unsigned short *fP2;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    //
    // pass in 2D LPF coefficients
    if(!PyArg_Parse(args, "(iiiiO)", &lowThreshold, &highThreshold, &closeWindow, &openWindow, &iArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    // this is int type and hard-wirred. pass this in from Python code
    itype  = 4; // unsigned short
    eArray = (PyObject*)PyArray_SimpleNew(nd, dims, itype);
    fP2    = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    
    if(!NI_RegionGrow(num, (int)dims[0], (int)dims[1], lowThreshold, highThreshold, closeWindow, openWindow,
		      fP1, fP2, &groups))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("Oi", eArray, groups-1);

}

static PyMethodDef NDI_SegmenterMethods[] =
{
    { "canny_edges",       NDI_Segmenter_CannyEdges,      METH_VARARGS },
    { "shen_castan_edges", NDI_Segmenter_ShenCastanEdges, METH_VARARGS },
    { "sobel_edges",       NDI_Segmenter_SobelEdges,      METH_VARARGS },
    { "get_object_stats",  NDI_Segmenter_GetObjectStats,  METH_VARARGS },
    { "morpho_thin_filt",  NDI_Segmenter_MorphoThinFilt,  METH_VARARGS },
    { "build_boundary",    NDI_Segmenter_BuildBoundary,   METH_VARARGS },
    { "voxel_measures",    NDI_Segmenter_VoxelMeasures,   METH_VARARGS },
    { "texture_measures",  NDI_Segmenter_TextureMeasures, METH_VARARGS },
    { "region_grow",       NDI_Segmenter_RegionGrow,      METH_VARARGS },
    {  NULL, NULL },
};

void init_segmenter(void)
{
    Py_InitModule("_segmenter", NDI_SegmenterMethods);
    import_array();
}

/*
static PyMethodDef NDI_SegmenterMethods[] =
{
    { "CannyEdges",       NDI_Segmenter_CannyEdges,      METH_VARARGS },
    { "ShenCastanEdges",  NDI_Segmenter_ShenCastanEdges, METH_VARARGS },
    { "SobelEdges",       NDI_Segmenter_SobelEdges,      METH_VARARGS },
    { "GetObjectStats",   NDI_Segmenter_GetObjectStats,  METH_VARARGS },
    { "MorphoThinFilt",   NDI_Segmenter_MorphoThinFilt,  METH_VARARGS },
    { "BuildBoundary",    NDI_Segmenter_BuildBoundary,   METH_VARARGS },
    { "VoxelMeasures",    NDI_Segmenter_VoxelMeasures,   METH_VARARGS },
    { "TextureMeasures",  NDI_Segmenter_TextureMeasures, METH_VARARGS },
    { "RegionGrow",       NDI_Segmenter_RegionGrow,      METH_VARARGS },
    {  NULL, NULL },
};

void initNDI_Segmenter()
{
    Py_InitModule("NDI_Segmenter", NDI_SegmenterMethods);
    import_array();
}
*/

