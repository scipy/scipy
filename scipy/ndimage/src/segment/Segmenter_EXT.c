#include "ndImage_Segmenter_structs.h"
#include "Python.h"
#include "numpy/arrayobject.h"

static PyObject *Segmenter_EdgePreFilter(PyObject *self, PyObject *args)
{

    int lowThreshold;
    int highThreshold;
    int num;
    int nd;
    int type;
    int aperature;
    int half_taps;
    npy_intp *dims;
    unsigned short *fP1;
    double *fP2;
    double *pKernel;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;
    PyObject *kernel = NULL;

    //
    // pass in 2D LPF coefficients
    if(!PyArg_ParseTuple(args, "iiiOOO", &lowThreshold, &highThreshold, 
			 &half_taps, &kernel, &iArray, &eArray))
	    goto exit;

    fP1  = (unsigned short *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);
    fP2     = (double *)PyArray_DATA(eArray);
    pKernel = (double *)PyArray_DATA(kernel);
    aperature = PyArray_SIZE(kernel);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    
    if(!NI_EdgePreFilter(num, (int)dims[0], (int)dims[1], lowThreshold, 
		         highThreshold, aperature, half_taps, fP1, fP2, pKernel))
		      
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}


static PyObject *Segmenter_SobelImage(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    double *fP1;
    double *fP2;
    double pAve;
    int minValue;
    int maxValue;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    if(!PyArg_ParseTuple(args, "OO", &iArray, &eArray))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);
    fP2  = (double *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    
    if(!NI_SobelImage(num, (int)dims[0], (int)dims[1], fP1, fP2, &pAve, &minValue, &maxValue))
		      
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("dii", pAve, minValue, maxValue);

}


static PyObject *Segmenter_SobelEdges(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    double *fP1;
    unsigned short *fP2;
    double pAve;
    int minValue;
    int maxValue;
    int mode;
    double sobelLow;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    if(!PyArg_ParseTuple(args, "OOdiiid", &iArray, &eArray, &pAve, &minValue, &maxValue, &mode,
			                  &sobelLow))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);
    fP2  = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    
    if(!NI_SobelEdge(num, (int)dims[0], (int)dims[1], fP1, fP2, mode, pAve, minValue, maxValue, sobelLow))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}



static PyObject *Segmenter_GetBlobs(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int mask;
    npy_intp *dims;
    unsigned short *fP1;
    unsigned short *fP2;
    int groups;
    PyObject *iArray = NULL;
    PyObject *eArray = NULL;

    if(!PyArg_ParseTuple(args, "OOi", &iArray, &eArray, &mask))
	    goto exit;

    fP1  = (unsigned short *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);
    fP2  = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(iArray) || !PyArray_ISCONTIGUOUS(eArray))
	    goto exit;

    
    if(nd == 2){ 
        if(!NI_GetBlobs2D(num, (int)dims[0], (int)dims[1], fP1, fP2, &groups, mask))
	    goto exit;
    }
    else if(nd == 3){ 
        if(!NI_GetBlobs3D(num, (int)dims[0], (int)dims[1], (int)dims[2], fP1, fP2, 
			  &groups, mask))
	    goto exit;
    }

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("i", groups);

}

static PyObject *Segmenter_GetBlobRegions(PyObject *self, PyObject *args)
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

    if(!PyArg_ParseTuple(args, "OO", &iArray, &nArray))
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

    /* need to pass in 2D/3D flag and mask. NI_GetBlobRegions will call
     * 2D or 3D blob_extraction  */

    if(nd == 2){ 
        if(!NI_GetBlobRegions2D((int)dims[0], (int)dims[1], (int)objNumber[0], fP1, myData))
	        goto exit;
    }
    else if(nd == 3){ 
        if(!NI_GetBlobRegions3D((int)dims[0], (int)dims[1], (int)dims[2],
			       	(int)objNumber[0], fP1, myData))
	        goto exit;
    }

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}


static PyObject *Segmenter_ThinFilter(PyObject *self, PyObject *args)
{

    int number_masks;
    int roi_rows;
    int roi_cols;
    int cols;
    unsigned char *input;
    unsigned char *cinput;
    unsigned char *erosion;
    unsigned char *dialation;
    unsigned char *hmt;
    unsigned char *copy;
    unsigned short *j_mask;
    unsigned short *k_mask;
    PyObject  *jArray = NULL;
    PyObject  *kArray = NULL;
    PyObject  *iArray = NULL;
    PyObject  *cArray = NULL;
    PyObject  *eArray = NULL;
    PyObject  *dArray = NULL;
    PyObject  *hArray = NULL;
    PyObject  *pArray = NULL;

    if(!PyArg_ParseTuple(args, "OOiiiiOOOOOO", &jArray, &kArray, &number_masks, &roi_rows,
			 &roi_cols, &cols, &iArray, &cArray, &eArray, &dArray, &hArray, &pArray))
	    goto exit;


    j_mask = (unsigned short *)PyArray_DATA(jArray);
    k_mask = (unsigned short *)PyArray_DATA(kArray);

    input     = (unsigned char *)PyArray_DATA(iArray);
    cinput    = (unsigned char *)PyArray_DATA(cArray);
    erosion   = (unsigned char *)PyArray_DATA(eArray);
    dialation = (unsigned char *)PyArray_DATA(dArray);
    hmt       = (unsigned char *)PyArray_DATA(hArray);
    copy      = (unsigned char *)PyArray_DATA(pArray);

    if(!PyArray_ISCONTIGUOUS(iArray))
	    goto exit;

    if(!NI_ThinMorphoFilter(roi_rows, roi_cols, cols, number_masks, j_mask, k_mask,
		       input, cinput, erosion, dialation, hmt, copy))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}


static PyObject *Segmenter_CannyFilter(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int aperature;
    npy_intp *dims;
    double *fP1;
    double *h_DG_image;
    double *v_DG_image;
    double *pKernel;
    float aveXValue;
    float aveYValue;
    PyObject  *iArray = NULL;
    PyObject  *hArray = NULL;
    PyObject  *vArray = NULL;
    PyObject  *kernel = NULL;

    if(!PyArg_ParseTuple(args, "OOOOi", &iArray, &hArray, &vArray, &kernel, &aperature))
	    goto exit;

    fP1  = (double *)PyArray_DATA(iArray);
    nd   = PyArray_NDIM(iArray);
    dims = PyArray_DIMS(iArray);
    type = PyArray_TYPE(iArray);
    num  = PyArray_SIZE(iArray);

    h_DG_image = (double *)PyArray_DATA(hArray);
    v_DG_image = (double *)PyArray_DATA(vArray);
    pKernel    = (double *)PyArray_DATA(kernel);

    if(!PyArray_ISCONTIGUOUS(iArray))
	    goto exit;

    if(!NI_CannyFilter(num, (int)dims[0], (int)dims[1], fP1, h_DG_image,  
	               v_DG_image, pKernel, aperature, &aveXValue, &aveYValue))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("dd", aveXValue, aveYValue);

}




static PyObject *Segmenter_CannyNonMaxSupress(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int aperature;
    int mode;
    npy_intp *dims;
    double *h_DG_image;
    double *v_DG_image;
    double *magnitude;
    double aveXValue;
    double aveYValue;
    double aveMagnitude;
    double canny_low;
    double canny_high;
    double canny_l;
    double canny_h;
    PyObject  *mArray = NULL;
    PyObject  *hArray = NULL;
    PyObject  *vArray = NULL;

    if(!PyArg_ParseTuple(args, "OOOidddd", &hArray, &vArray, &mArray, &mode, &aveXValue,
			                   &aveYValue, &canny_l, &canny_h))
	    goto exit;

    magnitude = (double *)PyArray_DATA(mArray);
    nd   = PyArray_NDIM(mArray);
    dims = PyArray_DIMS(mArray);
    type = PyArray_TYPE(mArray);
    num  = PyArray_SIZE(mArray);

    h_DG_image = (double *)PyArray_DATA(hArray);
    v_DG_image = (double *)PyArray_DATA(vArray);

    if(!PyArray_ISCONTIGUOUS(mArray))
	    goto exit;

    if(!NI_CannyNonMaxSupress(num, (int)dims[0], (int)dims[1], magnitude, h_DG_image,  
	                      v_DG_image, mode, aveXValue, aveYValue, &aveMagnitude,
		              &canny_low, &canny_high, canny_l, canny_h))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("ddd", aveMagnitude, canny_low, canny_high);

}



static PyObject *Segmenter_CannyHysteresis(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    int aperature;
    int mode;
    npy_intp *dims;
    double *magnitude;
    unsigned short *hys_image;
    double canny_low;
    double canny_high;
    PyObject  *mArray = NULL;
    PyObject  *hArray = NULL;

    if(!PyArg_ParseTuple(args, "OOdd", &mArray, &hArray, &canny_low, &canny_high)) 
	    goto exit;

    magnitude = (double *)PyArray_DATA(mArray);
    nd   = PyArray_NDIM(mArray);
    dims = PyArray_DIMS(mArray);
    type = PyArray_TYPE(mArray);
    num  = PyArray_SIZE(mArray);

    hys_image = (unsigned short *)PyArray_DATA(hArray);

    if(!PyArray_ISCONTIGUOUS(mArray))
	    goto exit;

    if(!NI_CannyHysteresis(num, (int)dims[0], (int)dims[1], magnitude, hys_image,  
		           canny_low, canny_high))
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}


static PyObject *Segmenter_BinaryEdge(PyObject *self, PyObject *args)
{

    int num;
    int nd;
    int type;
    npy_intp *dims;
    unsigned short *mask_image;
    unsigned short *edge_image;
    PyObject  *mArray = NULL;
    PyObject  *eArray = NULL;

    if(!PyArg_ParseTuple(args, "OO", &mArray, &eArray)) 
	    goto exit;

    mask_image = (unsigned short *)PyArray_DATA(mArray);
    nd   = PyArray_NDIM(mArray);
    dims = PyArray_DIMS(mArray);
    type = PyArray_TYPE(mArray);
    num  = PyArray_SIZE(mArray);
    edge_image = (unsigned short *)PyArray_DATA(eArray);

    if(!PyArray_ISCONTIGUOUS(mArray))
	    goto exit;

    if(!NI_BinaryEdge(num, (int)dims[0], (int)dims[1], mask_image, edge_image))  
	    goto exit;

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *Segmenter_LawsTextureMetric(PyObject *self, PyObject *args)
{

    int i;
    int num;
    int nd;
    int type;
    int mode;
    npy_intp *dims;
    npy_intp *laws_dims;
    float  *lawsImage;
    double *src_image;
    unsigned short *mask;
    double *L7;
    double *E7;
    double *S7;
    double *W7;
    double *R7;
    double *O7;
    int number_kernels;
    int kernel_size;
    int filters;
    LawsFilter7 lawsFilter;
    PyObject *lArray = NULL;
    PyObject *mArray = NULL;
    PyObject *sArray = NULL;
    PyObject *LArray = NULL;
    PyObject *EArray = NULL;
    PyObject *SArray = NULL;
    PyObject *WArray = NULL;
    PyObject *RArray = NULL;
    PyObject *OArray = NULL;

    if(!PyArg_ParseTuple(args, "OOOiiiOOOOOO", &mArray, &sArray, &lArray, &number_kernels, 
			                       &kernel_size, &filters, &LArray, &EArray,
					       &SArray, &WArray, &RArray, &OArray))
	    goto exit;

    src_image = (double*)PyArray_DATA(sArray);
    nd   = PyArray_NDIM(sArray);
    dims = PyArray_DIMS(sArray);
    type = PyArray_TYPE(sArray);
    num  = PyArray_SIZE(sArray);

    laws_dims = PyArray_DIMS(lArray);
    mask      = (unsigned short *)PyArray_DATA(mArray);
    lawsImage = (float*)PyArray_DATA(lArray);
    L7        = (double *)PyArray_DATA(LArray);
    E7        = (double *)PyArray_DATA(EArray);
    S7        = (double *)PyArray_DATA(SArray);
    W7        = (double *)PyArray_DATA(WArray);
    R7        = (double *)PyArray_DATA(RArray);
    O7        = (double *)PyArray_DATA(OArray);

    lawsFilter.numberKernels      = number_kernels;
    lawsFilter.kernelLength       = kernel_size;
    lawsFilter.numberFilterLayers = filters;
    for(i = 0; i < kernel_size; ++i){
        lawsFilter.lawsKernel[0][i] = L7[i];
        lawsFilter.lawsKernel[1][i] = E7[i];
        lawsFilter.lawsKernel[2][i] = S7[i];
        lawsFilter.lawsKernel[3][i] = W7[i];
        lawsFilter.lawsKernel[4][i] = R7[i];
        lawsFilter.lawsKernel[5][i] = O7[i];
    }

    if(!PyArray_ISCONTIGUOUS(sArray)){
            printf("PyArray_ISCONTIGUOUS error\n");
	    goto exit;
    }

    if(!NI_LawsTexture(num, (int)dims[0], (int)dims[1], src_image, mask, lawsImage,   
		             lawsFilter)){
	    goto exit;
    }

exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *Segmenter_RoiCoOccurence(PyObject *self, PyObject *args)
{
    int num;
    int nd;
    int type;
    int distance;
    int orientation;
    npy_intp *dims;
    npy_intp *dims_cocm;
    unsigned short *mask_image;
    unsigned short *raw_image;
    int *coc_matrix;
    PyObject *mArray = NULL;
    PyObject *rArray = NULL;
    PyObject *cArray = NULL;

    if(!PyArg_ParseTuple(args, "OOOii", &mArray, &rArray, &cArray, &distance, &orientation)) 
	    goto exit;

    mask_image = (unsigned short *)PyArray_DATA(mArray);
    nd   = PyArray_NDIM(mArray);
    dims = PyArray_DIMS(mArray);
    type = PyArray_TYPE(mArray);
    num  = PyArray_SIZE(mArray);
    raw_image  = (unsigned short *)PyArray_DATA(rArray);
    coc_matrix = (int *)PyArray_DATA(cArray);
    dims_cocm  = PyArray_DIMS(cArray);

    if(!PyArray_ISCONTIGUOUS(mArray) || !PyArray_ISCONTIGUOUS(rArray)){
            printf("PyArray_ISCONTIGUOUS error\n");
	    goto exit;
    }

    if(!NI_RoiCoOccurence(num, (int)dims[0], (int)dims[1], mask_image, raw_image,
			  coc_matrix, distance, orientation))  
	    goto exit;



exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyObject *Segmenter_GrowRegion(PyObject *self, PyObject *args)
{


    int num;
    int nd;
    int type;
    int Label;
    int N_connectivity; 
    double low_threshold;
    double high_threshold;
    npy_intp *dims;
    npy_intp *objNumber;
    unsigned short *label;
    double *section;
    PyObject  *sArray = NULL;
    PyObject  *lArray = NULL;
    PyObject  *eArray = NULL;
    PyObject  *nArray = NULL;
    objStruct *expanded_ROI;
    objStruct *newgrow_ROI;

    if(!PyArg_ParseTuple(args, "OOOOddii", &sArray, &lArray, &eArray, &nArray, &low_threshold,
			 &high_threshold, &Label, &N_connectivity)){
            printf("PyArg_ParseTuple error\n");
	    goto exit;
    }

    if(!PyArray_ISCONTIGUOUS(sArray) || !PyArray_ISCONTIGUOUS(lArray)){
            printf("PyArray_ISCONTIGUOUS error\n");
	    goto exit;
    }

    section = (double *)PyArray_DATA(sArray);
    nd      = PyArray_NDIM(sArray);
    dims    = PyArray_DIMS(sArray);
    type    = PyArray_TYPE(sArray);
    num     = PyArray_SIZE(sArray);

    label        = (unsigned short *)PyArray_DATA(lArray);
    expanded_ROI = (objStruct*)PyArray_DATA(eArray);
    newgrow_ROI  = (objStruct*)PyArray_DATA(nArray);
	
    if(nd == 2){ 
        if(!NI_GrowRegion2D((int)dims[0], (int)dims[1], section, label, expanded_ROI,
			    newgrow_ROI, low_threshold, high_threshold, Label, N_connectivity))
	    goto exit;
    }
    else if(nd == 3){ 
        if(!NI_GrowRegion3D((int)dims[0], (int)dims[1], (int)dims[2], section, label,
			    expanded_ROI, newgrow_ROI, low_threshold, high_threshold,
			    Label, N_connectivity))
	    goto exit;
    }


exit:

    return PyErr_Occurred() ? NULL : (PyObject*)Py_BuildValue("");

}

static PyMethodDef SegmenterMethods[] =
{
    { "region_grow",          Segmenter_GrowRegion,         METH_VARARGS, NULL },
    { "roi_co_occurence",     Segmenter_RoiCoOccurence,     METH_VARARGS, NULL },
    { "binary_edge",          Segmenter_BinaryEdge,         METH_VARARGS, NULL },
    { "laws_texture_metric",  Segmenter_LawsTextureMetric,  METH_VARARGS, NULL },
    { "canny_hysteresis",     Segmenter_CannyHysteresis,    METH_VARARGS, NULL },
    { "canny_nonmax_supress", Segmenter_CannyNonMaxSupress, METH_VARARGS, NULL },
    { "canny_filter",         Segmenter_CannyFilter,        METH_VARARGS, NULL },
    { "sobel_edges",          Segmenter_SobelEdges,         METH_VARARGS, NULL },
    { "sobel_image",          Segmenter_SobelImage,         METH_VARARGS, NULL },
    { "edge_prefilter",       Segmenter_EdgePreFilter,      METH_VARARGS, NULL },
    { "get_blobs",            Segmenter_GetBlobs,           METH_VARARGS, NULL },
    { "get_blob_regions",     Segmenter_GetBlobRegions,     METH_VARARGS, NULL },
    { "thin_filter",          Segmenter_ThinFilter,         METH_VARARGS, NULL },
    {  NULL, NULL, 0, NULL},
};

void init_segment(void)
{
    Py_InitModule("_segment", SegmenterMethods);
    import_array();
}





