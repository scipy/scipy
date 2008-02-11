#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ndImage_Segmenter_structs.h"

// these are for this standalone and come out with the full build
//
#define MAX(a, b) ((a) > (b) ? (a) : (b)) 
#define FALSE 0
#define TRUE  1

int NI_GetObjectStats(int rows, int cols, int numberObjects, unsigned short *labeledEdges,
                      objStruct objectMetrics[]){

	int i, j, k, m;
	int offset;
	int count;
	int LowX;
	int LowY;
	int HighX;
	int HighY;
	int status;
	float centerX;
	float centerY;

	for(k = 1; k < numberObjects; ++k){
	    offset     = cols;
	    LowX       = 32767;
	    LowY       = 32767;
	    HighX      = 0;
	    HighY      = 0;
	    count      = 0;
	    centerX    = (float)0.0;
	    centerY    = (float)0.0;
	    for(i = 1; i < (rows-1); ++i){
		for(j = 1; j < (cols-1); ++j){
		    m = labeledEdges[offset+j];
		    if(k == m){
			if(i < LowY)   LowY = i;
			if(j < LowX)   LowX = j;
			if(i > HighY) HighY = i;
			if(j > HighX) HighX = j;
	    		centerX += (float)j;
	    		centerY += (float)i;
	    		++count;
		    }
		}
		offset += cols;
	    }
	    /* the bounding box for the 2D blob */
	    objectMetrics[k-1].L     = LowX;
	    objectMetrics[k-1].R     = HighX;
	    objectMetrics[k-1].B     = LowY;
	    objectMetrics[k-1].T     = HighY;
	    objectMetrics[k-1].Area  = count;
	    objectMetrics[k-1].cX    = centerX/(float)count;
	    objectMetrics[k-1].cY    = centerY/(float)count;
	    objectMetrics[k-1].Label = k;
	}

	status = numberObjects;
	return status;

}


void buildKernel(double BPHigh, int HalfFilterTaps, int apearture, float *kernel){

	int i, j;
	float r, t1, t2, t3, t4;
	float LC, HC, tLOW, tHIGH;
	float pi = (float)3.14159, rad = (float)0.01745;

	LC = (float)0.0;
	HC = BPHigh * rad; 
	t2 = (float)2.0*pi; 
	t1 = (float)2.0*HalfFilterTaps + (float)1.0;
	/*
	// build the Filter Kernel 
	// the kernel starts at 1 only because it is linked to the internal filter2D routine
	// the code is not a Fortran code
	*/
	j = 1;
	for(i = -HalfFilterTaps; i <= HalfFilterTaps; ++i){
	    r = (float)i;
	    if(r == (float)0.0){
		tLOW  = LC;
	        tHIGH = HC;
	    }
	    else{
		tLOW  = (float)(sin(r*LC))/r;
	        tHIGH = (float)(sin(r*HC))/r;
	    }
	    t3 = (float)0.54 + (float)0.46*((float)cos(r*t2/t1));
	    t4 = t3*(tHIGH-tLOW);
	    kernel[j++] = t4;
	}

	/* normalize the kernel so unity gain (as is LP filter this is easy) */
	t1 = (float)0.0;
	for(j = 1; j <= apearture; ++j){  
	    t1 += kernel[j];
	}
	for(j = 1; j <= apearture; ++j){  
	    kernel[j] /= t1;
	}

	t1 = (float)0.0;
	for(j = 1; j <= apearture; ++j){  
	    t1 += kernel[j];
	}
	return;
}

void filter2D(int HalfFilterTaps, int rows, int cols, int lowThreshold, int highThreshold,
              float *kernel, double *Image){

	int i, j, k, n, num1;
    	int offset;
	float sum, value;
	float buffer[1024];

	num1 = HalfFilterTaps + 1;
	offset = 0;
	for(i = 0; i < rows; ++i){
	    /* copy image row to local buffer  */
	    for(j = 0; j < cols; ++j){
		buffer[num1+j] = Image[offset+j];
	    }
	    /* constant pad the ends of the buffer */
	    for(j = 0; j < num1; ++j){
		buffer[j] = buffer[num1];
	    }
	    for(j = cols+num1; j < cols+2*num1; ++j){
		buffer[j] = buffer[cols-1+num1];
	    }

	    /* Perform Symmetric Convolution in the X dimension. */
	    for(n = 0, j = num1; j < (cols+num1); ++j, ++n){
	        sum = buffer[j] * kernel[num1];
	        for(k = 1; k < num1; ++k){
	            sum += kernel[num1-k] * (buffer[j+k] + buffer[j-k]);
	        }
	        Image[offset+n] = sum;
	    }
	    offset += cols;
	}

	offset = 0;
	for(i = 0; i < cols; ++i){
	    /* copy image column to local buffer */
	    offset = 0;
	    for(j = 0; j < rows; ++j){
            buffer[num1+j] = Image[offset+i];
	        offset += cols;
	    }
	    /* constant pad the ends of the buffer */
	    for(j = 0; j < num1; ++j){
		buffer[j] = buffer[num1];
	    }
	    for(j = rows+num1; j < rows+2*num1; ++j){
	        buffer[j] = buffer[rows-1+num1];
	    }

	    /* Perform Symmetric Convolution in the Y dimension. */
	    offset = 0;
	    for(j = num1; j < (rows+num1); ++j){
	        sum = buffer[j] * kernel[num1];
	        for(k = 1; k < num1; ++k){
	            sum += kernel[num1-k] * (buffer[j+k] + buffer[j-k]);
	        }
	        Image[offset+i] = sum;
	        offset += cols;
	    }
	}

	/* threshold the image */
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		value = Image[offset+j];
		if(value < (float)lowThreshold)  value = (float)0.0;
		if(value > (float)highThreshold) value = (float)0.0;
		Image[offset+j] = value;
	    }
	    offset += cols;
	}

	return;

}

void doPreProcess(int samples, int rows, int cols, double *rawImage, double BPHigh, 
                  int apearture, int lowThreshold, int highThreshold){

	/*
	// 2D low pass filter using bisinc and threshold 
	// this specific example is on cardiac CT and focuses on segmenting the
	// aorta and blood-filled chambers. for MRI the threshold will be different
	*/

	float *kernel;
	int HalfFilterTaps = (apearture-1)/2;
	kernel = calloc(apearture+16, sizeof(float));

	buildKernel(BPHigh, HalfFilterTaps, apearture, kernel);
	filter2D(HalfFilterTaps, rows, cols, lowThreshold, highThreshold, kernel, rawImage);

	free(kernel);

	return;

}


int ConnectedEdgePoints(int rows, int cols, unsigned short *connectedEdges){

	int            i, j, k, l, m;
	int            offset;
	int            Label;
	int            Classes[4096];
	bool           NewLabel;
	bool           Change;
	unsigned short T[12];

	/*
	// connected components labeling. pixels touch within 3x3 mask for edge connectedness. 
	*/
	Label  = 1;
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(connectedEdges[offset+j] == 1){
		    connectedEdges[offset+j] = Label++; 
		}
	    }
	    offset += cols;
	}

	while(1){
	    Change = FALSE;
	    /*
	    // TOP-DOWN Pass for labeling
	    */
	    offset = cols;
	    for(i = 1; i < rows-1; ++i){
		for(j = 1; j < cols-1; ++j){
		    if(connectedEdges[offset+j] != 0){
			T[0] = connectedEdges[offset+j];
			T[1] = connectedEdges[offset+j+1];
			T[2] = connectedEdges[offset-cols+j+1];
			T[3] = connectedEdges[offset-cols+j];
			T[4] = connectedEdges[offset-cols+j-1];
			T[5] = connectedEdges[offset+j-1];
			T[6] = connectedEdges[offset+cols+j-1];
			T[7] = connectedEdges[offset+cols+j];
			T[8] = connectedEdges[offset+cols+j+1];
			m = T[0];
			for(l = 1; l < 9; ++l){
			    if(T[l] != 0){
				if(T[l] < m) m = T[l];
			    }
			}
			if(m != connectedEdges[offset+j]){
			    Change = TRUE;
			    connectedEdges[offset+j] = m;
			}
		    }
		}
		offset += cols;
	    }
	    /*
	    // BOTTOM-UP Pass for labeling
	    */
	    offset = (rows-1)*cols;
	    for(i = (rows-1); i > 1; --i){
		for(j = (cols-1); j > 1; --j){
		    if(connectedEdges[offset+j] != 0){
			T[0] = connectedEdges[offset+j];
			T[1] = connectedEdges[offset+j+1];
			T[2] = connectedEdges[offset-cols+j+1];
			T[3] = connectedEdges[offset-cols+j];
			T[4] = connectedEdges[offset-cols+j-1];
			T[5] = connectedEdges[offset+j-1];
			T[6] = connectedEdges[offset+cols+j-1];
			T[7] = connectedEdges[offset+cols+j];
			T[8] = connectedEdges[offset+cols+j+1];
			m = T[0];
			for(l = 1; l < 9; ++l){
			    if(T[l] != 0){
				if(T[l] < m) m = T[l];
			    }
			}
			if(m != connectedEdges[offset+j]){
			    Change = TRUE;
			    connectedEdges[offset+j] = m;
			}
		    }
		}
		offset -= cols;
	    }
	    if(!Change) break;
	}   /* end while loop */

	Classes[0] = 0;
	Label      = 1;
	offset     = cols;
	for(i = 1; i < (rows-1); ++i){
	    for(j = 1; j < (cols-1); ++j){
		m = connectedEdges[offset+j];
		if(m > 0){
		    NewLabel = TRUE;
		    for(k = 1; k < Label; ++k){
			if(Classes[k] == m) NewLabel = FALSE;
		    }
		    if(NewLabel){
			Classes[Label++] = m;
			if(Label > 4000){
			    return 0; /* too many labeled regions. this is a pathology */
			}
		    }
		}
	    }
	    offset += cols;
	}

	/*
	// re-label the connected blobs in continuous label order
	*/
	offset = cols;
	for(i = 1; i < (rows-1); ++i){
	    for(j = 1; j < (cols-1); ++j){
		m = connectedEdges[offset+j];
		if(m > 0){
		    for(k = 1; k < Label; ++k){
			if(Classes[k] == m){
			    connectedEdges[offset+j] = (unsigned short)k;
			    break;
			}
		    }
		}
	    }
	    offset += cols;
	}

	return Label;
}

float magnitude(float X, float Y){

	return (float)sqrt(X*X + Y*Y);
}

int traceEdge(int i, int j, int rows, int cols, double cannyLow, float *magImage,
              float *HYSImage){

	int n, m;
	int ptr;
	int flag;

	ptr = i * cols;
	if(HYSImage[ptr+j] == (float)0.0){
	    /*
	    // this point is above high threshold
	    */
	    HYSImage[ptr+j] = (float)1.0;
	    flag = 0;
	    for(n = -1; n <= 1; ++n){
		for(m = -1; m <= 1; ++m){
		    if(n == 0 && m == 0) continue;
		    if(((i+n) > 0) && ((j+m) > 0) && ((i+n) < rows) && ((j+m) < cols)){
			ptr = (i+n) * cols;
			if(magImage[ptr+j+m] > cannyLow){
	    		    /*
	    		    // this point is above low threshold
	    		    */
			    if(traceEdge(i+n, j+m, rows, cols, cannyLow, magImage, HYSImage)){
				flag = 1;
				break;
			    }
			}
		    }
		}
		if(flag) break;
	    }
	    return(1);
	}

	return(0);

}


void edgeThreshold(int rows, int cols, double cannyLow, float *magImage, 
                   float *HYSImage){

	int i, j;
	int ptr;

	for(i = 0; i < rows; ++i){
	    ptr = i * cols;
	    for(j = 0; j < cols; ++j){
		if(magImage[ptr+j] > cannyLow){
		    HYSImage[ptr+j] = (float)1.0;
		}
	    }
	}

	return;

}

void edgeHysteresis(int rows, int cols, double cannyLow, double cannyHigh,
                    float *magImage, float *HYSImage){

	int i, j;
	int ptr;

	for(i = 0; i < rows; ++i){
	    ptr = i * cols;
	    for(j = 0; j < cols; ++j){
		if(magImage[ptr+j] > cannyHigh){
		    traceEdge(i, j, rows, cols, cannyLow, magImage, HYSImage);
		}
	    }
	}

	return;

}

void nonMaxSupress(int rows, int cols, float aveXValue, float aveYValue,
                   double *cannyLow, double *cannyHigh, int mode, 
                   float *hDGImage, float *vDGImage, float *magImage){

	int i, j;
	int ptr, ptr_m1, ptr_p1;
	float xSlope, ySlope, G1, G2, G3, G4, G, xC, yC;
	float scale;
	float maxValue = (float)0.0;
	float minValue = (float)-1.0;
	int histogram[256];
	int value;
	int mValue;
	int mIndex;
	int count;
	double step;
	double tAve;

	for(i = 1; i < rows-1; ++i){
	    ptr = i * cols;
	    ptr_m1 = ptr - cols;
	    ptr_p1 = ptr + cols;
	    for(j = 1; j < cols; ++j){
		magImage[ptr+j] = (float)0.0;
		xC = hDGImage[ptr+j];
		yC = vDGImage[ptr+j];
		if((fabs(xC) < aveXValue) && (fabs(yC) < aveYValue)) continue;
		G = magnitude(xC, yC);
		if(fabs(yC) > fabs(xC)){
		    /* vertical gradient */
		    xSlope = (float)(fabs(xC) / fabs(yC));
		    ySlope = (float)1.0;
		    G2 = magnitude(hDGImage[ptr_m1+j], vDGImage[ptr_m1+j]);
		    G4 = magnitude(hDGImage[ptr_p1+j], vDGImage[ptr_p1+j]);	
		    if((xC*yC) > (float)0.0){
			G1 = magnitude(hDGImage[ptr_m1+j-1], vDGImage[ptr_m1+j-1]);
			G3 = magnitude(hDGImage[ptr_p1+j+1], vDGImage[ptr_p1+j+1]);
		    }
		    else{
			G1 = magnitude(hDGImage[ptr_m1+j+1], vDGImage[ptr_m1+j+1]);
			G3 = magnitude(hDGImage[ptr_p1+j-1], vDGImage[ptr_p1+j-1]);
		    }
		}
		else{
		    /* horizontal gradient */
		    xSlope = (float)(fabs(yC) / fabs(xC));
		    ySlope = (float)1.0;
		    G2 = magnitude(hDGImage[ptr+j+1], vDGImage[ptr+j+1]);
		    G4 = magnitude(hDGImage[ptr+j-1], vDGImage[ptr+j-1]);	
		    if((xC*yC) > (float)0.0){
			G1 = magnitude(hDGImage[ptr_p1+j+1], vDGImage[ptr_p1+j+1]);
			G3 = magnitude(hDGImage[ptr_m1+j-1], vDGImage[ptr_m1+j-1]);
		    }
		    else{
			G1 = magnitude(hDGImage[ptr_m1+j+1], vDGImage[ptr_m1+j+1]);
			G3 = magnitude(hDGImage[ptr_p1+j-1], vDGImage[ptr_p1+j-1]);
		    }
		}
		if((G > (xSlope*G1+(ySlope-xSlope)*G2))&&(G > (xSlope*G3+(ySlope-xSlope)*G4))){
		    magImage[ptr+j] = G;	
		}
		if(magImage[ptr+j] > maxValue) maxValue = magImage[ptr+j];
		if(magImage[ptr+j] < minValue) minValue = magImage[ptr+j];
	    }
	}

	scale = (float)1.0 / (maxValue-minValue);
	ptr   = 0;
	count = 0;
	tAve  = 0.0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		magImage[ptr] = scale * (magImage[ptr]-minValue);
		if(magImage[ptr] > 0.0){
		    tAve += magImage[ptr];
		    ++count;
		}
		++ptr;
	    }
	}
	tAve /= (float)count;

	step = 255.0;
	for(i = 0; i < 256; ++i){
	    histogram[i] = 0;
	}
	ptr = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		value = (int)(step*(magImage[ptr]));
	        ++histogram[value];
		++ptr;
	    }
	}
	/*
	// now get the max after skipping the low values
	*/
	mValue = -1;
	mIndex = 0;
	for(i = 10; i < 256; ++i){
	    if(histogram[i] > mValue){
		mValue = histogram[i];
		mIndex = i;
	    }
	}

	if(mode == 1){
	    /* based on the mean value of edge energy */
	    *cannyLow  = ((*cannyLow)  * tAve);
	    *cannyHigh = ((*cannyHigh) * tAve);
	}
	else{
	    /* based on the mode value of edge energy */
	    *cannyLow  = ((*cannyLow)  * ((float)mIndex/step));
	    *cannyHigh = ((*cannyHigh) * ((float)mIndex/step));
	}

	return;

}

void DGFilters(int samples, int rows, int cols, double cannySigma, int gWidth,
               float *aveXValue, float *aveYValue, double *rawImage,
               double *dgKernel, float *hDGImage, float *vDGImage){

	/*
	// implements the derivative of Gaussian filter. kernel set by CannyEdges
	*/
	int i, j, k;
	int ptr;
	int mLength;
	int count;
	float *tBuffer = NULL;
	double sum;

	*aveXValue = (float)0.0;
	*aveYValue = (float)0.0;	

	mLength = MAX(rows, cols) + 64;
	tBuffer = calloc(mLength, sizeof(float));

	/*
	// filter X 
	*/
	count = 0;
	for(i = 0; i < rows; ++i){
	    ptr = i * cols;
	    for(j = gWidth; j < cols-gWidth; ++j){
		sum = dgKernel[0] * rawImage[ptr+j];
		for(k = 1; k < gWidth; ++k){
		    sum += dgKernel[k] * (-rawImage[ptr+j+k] + rawImage[ptr+j-k]);
		}
		hDGImage[ptr+j] = (float)sum;
		if(sum != (float)0.0){
		    ++count;
		    *aveXValue += (float)fabs(sum);
		}
	    }
	}
	if(count){
	    *aveXValue /= (float)count;
	    *aveXValue = (float)0.5 * (*aveXValue);
	    /* this is 50% of the max, hardwirred for now, and is part of the threshold */
	}
	/*
	// filter Y 
	*/
	count = 0;
	for(i = 0; i < cols; ++i){
	    for(j = 0; j < rows; ++j){
		ptr = j * cols;
		tBuffer[j] = rawImage[ptr+i];
	    }
	    for(j = gWidth; j < rows-gWidth; ++j){
		ptr = j * cols;
		sum = dgKernel[0] * tBuffer[j];
		for(k = 1; k < gWidth; ++k){
		    sum += dgKernel[k] * (-tBuffer[j+k] + tBuffer[j-k]);
		}
		vDGImage[ptr+i] = sum;
		if(sum != (float)0.0){
		    ++count;
		    *aveYValue += (float)fabs(sum);
		}
	    }
	}
	if(count){
	    *aveYValue /= (float)count;
	    *aveYValue = (float)0.5 * (*aveYValue);
	    /* this is 50% of the max, hardwirred for now, and is part of the threshold */
	}

	free(tBuffer);

	return;

}


int NI_CannyEdges(int samples, int rows, int cols, double cannySigma, 
                  double cannyLow, double cannyHigh, int mode, 
                  int lowThreshold, int highThreshold, double BPHigh,
                  int apearture, double *rawImage,
		  unsigned short *edgeImage, int *groups){

	int i, j;
	int offset;
	int doHysteresis = 0;
	int gWidth;
	int mLength;
	int status;
	float aveXValue;
	float aveYValue;
	double t;
	double dgKernel[20];
	float *HYSImage = NULL;
	float *hDGImage = NULL;
	float *vDGImage = NULL;
	float *magImage = NULL;
	float *tBuffer  = NULL;

	/* filter */
	doPreProcess(samples, rows, cols, rawImage, BPHigh, apearture, lowThreshold, highThreshold);

	/*
	// memory for magnitude, horizontal and vertical derivative of Gaussian filter
	*/
	mLength  = MAX(rows, cols) + 64;
	HYSImage = calloc(samples, sizeof(float));
	hDGImage = calloc(samples, sizeof(float));
	vDGImage = calloc(samples, sizeof(float));
	magImage = calloc(samples, sizeof(float));
	tBuffer  = calloc(mLength, sizeof(float));

	/*
	// build derivative of Gaussian filter kernel
	// kernel is anti-symmetric so convolution is k[j]*(v[i+j] - v[i-j]) 
	*/
	gWidth = 20;
	for(i = 0; i < gWidth; ++i){
	    t = (float)i;
	    dgKernel[i]  = (float)exp((double)((-t*t)/((float)2.0 * cannySigma * cannySigma)));
	    dgKernel[i] *= -(t / (cannySigma * cannySigma));
	}
	for(i = 0; i < samples; ++i){
	    HYSImage[i] = (float)0.0;
	}

	DGFilters(samples, rows, cols, cannySigma, gWidth, &aveXValue, &aveYValue,
	          rawImage, dgKernel, hDGImage, vDGImage); 
	nonMaxSupress(rows, cols, aveXValue, aveYValue, &cannyLow, &cannyHigh,
	              mode, hDGImage, vDGImage, magImage);
	if(doHysteresis){
	    edgeHysteresis(rows, cols, cannyLow, cannyHigh, magImage, HYSImage);
	}
	else{
	    edgeThreshold(rows, cols, cannyLow, magImage, HYSImage);
	}

	/*
	// edge image
	*/
	for(i = 0; i < samples; ++i){
	    edgeImage[i] = (unsigned short)HYSImage[i];
	}
	*groups = ConnectedEdgePoints(rows, cols, edgeImage);

	/*
	// prune the isolated pixels
	*/
	offset  = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(edgeImage[offset+j] > (*groups)){
		    edgeImage[offset+j] = 0;
		}	
	    }
	    offset  += cols;
	}


	free(tBuffer);
	free(hDGImage);
	free(vDGImage);
	free(magImage);
	free(HYSImage);

	status = *groups;
	return status;

}

void doSobel(int samples, int rows, int cols, double sobelLow, int mode, 
             double *rawImage, unsigned short *edgeImage){

	int i, j;
	int p, m, n;
	int offset;
	int offsetM1;
	int offsetP1;
	int minValue, maxValue;
	int pAve  = 0;
	int count = 0;
	int histogram[256];
	int value;
	int maxIndex;
	float pThreshold;
	double scale;
	double step;
	float *filteredImage = NULL;

	filteredImage = calloc(samples, sizeof(float));

	minValue = 10000;
	maxValue = -10000;

	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		filteredImage[offset+j] = 0;
		edgeImage[offset+j]     = 0;
	    }
	    offset += cols;
	}

	/*
	// Sobel
	*/
	offset = cols;
	for(i = 1; i < rows-1; ++i){
	    offsetM1 = offset - cols;
	    offsetP1 = offset + cols;
	    for(j = 1; j < cols-1; ++j){
	        n = 2*rawImage[offsetM1+j] + rawImage[offsetM1+j-1] + rawImage[offsetM1+j+1] -
	            2*rawImage[offsetP1+j] - rawImage[offsetP1+j-1] - rawImage[offsetP1+j+1];
	        m = 2*rawImage[offset+j-1] + rawImage[offsetM1+j-1] + rawImage[offsetP1+j-1] -
	            2*rawImage[offset+j+1] - rawImage[offsetM1+j+1] - rawImage[offsetP1+j+1];
	        p = (int)sqrt((float)(m*m) + (float)(n*n));
		if(p > 0){
		    pAve += p;
		    ++count;
		    if(p > maxValue) maxValue = p;
		    if(p < minValue) minValue = p;
		}
	        filteredImage[offset+j] = p;
	    }
	    offset += cols;
	}

	/* threshold based on ave */
	pAve /= count;
	scale = 1.0 / maxValue;

	step = 255.0/(maxValue-minValue);
	for(i = 0; i < 256; ++i){
	    histogram[i] = 0;
	}
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		value = (int)(step*(filteredImage[offset+j]-minValue));
	        ++histogram[value];
	    }
	    offset += cols;
	}
	/*
	// now get the max after skipping the low values
	*/
	maxValue = -1;
	maxIndex = 0;
	for(i = 10; i < 256; ++i){
	    if(histogram[i] > maxValue){
		maxValue = histogram[i];
		maxIndex = i;
	    }
	}

	if(mode == 1){
	    /* based on the mean value of edge energy */
	    pThreshold = (int)(sobelLow * (float)pAve);
	}
	else{
	    /* based on the mode value of edge energy */
	    pThreshold = (sobelLow * (minValue + ((float)maxIndex/step)));
	}

	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(filteredImage[offset+j] > pThreshold){
		    edgeImage[offset+j] = 1;
		}
		else{
		    edgeImage[offset+j] = 0;
		}
		filteredImage[offset+j] *= scale; 
	    }
	    offset += cols;
	}

	free(filteredImage);

	return;


}

void estimateThreshold(float *lowThreshold, float *highThreshold, float ShenCastanLow, 
                       int rows, int cols, float *SourceImage){

	int i, j;
	int offset;
	int value;
	int mIndex;
	int histogram[256];
	float low, high;
	float scale;

	low  = (float)1000.0;
	high = (float)-1000.0;

	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
	        if(fabs(SourceImage[offset+j]) > high) high = fabs(SourceImage[offset+j]);
	        if(fabs(SourceImage[offset+j]) < low)  low  = fabs(SourceImage[offset+j]);
	    }
	    offset += cols;
	}

	scale = (float)255.0 / (high-low);
	for(i = 0; i < 256; ++i){
	    histogram[i] = 0;
	}
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
	        value = (int)(scale*(fabs(SourceImage[offset+j]) - low)); 
	        ++histogram[value];
	    }
	    offset += cols;
	}

	/*
	// now get the edge energy mode
	*/
	value  = 0;
	mIndex = 10;
	for(i = 10; i < 256; ++i){
	    if(histogram[i] > value){
	        value  = histogram[i];
	        mIndex = i;
	    }
	}

	*highThreshold = ((float)mIndex / scale) + low;
	*lowThreshold  = ((float)mIndex / scale) + low;

	*highThreshold *= ShenCastanLow;
	*lowThreshold  *= ShenCastanLow;

	return;

}

void thresholdEdges(float *SourceImage, unsigned short *EdgeImage, double ShenCastanLow,
                    int rows, int cols){

	int i, j;
	int offset;
	float tLow, tHigh;

	/*
	// SourceImage contains the adaptive gradient
	// get threshold from the mode of the edge energy
	*/
	estimateThreshold(&tLow, &tHigh, ShenCastanLow, rows, cols, SourceImage);

	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(SourceImage[offset+j] > tLow){
		    EdgeImage[offset+j] = 1;
		}
		else{
		    EdgeImage[offset+j] = 0;
		}
	    }
	    offset += cols;
	}

	return;

}

float adaptiveGradient(float *BLImage, float *FilterImage, int nrow, int ncol, 
                       int cols, int window){

	int i, j;
	int offset;
	int numOn, numOff;
	int hWindow = window/2;
	float sumOn, sumOff;
	float aveOn, aveOff;

	numOn  = 0;
       	numOff = 0;

	sumOn  = (float)0.0;
       	sumOff = (float)0.0;

	aveOn  = (float)0.0;
       	aveOff = (float)0.0;

	offset = nrow * cols;
	for(i = -hWindow; i < hWindow; ++i){
	    for(j = -hWindow; j < hWindow; ++j){
		if(BLImage[offset+(i*cols)+(j+ncol)] == 1){
		    sumOn += FilterImage[offset+(i*cols)+(j+ncol)]; 
		    ++numOn;
		}
		else{
		    sumOff += FilterImage[offset+(i*cols)+(j+ncol)]; 
		    ++numOff;
		}
	    }
	}

	if(numOn){
	    aveOn = sumOn / numOn;
	}

	if(numOff){
	    aveOff = sumOff / numOff;
	}

	return (aveOff-aveOn);

}

void getZeroCrossings(float *SourceImage, float *FilterImage, float *BLImage, 
                      int rows, int cols, int window){

	int i, j;
	int offset;
	bool validEdge;

	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		SourceImage[offset+j] = 0.0; 
	    }
	    offset += cols;
	}

	offset = window*cols;
	for(i = window; i < rows-window; ++i){
	    for(j = window; j < cols-window; ++j){
		validEdge = FALSE;
		if((BLImage[offset+j] == 1) && (BLImage[offset+cols+j] == 0)){
		    if((FilterImage[offset+cols+j] - FilterImage[offset-cols+j]) > 0.0){
			validEdge = TRUE;
		    } 
		}
		else if((BLImage[offset+j] == 1) && (BLImage[offset+j+1] == 0)){
		    if((FilterImage[offset+j+1] - FilterImage[offset+j-1]) > 0.0){
			validEdge = TRUE;
		    } 
		}
		else if((BLImage[offset+j] == 1) && (BLImage[offset-cols+j] == 0)){
		    if((FilterImage[offset+cols+j] - FilterImage[offset-cols+j]) < 0.0){
			validEdge = TRUE;
		    } 
		}
		else if((BLImage[offset+j] == 1) && (BLImage[offset+j-1] == 0)){
		    if((FilterImage[offset+j+1] - FilterImage[offset+j-1]) < 0.0){
			validEdge = TRUE;
		    } 
		}
		if(validEdge){
		    /* adaptive gradeint is signed */
		    SourceImage[offset+j] = (float)fabs(adaptiveGradient(BLImage, FilterImage, i, j, cols, window));
		}
	    }
	    offset += cols;
	}

	return;

}


void computeBandedLaplacian(float *image1, float *image2, float *BLImage, int rows, int cols){

	int i, j;
	int offset;
	float t;

	/*
	// like an unsharp mask
	*/
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
	        t = image1[offset+j] - image2[offset+j];
		if(t < (float)0.0){
		    t = (float)0.0;
		}
		else{
		    t = (float)1.0;
		}
		BLImage[offset+j] = t;
	    }
	    offset += cols;
	}

	return;

}

void thresholdImage(float *Raw, float *Filtered, int rows, int cols, int tLow, int tHigh){

	int i, j;
	int ptr;

	ptr = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(Raw[ptr] > tHigh){
		    Raw[ptr]      = 0.0;
		    Filtered[ptr] = 0.0;
		}
		if(Raw[ptr] < tLow){
		    Raw[ptr]      = 0.0;
		    Filtered[ptr] = 0.0;
		}
		++ptr;
	    }
	}

	return;

}

void ISEF_Vertical(float *SourceImage, float *FilterImage, float *A, float *B, 
                   int rows, int cols, double b){


	int i, j;
	int offset;
	float b1, b2;

	b1 = ((float)1.0 - b)/((float)1.0 + b);
	b2 = b * b1;

	/*
	// set the boundaries
	*/
	offset = (rows-1)*cols;
	for(i = 0; i < cols; ++i){
	    /* process row 0 */
	    A[i] = b1 * SourceImage[i];
	    /* process row N-1 */
	    B[offset+i] = b2 * SourceImage[offset+i];
	}

	/*
	// causal component of IIR filter
	*/
	offset = cols;
	for(i = 1; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		/*
	        // IIR ISEF filter applied across rows
		*/
	        A[offset+j] = (b * A[offset-cols+j]) + (b1 * SourceImage[offset+j]);
	    }
	    offset += cols;
	}

	/*
	// anti-causal component of IIR filter
	*/
	offset = (rows-2)*cols;
	for(i = rows-2; i >= 0; --i){
	    for(j = 0; j < cols; ++j){
		/*
	        // IIR ISEF filter applied across rows
		*/
	        B[offset+j] = (b * B[offset+cols+j]) + (b2 * SourceImage[offset+j]); 
	    }
	    offset -= cols;
	}

	offset = (rows-1)*cols;
	for(j = 0; j < cols-1; ++j){
	    FilterImage[offset+j] = A[offset+j];
	}

	/*
	// add causal and anti-causal IIR parts
	*/
	offset = 0;
	for(i = 1; i < rows-2; ++i){
	    for(j = 0; j < cols-1; ++j){
	        FilterImage[offset+j] = A[offset+j] + B[offset+cols+j];
	    }
	    offset += cols;
	}

	return;

}

void ISEF_Horizontal(float *SourceImage, float *FilterImage, float *A, float *B,
                     int rows, int cols, double b){


	/*
	// source and smooth are the same in this pass of the 2D IIR
	*/

	int i, j;
	int offset;
	float b1, b2;

	b1 = ((float)1.0 - b)/((float)1.0 + b);
	b2 = b * b1;

	/*
	// columns boundaries
	*/
	offset = 0;
	for(i = 0; i < rows; ++i){
	    // col 0
	    A[offset] = b1 * SourceImage[offset];
	    // col N-1
	    B[offset+cols-1] = b2 * SourceImage[offset+cols-1];
	}

	/*
	// causal IIR part
	*/
	offset = 0;
	for(j = 1; j < cols; ++j){
	    for(i = 0; i < rows; ++i){
		A[offset+j] = (b * A[offset+j-1]) + (b1 * SourceImage[offset+j]);
	    }
	    offset += cols;
	}

	/*
	// anti-causal IIR part
	*/
	offset = 0;
	for(j = cols-2; j > 0; --j){
	    for(i = 0; i < rows; ++i){
		B[offset+j] = (b * B[offset+j+1]) + (b2 * SourceImage[offset+j]);
	    }
	    offset += cols;
	}

	/*
	// filtered output. this is 2-pass IIR and pass 1 is vertical
	*/
	offset = 0;
	for(i = 0; i < rows; ++i){
	    FilterImage[offset+cols-1] = A[offset+cols-1];
	}

	/*
	// add causal and anti-causal IIR parts
	*/
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols-1; ++j){
	        FilterImage[offset+j] = A[offset+j] + B[offset+j+1];
	    }
	    offset += cols;
	}

	return;

}


void computeISEF(float *SourceImage, float *FilterImage, int rows, int cols, double b){

	int imageSize = rows*cols;
	float *A;
	float *B;

	A = calloc(imageSize, sizeof(float));
	B = calloc(imageSize, sizeof(float));

	ISEF_Vertical(SourceImage, FilterImage, A, B, rows, cols, b);
	ISEF_Horizontal(FilterImage, FilterImage, A, B, rows, cols, b);

	free(A);
	free(B);

	return;

}

void Shen_Castan(double b, double ShenCastanLow, int rows, int cols, int window,
                 int lowThreshold, int highThreshold,
	       	 double *RawImage, unsigned short *EdgeImage){

	int i;
	int imageSize = rows*cols;
	float *FilterImage;
	float *BinaryLaplacianImage;
	float *SourceImage;

	FilterImage          = calloc(imageSize, sizeof(float));
	BinaryLaplacianImage = calloc(imageSize, sizeof(float));
	SourceImage          = calloc(imageSize, sizeof(float));

	for(i = 0; i < imageSize; ++i){
	    SourceImage[i] = RawImage[i];
	}
	computeISEF(SourceImage, FilterImage, rows, cols, b);
	/* optional thresholding based on low, high */
	thresholdImage(SourceImage, FilterImage, rows, cols, lowThreshold, highThreshold);
	computeBandedLaplacian(FilterImage, SourceImage, BinaryLaplacianImage, rows, cols);
	/* the new source image is now the adaptive gradient */
	getZeroCrossings(SourceImage, FilterImage, BinaryLaplacianImage, rows, cols, window);
	thresholdEdges(SourceImage, EdgeImage, ShenCastanLow, rows, cols);

	free(FilterImage);
	free(BinaryLaplacianImage);
	free(SourceImage);

	return;

}

int NI_ShenCastanEdges(int samples, int rows, int cols, double b, double ShenCastanLow,
                       int window, int lowThreshold, int highThreshold, 
                       double *rawImage, unsigned short *edgeImage, int *groups){


	int i, j;
	int offset;
	int status = 0;

	Shen_Castan(b, ShenCastanLow, rows, cols, window, lowThreshold, highThreshold, rawImage, edgeImage);
	*groups = ConnectedEdgePoints(rows, cols, edgeImage);


	//
	// prune the isolated pixels
	//
	offset  = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(edgeImage[offset+j] > (*groups)){
		    edgeImage[offset+j] = 0;
		}	
	    }
	    offset  += cols;
	}

	status = *groups;

	return status;

}

void buildBinaryImage(int rows, int cols, double *rawImage, unsigned short *edgeImage,
                      int lowThreshold, int highThreshold){

	int i, j;
	int offset;
	double value;
	int maskValue;

	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		value = rawImage[offset+j];
		maskValue = 1;
		if(value < (double)lowThreshold)  maskValue = 0;
		if(value > (double)highThreshold) maskValue = 0;
		edgeImage[offset+j] = maskValue;
	    }
	    offset += cols;
	}

	return;

}



void morphoFilterBinaryImage(int rows, int cols, unsigned short *edgeImage,
                             int CloseSize, int OpenSize){


	int i, j;
	int offset, offset2;
	unsigned short *cmask;
	unsigned short *omask;
	int olapValuesC[4];
	int olapValuesO[4];
	int CloseMaskSize = 1;
	int OpenMaskSize = 1;
	int LowValue1, HighValue1;   
	int LowValue2, HighValue2;  
	int spadSize;
	int maskSize = 11;
	unsigned char *ImageE;
	unsigned char *ImageC;

	spadSize = MAX(rows, cols);

	ImageE = calloc(spadSize*spadSize, sizeof(unsigned char));
	ImageC = calloc(spadSize*spadSize, sizeof(unsigned char));

	cmask = calloc(11*11, sizeof(unsigned short));
	omask = calloc(11*11, sizeof(unsigned short));

	//
	// Close filter
	//
	if(CloseSize){
	    CloseMaskSize = (CloseSize-1)/2;
	    for(i = 0; i < 2*CloseMaskSize+1; ++i){
	        for(j = 0; j < 2*CloseMaskSize+1; ++j){
	            cmask[i*maskSize+j] = 1;
	        }
	    }
	    LowValue1      = 0;   
	    HighValue1     = 1;   
	    LowValue2      = 1;   
	    HighValue2     = 0;   
	    olapValuesC[0] = LowValue1;
	    olapValuesC[1] = HighValue1;
	    olapValuesC[2] = LowValue2;
	    olapValuesC[3] = HighValue2;
	}

	/*
	// Open filter
	*/
	if(OpenSize){
	    OpenMaskSize = (OpenSize-1)/2;
	    for(i = 0; i < 2*OpenMaskSize+1; ++i){
	        for(j = 0; j < 2*OpenMaskSize+1; ++j){
	            omask[i*maskSize+j] = 1;
	        }
	    }
	    LowValue1      = 1;   
	    HighValue1     = 0;   
	    LowValue2      = 0;   
	    HighValue2     = 1;   
	    olapValuesO[0] = LowValue1;
	    olapValuesO[1] = HighValue1;
	    olapValuesO[2] = LowValue2;
	    olapValuesO[3] = HighValue2;
	}

	offset  = 0;
	offset2 = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		ImageE[offset2+j] = (unsigned char)edgeImage[offset+j]; 
	    }
	    offset2 += spadSize;
	    offset  += cols;
	}

	if(OpenSize){
	    OpenCloseFilter(olapValuesO, OpenMaskSize, rows, cols, spadSize, ImageE, ImageC, omask);
	}

	if(CloseSize){
	    OpenCloseFilter(olapValuesC, CloseMaskSize, rows, cols, spadSize, ImageE, ImageC, cmask);
	}

	offset  = 0;
	offset2 = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(ImageE[offset2+j] == 1){
		    /* this will activate some original off-pixels */
		    edgeImage[offset+j] = 1;
		}
		else{
		    /* this will zero some original on-pixels */
		    edgeImage[offset+j] = 0;
		}
	    }
	    offset2 += spadSize;
	    offset  += cols;
	}

	free(ImageE);
	free(ImageC);

	free(cmask);
	free(omask);

	return;

}

void doRegionGrow(int samples, int rows, int cols, double *rawImage,
                  unsigned short *edgeImage, int lowThreshold, 
		  int highThreshold, int closeWindow, int openWindow){

	buildBinaryImage(rows, cols, rawImage, edgeImage, lowThreshold, highThreshold);
	morphoFilterBinaryImage(rows, cols, edgeImage, closeWindow, openWindow);

	return;

}

int NI_RegionGrow(int samples, int rows, int cols, int lowThreshold, int highThreshold,
                  int closeWindow, int openWindow, double *rawImage, 
                  unsigned short *edgeImage, int *groups){

	int i, j;
	int offset;
	int status;

	doRegionGrow(samples, rows, cols, rawImage, edgeImage, lowThreshold,
	             highThreshold, closeWindow, openWindow);
	*groups = ConnectedEdgePoints(rows, cols, edgeImage);

	//
	// prune the isolated pixels
	//
	offset  = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(edgeImage[offset+j] > (*groups)){
		    edgeImage[offset+j] = 0;
		}	
	    }
	    offset  += cols;
	}

	status = *groups;
	return status;

}

int NI_SobelEdges(int samples, int rows, int cols, double sobelLow, int mode,
                  int lowThreshold, int highThreshold, double BPHigh,   
                  int apearture, double *rawImage, unsigned short *edgeImage, int *groups){


	int i, j;
	int offset;
	int status;

	doPreProcess(samples, rows, cols, rawImage, BPHigh, apearture, lowThreshold, highThreshold);
	doSobel(samples, rows, cols, sobelLow, mode, rawImage, edgeImage);
	*groups = ConnectedEdgePoints(rows, cols, edgeImage);
	
	
	/*
	// prune the isolated pixels
	*/
	offset  = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(edgeImage[offset+j] > (*groups)){
		    edgeImage[offset+j] = 0;
		}	
	    }
	    offset  += cols;
	}

	status = *groups;
	return status;

}

void initThinFilter(int *J_mask, int *K_mask){

	int i, j;
	int Column;
	int maskCols = 3;

	for(i = 0; i < 3; ++i){
	    for(j = 0; j < 30; ++j){
		J_mask[i+j*maskCols] = 0;
		K_mask[i+j*maskCols] = 0;
	    }
	}

	Column = 0;
   	J_mask[0+maskCols*(Column+0)] = 1;
   	J_mask[0+maskCols*(Column+1)] = 1;
   	J_mask[0+maskCols*(Column+2)] = 1;
   	J_mask[1+maskCols*(Column+1)] = 1;

	Column += 3;
   	J_mask[0+maskCols*(Column+1)] = 1;
   	J_mask[1+maskCols*(Column+1)] = 1;
   	J_mask[1+maskCols*(Column+2)] = 1;

	Column += 3;
   	J_mask[0+maskCols*(Column+0)] = 1;
   	J_mask[1+maskCols*(Column+0)] = 1;
   	J_mask[2+maskCols*(Column+0)] = 1;
   	J_mask[1+maskCols*(Column+1)] = 1;

	Column += 3;
   	J_mask[0+maskCols*(Column+1)] = 1;
   	J_mask[1+maskCols*(Column+0)] = 1;
   	J_mask[1+maskCols*(Column+1)] = 1;

	Column += 3;
   	J_mask[0+maskCols*(Column+2)] = 1;
   	J_mask[1+maskCols*(Column+1)] = 1;
   	J_mask[1+maskCols*(Column+2)] = 1;
   	J_mask[2+maskCols*(Column+2)] = 1;

	Column += 3;
   	J_mask[1+maskCols*(Column+0)] = 1;
   	J_mask[1+maskCols*(Column+1)] = 1;
   	J_mask[2+maskCols*(Column+1)] = 1;

	Column += 3;
   	J_mask[1+maskCols*(Column+1)] = 1;
   	J_mask[2+maskCols*(Column+0)] = 1;
   	J_mask[2+maskCols*(Column+1)] = 1;
   	J_mask[2+maskCols*(Column+2)] = 1;

	Column += 3;
   	J_mask[1+maskCols*(Column+1)] = 1;
   	J_mask[1+maskCols*(Column+2)] = 1;
   	J_mask[2+maskCols*(Column+1)] = 1;

	Column = 0;
   	K_mask[2+maskCols*(Column+0)] = 1;
   	K_mask[2+maskCols*(Column+1)] = 1;
   	K_mask[2+maskCols*(Column+2)] = 1;

	Column += 3;
   	K_mask[1+maskCols*(Column+0)] = 1;
   	K_mask[2+maskCols*(Column+0)] = 1;
   	K_mask[2+maskCols*(Column+1)] = 1;

	Column += 3;
   	K_mask[0+maskCols*(Column+2)] = 1;
   	K_mask[1+maskCols*(Column+2)] = 1;
   	K_mask[2+maskCols*(Column+2)] = 1;

	Column += 3;
   	K_mask[1+maskCols*(Column+2)] = 1;
   	K_mask[2+maskCols*(Column+1)] = 1;
   	K_mask[2+maskCols*(Column+2)] = 1;

	Column += 3;
   	K_mask[0+maskCols*(Column+0)] = 1;
   	K_mask[1+maskCols*(Column+0)] = 1;
   	K_mask[2+maskCols*(Column+0)] = 1;

	Column += 3;
   	K_mask[0+maskCols*(Column+1)] = 1;
   	K_mask[0+maskCols*(Column+2)] = 1;
   	K_mask[1+maskCols*(Column+2)] = 1;

	Column += 3;
   	K_mask[0+maskCols*(Column+0)] = 1;
   	K_mask[0+maskCols*(Column+1)] = 1;
   	K_mask[0+maskCols*(Column+2)] = 1;

	Column += 3;
   	K_mask[0+maskCols*(Column+0)] = 1;
   	K_mask[0+maskCols*(Column+1)] = 1;
   	K_mask[1+maskCols*(Column+0)] = 1;

	return;

}

void ThinningFilter(int regRows, int regColumns, int spadSize, int *J_mask, int *K_mask,
	            unsigned char *Input, unsigned char *CInput, unsigned char *ErosionStage,
	            unsigned char *DialationStage, unsigned char *HMT, unsigned char *Copy){

	int i, j, k, l, m, n, overlap, hit;
	int LowValue1, HighValue1;   
	int LowValue2, HighValue2;   
	int Column, T, nloop;
	int Offset;
	int N, M;
	int maskCols = 3;
	int j_mask[3][3];
	int k_mask[3][3];

	N = regRows;
	M = regColumns;

	LowValue1  = 1;   
	HighValue1 = 0;   

	LowValue2  = 0;   
	HighValue2 = 1;   

	Offset = 0;
	for(i = 0; i < N; ++i){
	    for(j = 0; j < M; ++j){
		Copy[Offset+j] = Input[Offset+j];
	    }
	    Offset += spadSize;
	}

	nloop = 0;
	while(1){
	    /* erode */
	    Column = 0;
	    for(n = 0; n < 8; ++n){
		for(i = 0; i < 3; ++i){
		    for(j = 0; j < 3; ++j){
			j_mask[i][j] = J_mask[i+maskCols*(Column+j)];
		    }
		}
		for(i = 0; i < 3; ++i){
		    for(j = 0; j < 3; ++j){
			k_mask[i][j] = K_mask[i+maskCols*(Column+j)];
		    }
		}
		Column += 3;

		Offset = spadSize;
		for(i = 1; i < N-1; ++i){
		    for(j = 1; j < M-1; ++j){
			hit = LowValue1; 
			for(k = -1; k < 2; ++k){
			    for(l = -1; l < 2; ++l){
				T = j_mask[k+1][l+1];
				if(T == 1){
				    overlap = T*Input[Offset+(k*spadSize)+j+l];
				    if(overlap == HighValue1) hit = HighValue1;
				}
			    }
			}
			ErosionStage[Offset+j] = hit;
		    }
		    Offset += spadSize;
		}

		/* dialate */
		Offset = 0;
		for(i = 0; i < N; ++i){
		    for(j = 0; j < M; ++j){
			CInput[Offset+j] = (~Input[Offset+j]) & 0x1; 
		    }
		    Offset += spadSize;
		}

		Offset = spadSize;
		for(i = 1; i < N-1; ++i){
		    for(j = 1; j < M-1; ++j){
			hit = LowValue1; 
			for(k = -1; k < 2; ++k){
			    for(l = -1; l < 2; ++l){
				T = k_mask[k+1][l+1];
				if(T == 1){
				    overlap = T*CInput[Offset+(k*spadSize)+j+l];
				    if(overlap == HighValue1) hit = HighValue1;
			        }
			    }
			}
			DialationStage[Offset+j] = hit;
		    }
		    Offset += spadSize;
		}

		/* form the HMT */
		Offset = 0;
		for(i = 0; i < N; ++i){
		    for(j = 0; j < M; ++j){
			m = (ErosionStage[Offset+j]*DialationStage[Offset+j]);
			HMT[Offset+j] = m;
		    }
		    Offset += spadSize;
		}

		/* Thin for stage n */

		Offset = 0;
		for(i = 0; i < N; ++i){
		    for(j = 0; j < M; ++j){
			HMT[Offset+j] = (~HMT[Offset+j]) & 0x1; 
		    }
		    Offset += spadSize;
		}

		Offset = 0;
		for (i = 0; i < N; ++i){
		    for (j = 0; j < M; ++j){
			m = (Input[Offset+j]*HMT[Offset+j]);
			Input[Offset+j] = m;
		    }
		    Offset += spadSize;
		}
	    }

	    /* check for no change */
	    hit = 0;
	    Offset = 0;
	    for(i = 0; i < N; ++i){
		for(j = 0; j < M; ++j){
		    hit += abs(Copy[Offset+j]-Input[Offset+j]);
		}
		Offset += spadSize;
	    }
	    if(!hit) break;

	    hit = 0;
	    Offset = 0;
	    for(i = 0; i < N; ++i){
		for(j = 0; j < M; ++j){
		    Copy[Offset+j] = Input[Offset+j];
		    if(Input[Offset+j]) ++hit;
		}
		Offset += spadSize;
	    }
	    /* nloop is data dependent. */
	    ++nloop;
	}


	return;

}


int NI_ThinFilter(int samples, int rows, int cols, int numberObjects,
                  unsigned short *edgeImage, objStruct objectMetrics[]){

	int i, j;
	int loop;
	int label;
	int left, right, top, bottom;
	int roiRows, roiCols;
	int srcOffset;
	int dstOffset;
	int status;
	int inflate = 1;
	int *J_mask;
	int *K_mask;

	unsigned char *Input;
	unsigned char *CInput;
	unsigned char *ErosionStage;
	unsigned char *DialationStage;
	unsigned char *HMT;
	unsigned char *Copy;
	unsigned short *thinEdgeImage;

	/*
	// scratch pad (spad) memory
	*/
	Input          = calloc(samples, sizeof(unsigned char));
	CInput         = calloc(samples, sizeof(unsigned char));
	ErosionStage   = calloc(samples, sizeof(unsigned char));
	DialationStage = calloc(samples, sizeof(unsigned char));
	HMT            = calloc(samples, sizeof(unsigned char));
	Copy           = calloc(samples, sizeof(unsigned char));
	thinEdgeImage  = calloc(samples, sizeof(unsigned short));
	J_mask         = calloc(3*30,    sizeof(int));
	K_mask         = calloc(3*30,    sizeof(int));

	initThinFilter(J_mask, K_mask);
	for(loop = 0; loop < numberObjects; ++loop){
	    label   = objectMetrics[loop].Label;
	    left    = objectMetrics[loop].L;
	    right   = objectMetrics[loop].R;
	    top     = objectMetrics[loop].T;
	    bottom  = objectMetrics[loop].B;
	    roiRows = top-bottom+2*inflate;
	    roiCols = right-left+2*inflate;

	    /*
	    // clear the scratch pad
	    */
	    srcOffset = 0;
	    for(i = 0; i < roiRows; ++i){
	        for(j = 0; j < roiCols; ++j){
		    Input[srcOffset+j] = 0; 
	        }
	        srcOffset += cols;
	    }

	    /*
	    // copy the ROI for MAT (medial axis transformation) filter
	    */
	    dstOffset = inflate*rows;
	    for(i = bottom; i < top; ++i){
		srcOffset = i*cols;
		for(j = left; j < right; ++j){
		    if(edgeImage[srcOffset+j] == label){
			Input[dstOffset+j-left+inflate] = 1;
		    }
		}
		dstOffset += cols;
	    }
	    ThinningFilter(roiRows, roiCols, cols, J_mask, K_mask, Input, CInput,
	                   ErosionStage, DialationStage, HMT, Copy);

	    /*
	    // copy the MAT roi to the new edgeImage (clip the inflate border)
	    */
	    dstOffset = inflate*rows;
	    for(i = bottom; i < top; ++i){
		srcOffset = i*cols;
		for(j = left; j < right; ++j){
		    if(Input[dstOffset+j-left+inflate]){
		        thinEdgeImage[srcOffset+j] = label;
		    }
		}
		dstOffset += cols;
	    }
	}

	/*
	// copy the MAT edges and return the thinned edges
	// this will prune the isolated edge points from the edgeImage source
	*/
	for(i = 0; i < rows*cols; ++i){
	    edgeImage[i] = thinEdgeImage[i];
	}

	free(Input);
	free(CInput);
	free(ErosionStage);
	free(DialationStage);
	free(HMT);
	free(Copy);
	free(thinEdgeImage);
	free(J_mask);
	free(K_mask);

	status = 1;

	return status;

}


void generateMask(unsigned char *ImageH, bPOINT *boundary, int newSamples, int label, int cols){

	/*
	// get the boundary point pairs (left, right) for each line
	// if there is no pair, then the boundary is open
	// then fill the image in with the current label
	*/

	int i, j, k, m;
	int list[2048];
	int distance;
	int neighbor = 4;
	int index;
	int offset;
	int maxDistance = 1024;
	int x, y;
	int low, high;
	
	for(i = 0; i < newSamples; ++i){
	    boundary[i].haveLink  = FALSE;
	    boundary[i].linkIndex = -1;
	}
	
	for(i = 0; i < newSamples; ++i){
	    if(!boundary[i].haveLink){
		boundary[i].haveLink = TRUE;
		x = boundary[i].x;
		y = boundary[i].y;
		for(k = 0, j = 0; j < newSamples; ++j){
		    if((j != i)){
			if(boundary[j].y == y){
			    list[k] = j;
			    ++k;
			}
		    }
		}
		/* now get the closest boundary */
		if(k){
		    distance = maxDistance;
		    index    = -1;
		    for(j = 0; j < k; ++j){
			m = abs(x - boundary[list[j]].x);
			if((m < distance) && (m > neighbor)){
			    distance = m;
			    index = list[j];
			}
			else if(m <= neighbor){
			    boundary[list[j]].haveLink = TRUE;
			}
		    }
		    if(index != -1){
			boundary[i].linkIndex     = index;
			boundary[index].linkIndex = i;
			boundary[index].haveLink  = TRUE;
			if(boundary[i].x < boundary[index].x){
			    low  = boundary[i].x;
			    high = boundary[index].x;
			}
			else{
			    low  = boundary[index].x;
			    high = boundary[i].x;
			}
			/*
			// do the fill
			*/
			offset = y * cols;
			for(j = low; j <= high; ++j){
			    ImageH[offset+j] = label;
			}
		    }
		}
		else{
		    /* boundary point is isolated */
		    boundary[i].linkIndex = i;
		}
	    }
	}

	return;

}

void getBoundaryMetrics(bPOINT *boundary, float *length, float *minRadius,
                        float *maxRadius, float *aveRadius,
	         	float Xcenter, float Ycenter, int newSamples){

	int j;
	float dX, dY;
	float distance;

	if(newSamples < 2){
	    *length    = (float)0.0;
	    *minRadius = (float)0.0;
	    *maxRadius = (float)0.0;
	    *aveRadius = (float)0.0;
	    return;
	}

	*length = (float)0.0;
	for(j = 1; j < newSamples; ++j){
	    dX = (float)(boundary[j].x - boundary[j-1].x);
	    dY = (float)(boundary[j].y - boundary[j-1].y);
	    distance = (float)sqrt(dX*dX + dY*dY);
	    *length += distance;
	}

	*minRadius = (float)10000.0;
	*maxRadius = (float)-10000.0;
	*aveRadius = (float)0.0;
	for(j = 0; j < newSamples; ++j){
	    dX = (float)(boundary[j].x - Xcenter);
	    dY = (float)(boundary[j].y - Ycenter);
	    distance = (float)sqrt(dX*dX + dY*dY);
	    *aveRadius += distance;
	    if(distance < *minRadius) *minRadius = distance;
	    if(distance > *maxRadius) *maxRadius = distance;
	}

	if(newSamples){
	    *aveRadius /= (float)newSamples;
	}

	return;

}

void trackBoundary(unsigned char *Input, blobBoundary lBoundary[], int mcount, int spadSize, 
		   blobBoundary seedValue, int searchWindow){


	int i, j, k, m, p;
	int offset;
	int CurI;
	int CurJ;
	int StrI;
	int StrJ;
	int NewI;
	int NewJ;
	int MinD;
	int inflate = searchWindow;

    	CurI = seedValue.xy.x;
    	CurJ = seedValue.xy.y;
    	StrI = CurI;
    	StrJ = CurJ;

	p = 0;
	lBoundary[p].xy.x = StrI;
	lBoundary[p].xy.y = StrJ;
	offset = StrI * spadSize;

	p = 1;
	while(p < mcount){
	    offset = (CurI-inflate)*spadSize;
	    MinD = 1024;
	    NewI = -1;
	    NewJ = -1;
	    for(i = CurI-inflate; i < CurI+inflate; ++i){
		for(j = CurJ-inflate; j < CurJ+inflate; ++j){
		    m = Input[offset+j];
		    if(m == 1){
			/* city block distance */
			k = abs(i-CurI) + abs(j-CurJ);
			if(k < MinD){
			    MinD = k;
			    NewI = i;
			    NewJ = j;
			}
		    }
		}
		offset += spadSize;
	    }
	    if(NewI != -1) CurI = NewI;
	    if(NewJ != -1) CurJ = NewJ;
	    offset = CurI * spadSize;
	    Input[offset+CurJ] = 0;
	    lBoundary[p].xy.x = CurJ;
	    lBoundary[p].xy.y = CurI;
            ++p;
	}

	return;

}


void OpenCloseFilter(int olapValues[], int maskSize, int rows, int columns, int spadSize, 
                     unsigned char *input, unsigned char *output, unsigned short *mask){


	/*
	// do morphological open/close image filtering. the olapValues array determines
    	// if the filter is Open or Close. 
	*/
	int i, j, k, l, m, overlap, hit;
	int offset;
	int LowValue1, HighValue1;   
	int LowValue2, HighValue2;  
	int morphoMaskSize = 11;

	LowValue1  = olapValues[0];
	HighValue1 = olapValues[1];
	LowValue2  = olapValues[2];
	HighValue2 = olapValues[3];

	/* close - step 1 is dialate 
	   open  - step 1 is erode */
	offset = maskSize*spadSize;
	for(i = maskSize; i < rows-maskSize; ++i){
	    for(j = maskSize; j < columns-maskSize; ++j){
	        hit = LowValue1; 
		for(k = -maskSize; k < maskSize; ++k){
	    	    m = k*spadSize;
		    for(l = -maskSize; l < maskSize; ++l){
	    		overlap = mask[morphoMaskSize*(k+maskSize)+(l+maskSize)]*input[offset+m+j+l];
			if(overlap == HighValue1){
			    hit = HighValue1;
			}
		    }
		}
	    	output[offset+j] = hit;
	    }
	    offset += spadSize;
	}

	/* close - step 2 is erode
	   open -  step 2 is dialate */
	offset = maskSize*spadSize;
	for(i = maskSize; i < rows-maskSize; ++i){
	    for(j = maskSize; j < columns-maskSize; ++j){
	        hit = LowValue2; 
		for(k = -maskSize; k < maskSize; ++k){
	    	    m = k*spadSize;
		    for(l = -maskSize; l < maskSize; ++l){
	    		overlap = mask[morphoMaskSize*(k+maskSize)+(l+maskSize)]*output[offset+m+j+l];
			if(overlap == HighValue2){
			    hit = HighValue2;
			}
		    }
		}
	    	input[offset+j] = hit;
	    }
	    offset += spadSize;
	}

	return;
}

void getCompactness(unsigned char *Input, RECT roi, int label, int spadSize,
                    float *vCompactness, float length){

	int i, j;
	int maskOffset;
	int area;
	static float fpi = (float)(4.0 * 3.14159);

	area = 0;
	for(i = roi.bottom; i < roi.top; ++i){
	    maskOffset = i*spadSize;
	    for(j = roi.left; j < roi.right; ++j){
		if(Input[maskOffset+j] == label){
		    ++area;
		}
	    }
	}
	if(area && (length != (float)0.0)){
	    *vCompactness = (fpi * (float)area) / (length*length);
	}
	else{
	    *vCompactness = (float)0.0;
	}

	return;
}


void doMorphology(unsigned char *Input, unsigned char *ImageE, unsigned char *ImageC,
                  unsigned char *ImageH, int olapValuesC[], int olapValuesO[], 
       	          unsigned short *cmask, unsigned short *omask,
	          RECT roi, int label, int CloseMaskSize, int OpenMaskSize, int spadSize){

	int i, j;
	int rows, cols;
	int srcOffset;
	int dstOffset;
	int maskSize;

	cols = roi.right - roi.left;
	rows = roi.top - roi.bottom;

	for(i = 0; i < spadSize*spadSize; ++i){
	    ImageE[i] = 0;
	    ImageC[i] = 0;
	}

	/*
	// put the ROI in the ImageE array centered in ULC
	*/
	dstOffset = 0;
	for(i = roi.bottom; i < roi.top; ++i){
	    srcOffset = i*spadSize;
	    for(j = roi.left; j < roi.right; ++j){
		if(ImageH[srcOffset+j] == label){
		    ImageE[dstOffset+j-roi.left] = 1;
		}
	    }
	    dstOffset += spadSize;
	}

	/*
	// open
	*/
	maskSize = OpenMaskSize;
	OpenCloseFilter(olapValuesO, maskSize, rows, cols, spadSize, ImageE, ImageC, omask);
	/*
	// close
	*/
	maskSize = CloseMaskSize;
	OpenCloseFilter(olapValuesC, maskSize, rows, cols, spadSize, ImageE, ImageC, cmask);

	/*
	// put the closed ROI (in ImageE) back in its roi space
	*/

	srcOffset = 0;
	for(i = roi.bottom; i < roi.top+2*maskSize+1; ++i){
	    dstOffset = (i-(2*maskSize+1))*spadSize;
	    for(j = roi.left-maskSize-1; j < roi.right+maskSize+1; ++j){
		if(ImageE[srcOffset+j-roi.left] == 1){
		    Input[dstOffset+j-maskSize+1] = label;
		}
	    }
	    srcOffset += spadSize;
	}

	return;

}


void getBoundary(unsigned short *ThinEdgeImage, unsigned char *Input,
                 blobBoundary *pBoundary, blobBoundary *lBoundary, 
	         boundaryIndex *pBoundaryIndex, RECT boundBox, int label,
	         int bBox, int nextSlot, int memOffset,
		 int spadSize, int searchWindow){

	int i, j;
	int dstOffset;
	int srcOffset;
	int mcount;
	int rows;
	int columns;
	bool first;
	blobBoundary value;
	int inflate = searchWindow+1;
	int count;

	pBoundaryIndex[bBox+1].rectangle.left   = boundBox.left;
	pBoundaryIndex[bBox+1].rectangle.right  = boundBox.right;
	pBoundaryIndex[bBox+1].rectangle.top    = boundBox.top;
	pBoundaryIndex[bBox+1].rectangle.bottom = boundBox.bottom;

	for(i = 0; i < spadSize*spadSize; ++i){
	    Input[i] = 0;
	}

	/* copy to spad */

	count = 0;
	rows    = boundBox.top-boundBox.bottom+2*inflate;
	columns = boundBox.right-boundBox.left+2*inflate;
	dstOffset = inflate*spadSize;
	for(i = boundBox.bottom; i < boundBox.top; ++i){
	    srcOffset = i*spadSize;
	    for(j = boundBox.left; j < boundBox.right; ++j){
		if(ThinEdgeImage[srcOffset+j] == label){
		    Input[dstOffset+j-boundBox.left+inflate] = 1;
		    ++count;
		}
	    }
	    dstOffset += spadSize;
	}

	mcount    = 0;
	first     = TRUE;
	srcOffset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < columns; ++j){
		if(Input[srcOffset+j]){
		    if(first){
			first = FALSE;
			/* index of the seed sample */
			value.xy.x = i;
			value.xy.y = j;
		    }
		    ++mcount;
		}
	    }
	    srcOffset += spadSize;
	}

	trackBoundary(Input, lBoundary, mcount, spadSize, value, searchWindow);	

	pBoundaryIndex[nextSlot].numberPoints = mcount;
	for(i = 0; i < mcount; ++i){
	    value.xy.x = lBoundary[i].xy.x + boundBox.left   - inflate;
	    value.xy.y = lBoundary[i].xy.y + boundBox.bottom - inflate + 1;
	    pBoundary[memOffset].xy.x = value.xy.x;
	    pBoundary[memOffset].xy.y = value.xy.y;
	    ++memOffset;
	}

	return;

}


void buildBoundary(objStruct objectMetrics[], int searchWindow, unsigned short *ThinEdgeImage,
		   int numberObjects, int srcRows, int srcCols){

	int i, j, k;
	int count;
	int numBoundaries;
	int numSamples;
	int offset;
	int offset2;
	int end;
	int label;
	int distance;
	/* these will be user-setup parameters */
	int closureDistance = 12;
	int CloseSize       = 5;
	int OpenSize        = 5;
	int threshold       = 3;
	int newSamples;
	int spadSize;
	POINT rectPoint[4];
	int in[4];
	float length;
	float minRadius;
	float maxRadius;
	float aveRadius;
	float vCompactness;
	/* for morphological close of mask. max structuring element is 11x11 */
	unsigned short *cmask;
	unsigned short *omask;
	int maskSize = 11;
	int olapValuesC[4];
	int olapValuesO[4];
	int CloseMaskSize;
	int OpenMaskSize;
	int LowValue1, HighValue1;   
	int LowValue2, HighValue2;  
	RECT bBox;

	boundaryIndex *pBoundaryIndex;
	blobBoundary  *pBoundary;
	blobBoundary  *lBoundary;
	bPOINT        *boundary;
	unsigned char *Input;
	unsigned char *ImageE;
	unsigned char *ImageC;
	unsigned char *ImageH;

	spadSize = srcCols;
	pBoundaryIndex = calloc(srcRows+srcCols,   sizeof(boundaryIndex));
	Input          = calloc(spadSize*spadSize, sizeof(unsigned char));
	ImageE         = calloc(spadSize*spadSize, sizeof(unsigned char));
	ImageC         = calloc(spadSize*spadSize, sizeof(unsigned char));
	ImageH         = calloc(spadSize*spadSize, sizeof(unsigned char));
	pBoundary      = calloc(srcRows*srcCols,   sizeof(blobBoundary));
	lBoundary      = calloc(32767, sizeof(blobBoundary));
	boundary       = calloc(32767, sizeof(POINT));
	cmask          = calloc(11*11, sizeof(unsigned short));
	omask          = calloc(11*11, sizeof(unsigned short));

	/*
	// Close filter
	*/
	CloseMaskSize = (CloseSize-1)/2;
	for(i = 0; i < 2*CloseMaskSize+1; ++i){
	    for(j = 0; j < 2*CloseMaskSize+1; ++j){
	        cmask[i*maskSize+j] = 1;
	    }
	}
	LowValue1      = 0;   
	HighValue1     = 1;   
	LowValue2      = 1;   
	HighValue2     = 0;   
	olapValuesC[0] = LowValue1;
	olapValuesC[1] = HighValue1;
	olapValuesC[2] = LowValue2;
	olapValuesC[3] = HighValue2;

	/*
	// Open filter
	*/
	OpenMaskSize = (OpenSize-1)/2;
	for(i = 0; i < 2*OpenMaskSize+1; ++i){
	    for(j = 0; j < 2*OpenMaskSize+1; ++j){
	        omask[i*maskSize+j] = 1;
	    }
	}
	LowValue1      = 1;   
	HighValue1     = 0;   
	LowValue2      = 0;   
	HighValue2     = 1;   
	olapValuesO[0] = LowValue1;
	olapValuesO[1] = HighValue1;
	olapValuesO[2] = LowValue2;
	olapValuesO[3] = HighValue2;

	for(i = 0; i < (srcRows+srcCols); ++i){
	    pBoundaryIndex[i].numberPoints = 0;
	    pBoundaryIndex[i].curveClose   = 0;
	    pBoundaryIndex[i].isWithin     = FALSE;
	    pBoundaryIndex[i].criticalSize = FALSE;
	    pBoundaryIndex[i].closedCurve  = FALSE;
	}


	for(i = 0; i < numberObjects; ++i){
	    ++pBoundaryIndex[0].numberPoints;
	    count = 0;
	    j = 1;
	    while(pBoundaryIndex[j].numberPoints){
		count += pBoundaryIndex[j++].numberPoints;
	    }
	    bBox.left   = objectMetrics[i].L;
	    bBox.right  = objectMetrics[i].R;
	    bBox.top    = objectMetrics[i].T;
	    bBox.bottom = objectMetrics[i].B;
	    label       = objectMetrics[i].Label;
	    pBoundaryIndex[i+1].Label = label;
	    getBoundary(ThinEdgeImage, Input, pBoundary, lBoundary, pBoundaryIndex, bBox, label,
		        i, pBoundaryIndex[0].numberPoints, count, spadSize, searchWindow);
	}

	/*
	// Input will now be used in the fill. Copy the labeled edge image
	*/

	offset = 0;
	numBoundaries = pBoundaryIndex[0].numberPoints;
	for(i = 0; i < numBoundaries; ++i){
	    numSamples = pBoundaryIndex[i+1].numberPoints;
	    end        = numSamples-2; 
	    newSamples = numSamples-1;
	    for(j = 0; j < numSamples; ++j){
		boundary[j].x = pBoundary[offset+j+1].xy.x;
		boundary[j].y = pBoundary[offset+j+1].xy.y;
	    }

	    /*
	    // clip off the ends where stray boundary pixels were left over
	    */
	    while(1){
		distance = abs(boundary[end].x-boundary[end-1].x) + abs(boundary[end].y-boundary[end-1].y);
		if(distance > threshold){
		    --end;
		    --newSamples;
		}
		else{
		    break;
		}
	    }

	    distance = abs(boundary[0].x-boundary[end-2].x) + abs(boundary[0].y-boundary[end-2].y);
	    pBoundaryIndex[i+1].curveClose = distance;

	    if(pBoundaryIndex[i+1].curveClose < closureDistance){
		pBoundaryIndex[i+1].closedCurve = TRUE;
	    }
	    pBoundaryIndex[i+1].centroid.x = 0;
	    pBoundaryIndex[i+1].centroid.y = 0;
	    for(j = 0; j < newSamples; ++j){
	        pBoundaryIndex[i+1].centroid.x += boundary[j].x;
	        pBoundaryIndex[i+1].centroid.y += boundary[j].y;
	    }
	    if(newSamples){
	        pBoundaryIndex[i+1].centroid.x /= newSamples;
	        pBoundaryIndex[i+1].centroid.y /= newSamples;
	    }
	    getBoundaryMetrics(boundary, &length, &minRadius, &maxRadius, &aveRadius,
		       	      (float)pBoundaryIndex[i+1].centroid.x,
		       	      (float)pBoundaryIndex[i+1].centroid.y, newSamples);
	    pBoundaryIndex[i+1].boundaryLength = length;
	    pBoundaryIndex[i+1].minRadius      = minRadius;
	    pBoundaryIndex[i+1].maxRadius      = maxRadius;
	    pBoundaryIndex[i+1].aveRadius      = aveRadius;
	    if(minRadius != 0.0){
	        pBoundaryIndex[i+1].ratio = maxRadius / minRadius;
	    }
	    else{
	        pBoundaryIndex[i+1].ratio = -1.0;
	    }

	    /*
	    // augment the ROI boundary
	    */
	    pBoundaryIndex[i+1].rectangle.left   -= 2*CloseMaskSize;
	    pBoundaryIndex[i+1].rectangle.right  += 2*CloseMaskSize;
	    pBoundaryIndex[i+1].rectangle.bottom -= 2*CloseMaskSize;
	    pBoundaryIndex[i+1].rectangle.top    += 2*CloseMaskSize;
	    label = pBoundaryIndex[i+1].Label;

	    /*
	    // mask goes in ImageH. morpho filter the mask first
	    */
	    generateMask(ImageH, boundary, newSamples, label, spadSize);

	    /*
	    // open-close the mask 
	    */
	    doMorphology(Input, ImageE, ImageC, ImageH, olapValuesC, olapValuesO, cmask, omask,
		         pBoundaryIndex[i+1].rectangle, label, CloseMaskSize, OpenMaskSize, spadSize);

	    /*
	    // now get the compactness metrics
	    */
	    getCompactness(Input, pBoundaryIndex[i+1].rectangle, label, spadSize, &vCompactness, length);
	    pBoundaryIndex[i+1].compactness = vCompactness;

	    /*
	    // reset the ROI boundary
	    */
	    pBoundaryIndex[i+1].rectangle.left   += 2*CloseMaskSize;
	    pBoundaryIndex[i+1].rectangle.right  -= 2*CloseMaskSize;
	    pBoundaryIndex[i+1].rectangle.bottom += 2*CloseMaskSize;
	    pBoundaryIndex[i+1].rectangle.top    -= 2*CloseMaskSize;
	    offset += numSamples;
	}
	

	for(i = 0; i < numBoundaries; ++i){
	    for(j = 0; j < numBoundaries; ++j){
		if(j != i){
		    rectPoint[0].x = pBoundaryIndex[j+1].rectangle.left;
		    rectPoint[0].y = pBoundaryIndex[j+1].rectangle.bottom;
		    rectPoint[1].x = pBoundaryIndex[j+1].rectangle.left;
		    rectPoint[1].y = pBoundaryIndex[j+1].rectangle.top;
		    rectPoint[2].x = pBoundaryIndex[j+1].rectangle.right;
		    rectPoint[2].y = pBoundaryIndex[j+1].rectangle.bottom;
		    rectPoint[3].x = pBoundaryIndex[j+1].rectangle.right;
		    rectPoint[3].y = pBoundaryIndex[j+1].rectangle.top;
		    in[0] = 0;
		    in[1] = 0;
		    in[2] = 0;
		    in[3] = 0;
		    for(k = 0; k < 4; ++k){
			if((rectPoint[k].x > pBoundaryIndex[i+1].rectangle.left) &&
			   (rectPoint[k].x < pBoundaryIndex[i+1].rectangle.right)){
			    if((rectPoint[k].y > pBoundaryIndex[i+1].rectangle.bottom) &&
			       (rectPoint[k].y < pBoundaryIndex[i+1].rectangle.top)){
				in[k] = 1;
			    }
			}
		    }
		    if(in[0] && in[1] && in[2] && in[3]){
			pBoundaryIndex[j+1].isWithin = TRUE;
		    }
		}
	    }
	}

	/*
	// fill in the Python features
	*/
	for(i = 0; i < numBoundaries; ++i){
	    objectMetrics[i].curveClose     = pBoundaryIndex[i+1].curveClose;
	    objectMetrics[i].cXBoundary     = pBoundaryIndex[i+1].centroid.x;
	    objectMetrics[i].cYBoundary     = pBoundaryIndex[i+1].centroid.y;
	    objectMetrics[i].boundaryLength = pBoundaryIndex[i+1].boundaryLength;
	    objectMetrics[i].minRadius      = pBoundaryIndex[i+1].minRadius;
	    objectMetrics[i].maxRadius      = pBoundaryIndex[i+1].maxRadius;
	    objectMetrics[i].aveRadius      = pBoundaryIndex[i+1].aveRadius;
	    objectMetrics[i].ratio          = pBoundaryIndex[i+1].ratio;
	    objectMetrics[i].compactness    = pBoundaryIndex[i+1].compactness;
	} 

	// debug only
	if(0){
	for(i = 0; i < numBoundaries; ++i){
	    if(pBoundaryIndex[i+1].boundaryLength != (float)0.0){
	        printf("boundary %d:\n", i);
	        printf("\t\tRect (%d, %d, %d, %d)\n", pBoundaryIndex[i+1].rectangle.left,
	                                              pBoundaryIndex[i+1].rectangle.right,
	                                              pBoundaryIndex[i+1].rectangle.top,
	                                              pBoundaryIndex[i+1].rectangle.bottom);
	        printf("\t\tCentroid (%d, %d)\n",     pBoundaryIndex[i+1].centroid.x, pBoundaryIndex[i+1].centroid.y);
	        printf("\t\tLength (%f)\n",           pBoundaryIndex[i+1].boundaryLength);
	        printf("\t\tRatio (%f)\n",            pBoundaryIndex[i+1].ratio);
	        printf("\t\taveRadius (%f)\n",        pBoundaryIndex[i+1].aveRadius);
	        printf("\t\tLabel (%d)\n",            pBoundaryIndex[i+1].Label);
	        printf("\t\tCompactness (%f)\n",      pBoundaryIndex[i+1].compactness);
	        printf("\t\tCurveClose (%d)\n",       pBoundaryIndex[i+1].curveClose);
	        if(pBoundaryIndex[i+1].isWithin){
	            printf("\t\tContained (T)\n");
	        }
	        else{
	            printf("\t\tContained (F)\n");
	        }
	        if(pBoundaryIndex[i+1].closedCurve){
	            printf("\t\tclosedCurve (T)\n");
	        }
	        else{
	            printf("\t\tclosedCurve (F)\n");
	        }
	    }
	}
	}

	/*
	// need to return input which is now mask image
	*/

	offset  = 0;
	offset2 = 0;
	for(i = 0; i < srcRows; ++i){
	    for(j = 0; j < srcCols; ++j){
	        ThinEdgeImage[offset+j] = (unsigned short)Input[offset2+j];
	    }
	    offset  += srcCols;
	    offset2 += spadSize;
	}

	free(pBoundaryIndex);
	free(Input);
	free(ImageE);
	free(ImageC);
	free(ImageH);
	free(pBoundary);
	free(lBoundary);
	free(boundary);
	free(cmask);
	free(omask);

	return;

}


void initLaws(LawsFilter7 *lawsFilter){

	int i;
	float sum;
	float L7[7] = { 1.0,  6.0,  15.0, 20.0,  15.0,  6.0,  1.0};
	float E7[7] = {-1.0, -4.0,  -5.0,  0.0,   5.0,  4.0,  1.0};
	float S7[7] = {-1.0, -2.0,   1.0,  4.0,   1.0, -2.0, -1.0};
	float W7[7] = {-1.0,  0.0,   3.0,  0.0,  -3.0,  0.0,  1.0};
	float R7[7] = { 1.0, -2.0,  -1.0,  4.0,  -1.0, -2.0,  1.0};
	float O7[7] = {-1.0,  6.0, -15.0, 20.0, -15.0,  6.0, -1.0};
	
	lawsFilter->numberKernels      = 6;
	lawsFilter->kernelLength       = 7;
	lawsFilter->numberFilterLayers = 21;
	lawsFilter->name[0] = 'L';
	lawsFilter->name[1] = 'E';
	lawsFilter->name[2] = 'S';
	lawsFilter->name[3] = 'W';
	lawsFilter->name[4] = 'R';
	lawsFilter->name[5] = 'O';
	for(i = 0; i < 7; ++i){
	    lawsFilter->lawsKernel[0][i] = L7[i];
	    lawsFilter->lawsKernel[1][i] = E7[i];
	    lawsFilter->lawsKernel[2][i] = S7[i];
	    lawsFilter->lawsKernel[3][i] = W7[i];
	    lawsFilter->lawsKernel[4][i] = R7[i];
	    lawsFilter->lawsKernel[5][i] = O7[i];
	}

	/* L filter is unity gain */
	sum = (float)0.0;
	for(i = 0; i < 7; ++i){
	    sum += lawsFilter->lawsKernel[0][i];
	}
	for(i = 0; i < 7; ++i){
	    lawsFilter->lawsKernel[0][i] /= sum;
	}
	
	return;

}

float lawsConvolution(float *image, float *rowFilter, float *colFilter, int kernelSize){

	int i, j;
	int offset;
	float result[7];
	float sum;

	/* filter rows */
	for(i = 0; i < kernelSize; ++i){
	    sum = (float)0.0;
	    offset = i * kernelSize;
	    for(j = 0; j < kernelSize; ++j){
		sum += (rowFilter[j]*image[offset+j]);
	    }
	    result[i] = sum;
	}

	/* filter columns */
	sum = (float)0.0;
	for(j = 0; j < kernelSize; ++j){
	    sum += (rowFilter[j]*result[j]);
	}

	return(sum);

}


void getLawsTexture(LawsFilter7 lawsFilter, tTEM LawsFeatures[],
                    objStruct objectMetrics[], double *sourceImage, 
	            unsigned short *MaskImage, int numberObjects,
	            int srcRows, int srcCols){

	int i, j;
	int label;
	RECT bBox;
	int aperature = (lawsFilter.kernelLength-1)/2;
	unsigned char *ImageH;
	float *ImageT;
	float *lawsImage;

	ImageH    = calloc(srcRows*srcCols, sizeof(unsigned char));
	ImageT    = calloc(srcRows*srcCols, sizeof(float));
	lawsImage = calloc(lawsFilter.numberFilterLayers*srcRows*srcCols, sizeof(float));
	
	for(i = 0; i < numberObjects; ++i){
	    bBox.left   = objectMetrics[i].L;
	    bBox.right  = objectMetrics[i].R;
	    bBox.top    = objectMetrics[i].T;
	    bBox.bottom = objectMetrics[i].B;
	    label       = objectMetrics[i].Label;
	    if(objectMetrics[i].voxelMean != (float)0.0){
		/*
		// valid size region
		*/
	        computeLaws(lawsFilter, LawsFeatures, bBox, label, aperature, srcRows, srcCols, ImageH, ImageT,
		            MaskImage, lawsImage, sourceImage);
		for(j = 1; j < lawsFilter.numberFilterLayers; ++j){
		    objectMetrics[i].TEM[j-1] = LawsFeatures[j].Variance;
		}
	        /* -- later will need to return a view of the texture images
		int index;
		int offset;
		int layerStep = srcRows*srcCols;
	        if(label == debugBlob){ 
		    index = 0;
		    for(j = 1; j < lawsFilter.numberFilterLayers; ++j){
		        if(LawsFeatures[j].Variance == (float)1.0) index = j;
		    }
		    // overwrite the raw image
		    offset = index * layerStep;
		    for(j = 0; j < layerStep; ++j){
		        sourceImage[j] = lawsImage[offset+j];
	            }
	        }
		*/
	    }
	}

	free(ImageH);
	free(ImageT);
	free(lawsImage);

	return;

}

void computeLaws(LawsFilter7 lawsFilter, tTEM LawsFeatures[], RECT roi, int label,
                 int aperature, int srcRows, int srcCols, 
	         unsigned char *ImageH, float *ImageT, unsigned short *MaskImage,
	         float *lawsImage, double *sourceImage){

	/*
	// hard-wirred to Law's 7 kernels
	*/
	int i, j, k;
	int lawsLayer;
	int column, row;
	int offset;
	int maskOffset[7];
	int dataOffset[7];
	float myImage[49];
	int count;
	int outerKernelNumber;
	int innerKernelNumber;
	int rowNumber;
	int kernelSize = lawsFilter.kernelLength;
	int fullMask   = kernelSize*kernelSize;
	int layerStep  = srcRows*srcCols;
	float *rowFilter;
	float *colFilter;
	float filterResult1;
	float filterResult2;
	float lawsLL=1.0;
	float t;
	float maxValue;
	float scale;
	char I, J;
	char combo[24];
	char dual[24];


	/* zero the laws mask memory first */
	for(i = 0; i < srcRows*srcCols; ++i){
	    ImageH[i] = 0;
	}
	for(j = 0; j < lawsFilter.numberFilterLayers; ++j){
	    LawsFeatures[j].Mean     = (float)0.0;
	    LawsFeatures[j].Variance = (float)0.0;
	}

	for(i = roi.bottom+aperature; i < roi.top-aperature; ++i){
	    // get the row array offset for mask and data source. 
	    for(row = -aperature; row <= aperature; ++row){
		maskOffset[row+aperature] = (i+row)*srcCols;
		dataOffset[row+aperature] = maskOffset[row+aperature];
	    }
	    for(j = roi.left+aperature; j < roi.right-aperature; ++j){
		/*
		// get 7x7 segment and make sure have 100% mask coverage
		*/
		count = 0;
		for(row = -aperature; row <= aperature; ++row){
		    rowNumber = (row+aperature)*kernelSize;
		    for(column = -aperature; column <= aperature; ++column){
			if(MaskImage[maskOffset[row+aperature]+j+column] == label){
			    myImage[rowNumber+column+aperature] = sourceImage[dataOffset[row+aperature]+j+column];
			    ++count;
			}
		    }
		}
		if(count == fullMask){
		    /*
		    // 100% coverage. now do the Law's texture filters
		    */
		    ImageH[i*srcCols+j] = 1;
		    lawsLayer = 0;
		    for(outerKernelNumber = 0; outerKernelNumber < lawsFilter.numberKernels; ++outerKernelNumber){
			/*
			// outer loop pulls the i'th kernel. kernel 0 is the LP kernel
			// the outer loop is the iso-kernel
			*/
			I = lawsFilter.name[outerKernelNumber];
			sprintf(dual, "%c_%c", I, I);
			rowFilter = &lawsFilter.lawsKernel[outerKernelNumber][0];
			colFilter = &lawsFilter.lawsKernel[outerKernelNumber][0];
			filterResult1 = lawsConvolution(myImage, rowFilter, colFilter, kernelSize);
			/* lawsLayer 0 is the LP and needs to be used to scale. */
			if(outerKernelNumber){
			    lawsImage[lawsLayer*layerStep + i*srcCols + j] = (float)2.0 * filterResult1 / lawsLL;
			}
			else{
			    lawsLL = (float)2.0 * filterResult1;
			    lawsImage[lawsLayer*layerStep + i*srcCols + j] = (float)2.0 * filterResult1;
			}
			strcpy(&LawsFeatures[lawsLayer].filterName[0], dual);
			++lawsLayer;
			/*
			// now do the inner loop and get the column filters for the other laws kernels
			*/
			for(innerKernelNumber = outerKernelNumber+1;
			                        innerKernelNumber < lawsFilter.numberKernels;
			                        ++innerKernelNumber){
			    J = lawsFilter.name[innerKernelNumber];
			    sprintf(combo, "%c_%c", I, J);
			    strcpy(&LawsFeatures[lawsLayer].filterName[0], combo);
			    colFilter = &lawsFilter.lawsKernel[innerKernelNumber][0];
			    filterResult1 = lawsConvolution(myImage, rowFilter, colFilter, kernelSize);
			    filterResult2 = lawsConvolution(myImage, colFilter, rowFilter, kernelSize);
			    lawsImage[lawsLayer*layerStep + i*srcCols + j] =
			                        (filterResult1 / lawsLL) + (filterResult2 / lawsLL);
			    ++lawsLayer;
			}
		    }
		}
	    }
	}

	for(i = 0; i < lawsFilter.numberFilterLayers; ++i){
	    LawsFeatures[i].Mean     = (float)0.0;
	    LawsFeatures[i].Variance = (float)0.0;
	}

	count = 0;
	for(i = roi.bottom+aperature; i < roi.top-aperature; ++i){
	    row = i * srcCols;
	    for(j = roi.left+aperature; j < roi.right-aperature; ++j){
		if(ImageH[row+j]){
		    ++count;
		    for(k = 0; k < lawsFilter.numberFilterLayers; ++k){
			offset = k * layerStep + row;
			LawsFeatures[k].Mean += lawsImage[offset+j];
		    }
	        }
	    }
	}

	if(count == 0){
	    // debug statement
	    printf("no samples for texture\n");
	    return;
	}

	for(k = 0; k < lawsFilter.numberFilterLayers; ++k){
	    LawsFeatures[k].Mean /= (float)count;
	}
	for(i = roi.bottom+aperature; i < roi.top-aperature; ++i){
	    row = i * srcCols;
	    for(j = roi.left+aperature; j < roi.right-aperature; ++j){
		if(ImageH[row+j]){
		    for(k = 0; k < lawsFilter.numberFilterLayers; ++k){
			offset = k * layerStep + row;
			t = lawsImage[offset+j] - LawsFeatures[k].Mean;
			LawsFeatures[k].Variance += (t * t);
		    }
		}
	    }
	}
	for(k = 0; k < lawsFilter.numberFilterLayers; ++k){
	    LawsFeatures[k].Variance /= (float)count;
	    LawsFeatures[k].Variance = (float)(sqrt(LawsFeatures[k].Variance));
	}

	/*
	// now normalize the variance feature (TEM)
	*/
	maxValue = (float)0.0;
	for(i = 1; i < lawsFilter.numberFilterLayers; ++i){
	    if((LawsFeatures[i].Variance) > maxValue) maxValue = LawsFeatures[i].Variance;
	}
	scale = (float)1.0 / maxValue;
	for(i = 1; i < lawsFilter.numberFilterLayers; ++i){
	    LawsFeatures[i].Variance = scale * LawsFeatures[i].Variance;
	}


	return;

}

void getVoxelMeasures(objStruct objectMetrics[], double *sourceImage,
                      unsigned short *MaskImage, int numberObjects, 
		      int srcRows, int srcCols){

	int i, j, k;
	int label;
	int offset;
	int count;
	float mean, std, t;
	RECT bBox;

	for(i = 0; i < numberObjects; ++i){
	    bBox.left   = objectMetrics[i].L;
	    bBox.right  = objectMetrics[i].R;
	    bBox.top    = objectMetrics[i].T;
	    bBox.bottom = objectMetrics[i].B;
	    label       = objectMetrics[i].Label;
	    count = 0;
	    mean  = (float)0.0;
	    for(j = bBox.bottom; j < bBox.top; ++j){
	        offset = j * srcCols;
	        for(k = bBox.left; k < bBox.right; ++k){
		    if(MaskImage[offset+k] == label){
	    		mean += sourceImage[offset+k];
			++count;
		    }
	        }
	    }
	    if(count){
	    	mean /= (float)count; 
	        std = (float)0.0;
	        for(j = bBox.bottom; j < bBox.top; ++j){
	            offset = j * srcCols;
	            for(k = bBox.left; k < bBox.right; ++k){
		        if(MaskImage[offset+k] == label){
	    		    t = (sourceImage[offset+k]-mean);
			    std += (t * t);
		        }
	            }
	        }
	    }
	    if(count){
	        std /= (float)count; 
	        std = sqrt(std);
	        objectMetrics[i].voxelMean = mean;
	        objectMetrics[i].voxelVar  = std;
	    }
	    else{
	        objectMetrics[i].voxelMean = 0.0;
	        objectMetrics[i].voxelVar  = 0.0;
	    }
	}

	return;

}

int NI_BuildBoundary(int samples, int rows, int cols, int numberObjects, 
	             unsigned short *edgeImage, objStruct objectMetrics[]){

	int searchWindow = 5;  // 5 is good value for Sobel
	int status = 1;

	buildBoundary(objectMetrics, searchWindow, edgeImage, numberObjects, rows, cols);

	return status;

}

int NI_VoxelMeasures(int samples, int rows, int cols, int numberObjects, double *sourceImage,
	             unsigned short *maskImage, objStruct objectMetrics[]){

	int status = 1;
	getVoxelMeasures(objectMetrics, sourceImage, maskImage, numberObjects, rows, cols);

	return status;

}


int NI_TextureMeasures(int samples, int rows, int cols, int numberObjects, double *sourceImage,
	               unsigned short *maskImage, objStruct objectMetrics[]){

	int status = 1;
	LawsFilter7 lawsFilter;
	tTEM LawsFeatures[21];

	initLaws(&lawsFilter);
	getLawsTexture(lawsFilter, LawsFeatures, objectMetrics, sourceImage,
	               maskImage, numberObjects, rows, cols);

	return status;

}



