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


int NI_EdgePreFilter(int num, int rows, int cols, int lowThreshold, int highThreshold,
                     int aperature, int HalfFilterTaps, unsigned short *sImage, double *dImage,
		     double *kernel){

	int i, j, k, n, num1;
    	int offset;
	double sum, value;
	double *buffer;
	int max_buffer = MAX(rows, cols);
	int status;

	buffer = calloc(max_buffer+aperature+16, sizeof(double));

	num1 = HalfFilterTaps;
	offset = 0;
	for(i = 0; i < rows; ++i){
	    /* copy image row to local buffer  */
	    for(j = 0; j < cols; ++j){
		buffer[num1+j] = sImage[offset+j];
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
	        dImage[offset+n] = sum;
	    }
	    offset += cols;
	}

	offset = 0;
	for(i = 0; i < cols; ++i){
	    /* copy image column to local buffer */
	    offset = 0;
	    for(j = 0; j < rows; ++j){
            buffer[num1+j] = dImage[offset+i];
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
	        dImage[offset+i] = sum;
	        offset += cols;
	    }
	}

	/* threshold the image */
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		value = dImage[offset+j];
		if(value < (float)lowThreshold)  value = (float)0.0;
		if(value > (float)highThreshold) value = (float)0.0;
		dImage[offset+j] = value;
	    }
	    offset += cols;
	}

	free(buffer);

	status = 1;

	return(status);

}

int NI_SobelImage(int samples, int rows, int cols, double *rawImage, double *edgeImage, double *pAve,
	          int *minValue, int *maxValue){
             
	int i, j;
	int p, m, n;
	int offset;
	int offsetM1;
	int offsetP1;
	int status;
	int count = 0;

	/*
	// Sobel
	*/
	offset = cols;
	*pAve = 0.0;
	*minValue = 10000;
	*maxValue = -10000;
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
		    *pAve += p;
		    if(p > *maxValue) *maxValue = p;
		    if(p < *minValue) *minValue = p;
		    ++count;
		}
	        edgeImage[offset+j] = p;
	    }
	    offset += cols;
	}
	/* threshold based on ave */
	*pAve /= count;

	status = 1;

	return(status);

}


int NI_BinaryEdge(int samples, int rows, int cols, unsigned short *labelImage, unsigned short *edgeImage){ 

	int i, j, k;
	int maxValue;
	int offset;
	int offsetM1;
	int offsetP1;
	int values3x3[8];
	int status;

	offset = cols;
	for(i = 1; i < rows-1; ++i){
	    offsetM1 = offset - cols;
	    offsetP1 = offset + cols;
	    for(j = 1; j < cols-1; ++j){
		values3x3[0] = labelImage[offset+j] - labelImage[offset+j+1];
		values3x3[1] = labelImage[offset+j] - labelImage[offsetM1+j+1];
		values3x3[2] = labelImage[offset+j] - labelImage[offsetM1+j];
		values3x3[3] = labelImage[offset+j] - labelImage[offsetM1+j-1];
		values3x3[4] = labelImage[offset+j] - labelImage[offset+j-1];
		values3x3[5] = labelImage[offset+j] - labelImage[offsetP1+j-1];
		values3x3[6] = labelImage[offset+j] - labelImage[offsetP1+j];
		values3x3[7] = labelImage[offset+j] - labelImage[offsetP1+j+1];
		maxValue = -1;
		for(k = 0; k < 8; ++k){
		    maxValue = MAX(maxValue, values3x3[k]);
		}
	        edgeImage[offset+j] = maxValue;
	    }
	    offset += cols;
	}

	status = 1;

	return(status);

}

int NI_SobelEdge(int samples, int rows, int cols, double *edgeImage, unsigned short *edges, 
	         int mode, double pAve, int minValue, int maxValue, double sobelLow){

	int i, j;
	int offset;
	int value;
	int maxIndex;
	int status;
	int histogram[256];
	float pThreshold;
	double scale;
	double step;

	scale = 1.0 / maxValue;

	step = 255.0/(maxValue-minValue);
	for(i = 0; i < 256; ++i){
	    histogram[i] = 0;
	}
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		value = (int)(step*(edgeImage[offset+j]-minValue));
	        ++histogram[value];
	    }
	    offset += cols;
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
		if(edgeImage[offset+j] > pThreshold){
		    edges[offset+j] = 1;
		}
		else{
		    edges[offset+j] = 0;
		}
	    }
	    offset += cols;
	}

	status = 1;
	return(status);

}

int NI_GetBlobs3D(int samples, int layers, int rows, int cols, unsigned short *edges,
	       	   unsigned short *connectedEdges, int *groups, int mask){ 

	int  i, j, k, l, m;
	int  lOffset, rOffset, Label;
	int  lOffsetP, lOffsetN;
	int  rOffsetP, rOffsetN;
	int  Classes[4096];
	int  dwImageSize, ptr;
	bool NewLabel;
	bool Change;
	bool connected;
	int  T[27];
	int  *ccompImage;
	int  layerSize;
	int  count;
	int  status;

	layerSize   = rows * cols;
	dwImageSize = layers * rows * cols;
	ccompImage  = calloc(dwImageSize, sizeof(int ));

	Label = 1;
	for(i = 1; i < layers-1; ++i){
	    lOffset  = i * layerSize;
	    lOffsetP = lOffset+layerSize;
	    lOffsetN = lOffset-layerSize;
	    for(j = 1; j < rows-1; ++j){
		rOffset = j * cols;
		rOffsetP = rOffset+cols;
		rOffsetN = rOffset-cols;
		for(k = 1; k < cols-1; ++k){
		    if(edges[lOffset+rOffset+k]){
			/*
			 check 3x3x3 connectivity
			*/

			T[0]  = edges[lOffset+rOffset+k];
			T[1]  = edges[lOffset+rOffset+k+1];
			T[2]  = edges[lOffset+rOffsetN+k+1];
			T[3]  = edges[lOffset+rOffsetN+k];
			T[4]  = edges[lOffset+rOffsetN+k-1];
			T[5]  = edges[lOffset+rOffset+k-1];
			T[6]  = edges[lOffset+rOffsetP+k-1];
			T[7]  = edges[lOffset+rOffsetP+k];
			T[8]  = edges[lOffset+rOffsetP+k+1];

			T[9]  = edges[lOffsetN+rOffset+k];
			T[10] = edges[lOffsetN+rOffset+k+1];
			T[11] = edges[lOffsetN+rOffsetN+k+1];
			T[12] = edges[lOffsetN+rOffsetN+k];
			T[13] = edges[lOffsetN+rOffsetN+k-1];
			T[14] = edges[lOffsetN+rOffset+k-1];
			T[15] = edges[lOffsetN+rOffsetP+k-1];
			T[16] = edges[lOffsetN+rOffsetP+k];
			T[17] = edges[lOffsetN+rOffsetP+k+1];

			T[18] = edges[lOffsetP+rOffset+k];
			T[19] = edges[lOffsetP+rOffset+k+1];
			T[20] = edges[lOffsetP+rOffsetN+k+1];
			T[21] = edges[lOffsetP+rOffsetN+k];
			T[22] = edges[lOffsetP+rOffsetN+k-1];
			T[23] = edges[lOffsetP+rOffset+k-1];
			T[24] = edges[lOffsetP+rOffsetP+k-1];
			T[25] = edges[lOffsetP+rOffsetP+k];
			T[26] = edges[lOffsetP+rOffsetP+k+1];

			connected = FALSE;
			if(mask == 1){
			    count = 0;
			    for(l = 1; l < 27; ++l){
				count += T[l];
			    }
			    if(count){
				connected = TRUE;
			    }
			}
			else if(mask == 6){
			    count = (T[2] + T[4] + T[6] + T[8] + T[9] + T[18]);
			    if(count == 6){
				connected = TRUE;
			    }
			}
			else if(mask == 14){
			    count = (T[2] + T[4] + T[6] + T[8] + T[9] + T[18] + T[11] +
				     T[13] + T[15] + T[17] + T[20] + T[22] + T[24] + T[26]);
			    if(count == 14){
				connected = TRUE;
			    }
			}
			else if(mask == 26){
			    count = 0;
			    for(l = 1; l < 27; ++l){
				count += T[l];
			    }
			    if(count == 26){
				connected = TRUE;
			    }
			}
			if(connected){
			    ccompImage[lOffset+rOffset+k] = Label++;
			}
		    }
		}
	    }
	}


	while(1){
	Change = FALSE;
	    /*
	    // TOP-DOWN Pass for labeling
	    */
	    for(i = 1; i < layers-1; ++i){
		lOffset  = i * layerSize;
		lOffsetP = lOffset+layerSize;
		lOffsetN = lOffset-layerSize;
		for(j = 1; j < rows-1; ++j){
		    rOffset = j * cols;
		    rOffsetP = rOffset+cols;
		    rOffsetN = rOffset-cols;
		    for(k = 1; k < cols-1; ++k){
			if(ccompImage[lOffset+rOffset+k] != 0){

			    T[0]  = ccompImage[lOffset+rOffset+k];
			    T[1]  = ccompImage[lOffset+rOffset+k+1];
			    T[2]  = ccompImage[lOffset+rOffsetN+k+1];
			    T[3]  = ccompImage[lOffset+rOffsetN+k];
			    T[4]  = ccompImage[lOffset+rOffsetN+k-1];
			    T[5]  = ccompImage[lOffset+rOffset+k-1];
			    T[6]  = ccompImage[lOffset+rOffsetP+k-1];
			    T[7]  = ccompImage[lOffset+rOffsetP+k];
			    T[8]  = ccompImage[lOffset+rOffsetP+k+1];

			    T[9]  = ccompImage[lOffsetN+rOffset+k];
			    T[10] = ccompImage[lOffsetN+rOffset+k+1];
			    T[11] = ccompImage[lOffsetN+rOffsetN+k+1];
			    T[12] = ccompImage[lOffsetN+rOffsetN+k];
			    T[13] = ccompImage[lOffsetN+rOffsetN+k-1];
			    T[14] = ccompImage[lOffsetN+rOffset+k-1];
			    T[15] = ccompImage[lOffsetN+rOffsetP+k-1];
			    T[16] = ccompImage[lOffsetN+rOffsetP+k];
			    T[17] = ccompImage[lOffsetN+rOffsetP+k+1];

			    T[18] = ccompImage[lOffsetP+rOffset+k];
			    T[19] = ccompImage[lOffsetP+rOffset+k+1];
			    T[20] = ccompImage[lOffsetP+rOffsetN+k+1];
			    T[21] = ccompImage[lOffsetP+rOffsetN+k];
			    T[22] = ccompImage[lOffsetP+rOffsetN+k-1];
			    T[23] = ccompImage[lOffsetP+rOffset+k-1];
			    T[24] = ccompImage[lOffsetP+rOffsetP+k-1];
			    T[25] = ccompImage[lOffsetP+rOffsetP+k];
			    T[26] = ccompImage[lOffsetP+rOffsetP+k+1];
						
			    m = T[0];
			    for(l = 1; l < 27; ++l){
				if(T[l] != 0){
				    if(T[l] < m) m = T[l];
			        }
			    }
			    if(m != ccompImage[lOffset+rOffset+k]){
				Change = TRUE;
				ccompImage[lOffset+rOffset+k] = m;
			    }
			}
		    }
		}
	    }
	    /*
	    // BOTTOM-UP Pass for labeling
	    */
	    for(i = layers-1; i > 0; --i){
		lOffset  = i * layerSize;
		lOffsetP = lOffset+layerSize;
		lOffsetN = lOffset-layerSize;
		for(j = rows-1; j > 0; --j){
		    rOffset = j * cols;
		    rOffsetP = rOffset+cols;
		    rOffsetN = rOffset-cols;
		    for(k = cols-1; k > 0; --k){
			if(ccompImage[lOffset+rOffset+k] != 0){

			    T[0]  = ccompImage[lOffset+rOffset+k];
			    T[1]  = ccompImage[lOffset+rOffset+k+1];
			    T[2]  = ccompImage[lOffset+rOffsetN+k+1];
			    T[3]  = ccompImage[lOffset+rOffsetN+k];
			    T[4]  = ccompImage[lOffset+rOffsetN+k-1];
			    T[5]  = ccompImage[lOffset+rOffset+k-1];
			    T[6]  = ccompImage[lOffset+rOffsetP+k-1];
			    T[7]  = ccompImage[lOffset+rOffsetP+k];
			    T[8]  = ccompImage[lOffset+rOffsetP+k+1];

			    T[9]  = ccompImage[lOffsetN+rOffset+k];
			    T[10] = ccompImage[lOffsetN+rOffset+k+1];
			    T[11] = ccompImage[lOffsetN+rOffsetN+k+1];
			    T[12] = ccompImage[lOffsetN+rOffsetN+k];
			    T[13] = ccompImage[lOffsetN+rOffsetN+k-1];
			    T[14] = ccompImage[lOffsetN+rOffset+k-1];
			    T[15] = ccompImage[lOffsetN+rOffsetP+k-1];
			    T[16] = ccompImage[lOffsetN+rOffsetP+k];
			    T[17] = ccompImage[lOffsetN+rOffsetP+k+1];

			    T[18] = ccompImage[lOffsetP+rOffset+k];
			    T[19] = ccompImage[lOffsetP+rOffset+k+1];
			    T[20] = ccompImage[lOffsetP+rOffsetN+k+1];
			    T[21] = ccompImage[lOffsetP+rOffsetN+k];
			    T[22] = ccompImage[lOffsetP+rOffsetN+k-1];
			    T[23] = ccompImage[lOffsetP+rOffset+k-1];
			    T[24] = ccompImage[lOffsetP+rOffsetP+k-1];
			    T[25] = ccompImage[lOffsetP+rOffsetP+k];
			    T[26] = ccompImage[lOffsetP+rOffsetP+k+1];
						
			    m = T[0];
			    for(l = 1; l < 27; ++l){
				if(T[l] != 0){
				    if(T[l] < m) m = T[l];
			        }
			    }
			    if(m != ccompImage[lOffset+rOffset+k]){
				Change = TRUE;
				ccompImage[lOffset+rOffset+k] = m;
			    }
			}
		    }
		}
	    }

	    if(!Change) break;

	}   /* end while loop  */

	Label      = 1;
	Classes[0] = 0;
	ptr        = 0;
	for(i = 0; i < layers; ++i){
	    for(j = 0; j < rows; ++j){
		for(k = 0; k < cols; ++k){
		    m =	ccompImage[ptr];
		    ++ptr;
		    if(m > 0){
			NewLabel = TRUE;
			for(l = 1; l < Label; ++l){
			    if(Classes[l] == m) NewLabel = FALSE;
			}
			if(NewLabel){
			    Classes[Label++] = m;
			    if(Label > 4000){
				return 0;
			    }
			}
		    }
		}
	    }
	}

	*groups = Label;

	ptr = 0;
	for(i = 0; i < layers; ++i){
	    for(j = 0; j < rows; ++j){
		for(k = 0; k < cols; ++k){
		    m =	ccompImage[ptr];
		    for(l = 1; l < Label; ++l){
			if(Classes[l] == m){
			    connectedEdges[ptr] = l;
			    break;
			}
		    }
		    ++ptr;
		}
	    }
	}

	free(ccompImage);

	status = 1;
	return(status);

}

int NI_GetBlobs2D(int samples, int rows, int cols, unsigned short *edges, unsigned short *connectedEdges,
	       	  int *groups, int mask){ 

	int            i, j, k, l, m;
	int            offset;
	int            Label;
	int            status;
	int            Classes[4096];
	bool           NewLabel;
	bool           Change;
	bool           connected;
	int            count;
	unsigned short T[12];

	/*
	// connected components labeling. pixels with 1, 4 or 8 connectedness. 
	*/
	Label  = 1;
	offset = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		connectedEdges[offset+j] = 0; 
		if(edges[offset+j] == 1){
		    connected = FALSE;
		    if(mask == 1){
			count = 0;
			for(l = 1; l < 9; ++l){
			    count += T[l];
			}
			if(count){
			    connected = TRUE;
			}
		    }
		    else if(mask == 4){
			count = (T[2] + T[4] + T[6] + T[8]);
			if(count == 4){
			    connected = TRUE;
			}
		    }
		    else if(mask == 8){
			count = 0;
			for(l = 1; l < 9; ++l){
			    count += T[l];
			}
			if(count == 8){
			    connected = TRUE;
			}
		    }
		    if(connected){
		        connectedEdges[offset+j] = Label++; 
		    }
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

	*groups = Label;

	/*
	// prune the isolated pixels
	*/
	offset  = 0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		if(connectedEdges[offset+j] > (*groups)){
		    connectedEdges[offset+j] = 0;
		}	
	    }
	    offset  += cols;
	}

	status = 1;
	return(status);

}

int NI_GetBlobRegions3D(int layers, int rows, int cols, int numberObjects, 
                        unsigned short *labeledEdges, objStruct objectMetrics[]){


	int status;
	int i, j, k, l, m;
	int offset;
	int count;
	int LowX;
	int LowY;
	int LowZ;
	int HighX;
	int HighY;
	int HighZ;
	int ptr;
	float centerX;
	float centerY;
	float centerZ;

	for(l = 1; l < numberObjects; ++l){
	    offset     = cols;
	    LowX       = 32767;
	    LowY       = 32767;
	    LowZ       = 32767;
	    HighX      = 0;
	    HighY      = 0;
	    HighZ      = 0;
	    count      = 0;
	    centerX    = (float)0.0;
	    centerY    = (float)0.0;
	    centerZ    = (float)0.0;
	    ptr        = 0;
	    for(i = 0; i < layers; ++i){
		for(j = 0; j < rows; ++j){
		    for(k = 0; k < cols; ++k){
		        m = labeledEdges[ptr++];
			if(l == m){
			    if(i < LowZ)   LowZ = i;
			    if(j < LowY)   LowY = j;
			    if(k < LowX)   LowX = k;
			    if(i > HighZ) HighZ = i;
			    if(j > HighY) HighY = j;
			    if(k > HighX) HighX = k;
	    		    centerX += (float)k;
	    		    centerY += (float)j;
	    		    centerZ += (float)i;
	    		    ++count;
			}
		    }
		}
	    }
	    /* the bounding box for the 2D blob */
	    objectMetrics[l-1].Left   = LowX;
	    objectMetrics[l-1].Right  = HighX;
	    objectMetrics[l-1].Bottom = LowY;
	    objectMetrics[l-1].Top    = HighY;
	    objectMetrics[l-1].Front  = LowZ;
	    objectMetrics[l-1].Back   = HighZ;
	    objectMetrics[l-1].Mass   = count;
	    objectMetrics[l-1].cX     = centerX/(float)count;
	    objectMetrics[l-1].cY     = centerY/(float)count;
	    objectMetrics[l-1].cZ     = centerZ/(float)count;
	    objectMetrics[l-1].Label  = l;
	}

	status = numberObjects;

	return(status);

}

int NI_GetBlobRegions2D(int rows, int cols, int numberObjects, unsigned short *labeledEdges,
                        objStruct objectMetrics[]){

	int i, j, k, m;
	int count;
	int LowX;
	int LowY;
	int HighX;
	int HighY;
	int status;
	int ptr;
	float centerX;
	float centerY;

	for(k = 1; k < numberObjects; ++k){
	    LowX       = 32767;
	    LowY       = 32767;
	    HighX      = 0;
	    HighY      = 0;
	    count      = 0;
	    centerX    = (float)0.0;
	    centerY    = (float)0.0;
	    ptr        = 0;
	    for(i = 0; i < rows; ++i){
		for(j = 0; j < cols; ++j){
		    m = labeledEdges[ptr++];
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
	    }
	    /* the bounding box for the 2D blob */
	    objectMetrics[k-1].Left   = LowX;
	    objectMetrics[k-1].Right  = HighX;
	    objectMetrics[k-1].Bottom = LowY;
	    objectMetrics[k-1].Top    = HighY;
	    objectMetrics[k-1].Mass   = count;
	    objectMetrics[k-1].cX     = centerX/(float)count;
	    objectMetrics[k-1].cY     = centerY/(float)count;
	    objectMetrics[k-1].Label  = k;
	}

	status = numberObjects;
	return status;

}


int NI_ThinMorphoFilter(int regRows, int regColumns, int spadSize, int masks, unsigned short *J_mask, 
	                 unsigned short *K_mask, unsigned char *Input, unsigned char *CInput, 
	                 unsigned char *ErosionStage, unsigned char *DialationStage, 
		         unsigned char *HMT, unsigned char *Copy){

	int i, j, k, l, m, n, overlap, hit;
	int LowValue1, HighValue1;   
	int LowValue2, HighValue2;   
	int Column, T, nloop;
	int Offset;
	int N, M;
	int maskCols = 3;
	int j_mask[3][3];
	int k_mask[3][3];
	int status;

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
	    for(n = 0; n < masks; ++n){
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


	status = 1;
	return status;

}


int NI_CannyFilter(int samples, int rows, int cols, double *rawImage,
		   double *hDGImage, double *vDGImage, double *dgKernel, 
                   int gWidth, float *aveXValue, float *aveYValue){
               

	/*
	// implements the derivative of Gaussian filter. kernel set by CannyEdges
	*/
	int i, j, k;
	int ptr;
	int mLength;
	int count;
	int status;
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
	}

	free(tBuffer);

	status = 1;

	return status;

}

double tmagnitude(double X, double Y){
	return sqrt(X*X + Y*Y);
}

int NI_CannyNonMaxSupress(int num, int rows, int cols, double *magImage, double *hDGImage,
	                  double *vDGImage, int mode, double aveXValue, double aveYValue,
			  double *tAve, double *cannyLow, double *cannyHigh, 
			  double cannyL, double cannyH){
                   
	int i, j;
	int ptr, ptr_m1, ptr_p1;
	float xSlope, ySlope, G1, G2, G3, G4, G, xC, yC;
	float scale;
	float maxValue = (float)0.0;
	float minValue = (float)0.0;
	int value;
	int mValue;
	int mIndex;
	int count;
	int status;
	int histogram[256];
	double step;

	for(i = 1; i < rows-1; ++i){
	    ptr = i * cols;
	    ptr_m1 = ptr - cols;
	    ptr_p1 = ptr + cols;
	    for(j = 1; j < cols; ++j){
		magImage[ptr+j] = (float)0.0;
		xC = hDGImage[ptr+j];
		yC = vDGImage[ptr+j];
		if(!((fabs(xC) < aveXValue) && (fabs(yC) < aveYValue))){
		    G = tmagnitude(xC, yC);
		    if(fabs(yC) > fabs(xC)){
		        /* vertical gradient */
		        xSlope = (float)(fabs(xC) / fabs(yC));
		        ySlope = (float)1.0;
		        G2 = tmagnitude(hDGImage[ptr_m1+j], vDGImage[ptr_m1+j]);
		        G4 = tmagnitude(hDGImage[ptr_p1+j], vDGImage[ptr_p1+j]);	
		        if((xC*yC) > (float)0.0){
			    G1 = tmagnitude(hDGImage[ptr_m1+j-1], vDGImage[ptr_m1+j-1]);
			    G3 = tmagnitude(hDGImage[ptr_p1+j+1], vDGImage[ptr_p1+j+1]);
		        }
		        else{
			    G1 = tmagnitude(hDGImage[ptr_m1+j+1], vDGImage[ptr_m1+j+1]);
			    G3 = tmagnitude(hDGImage[ptr_p1+j-1], vDGImage[ptr_p1+j-1]);
		        }
		    }
		    else{
		        /* horizontal gradient */
		        xSlope = (float)(fabs(yC) / fabs(xC));
		        ySlope = (float)1.0;
		        G2 = tmagnitude(hDGImage[ptr+j+1], vDGImage[ptr+j+1]);
		        G4 = tmagnitude(hDGImage[ptr+j-1], vDGImage[ptr+j-1]);	
		        if((xC*yC) > (float)0.0){
			    G1 = tmagnitude(hDGImage[ptr_p1+j+1], vDGImage[ptr_p1+j+1]);
			    G3 = tmagnitude(hDGImage[ptr_m1+j-1], vDGImage[ptr_m1+j-1]);
		        }
		        else{
			    G1 = tmagnitude(hDGImage[ptr_m1+j+1], vDGImage[ptr_m1+j+1]);
			    G3 = tmagnitude(hDGImage[ptr_p1+j-1], vDGImage[ptr_p1+j-1]);
		        }
		    }
		    if((G > (xSlope*G1+(ySlope-xSlope)*G2))&&(G > (xSlope*G3+(ySlope-xSlope)*G4))){
		        magImage[ptr+j] = G;	
		    }
		    if(magImage[ptr+j] > maxValue) maxValue = magImage[ptr+j];
		    if(magImage[ptr+j] < minValue) minValue = magImage[ptr+j];
		}
	    }
	}

	scale = (float)1.0 / (maxValue-minValue);
	ptr   = 0;
	count = 0;
	*tAve  = 0.0;
	for(i = 0; i < rows; ++i){
	    for(j = 0; j < cols; ++j){
		magImage[ptr] = scale * (magImage[ptr]-minValue);
		if(magImage[ptr] > 0.0){
		    *tAve += magImage[ptr];
		    ++count;
		}
		++ptr;
	    }
	}
	*tAve /= (float)count;

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
	    *cannyLow  = ((cannyL)  * *tAve);
	    *cannyHigh = ((cannyH) * *tAve);
	}
	else{
	    /* based on the mode value of edge energy */
	    *cannyLow  = ((cannyL)  * ((float)mIndex/step));
	    *cannyHigh = ((cannyH) * ((float)mIndex/step));
	}
	status = 1;

	return status;

}

int trace_Edge(int i, int j, int rows, int cols, double cannyLow, double *magImage,
               unsigned short *hys_image){

	int n, m;
	int ptr;
	int flag;

	ptr = i * cols;
	if(hys_image[ptr+j] == 0){
	    /*
	    // this point is above high threshold
	    */
	    hys_image[ptr+j] = 1;
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
			    if(trace_Edge(i+n, j+m, rows, cols, cannyLow, magImage, hys_image)){
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

int NI_CannyHysteresis(int num, int rows, int cols, double *magImage, unsigned short *hys_image,
		       double cannyLow, double cannyHigh){ 


	int status;
	int i, j;
	int ptr;

	for(i = 0; i < rows; ++i){
	    ptr = i * cols;
	    for(j = 0; j < cols; ++j){
		if(magImage[ptr+j] > cannyHigh){
		    trace_Edge(i, j, rows, cols, cannyLow, magImage, hys_image);
		}
	    }
	}

	status = 1;

	return status;

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

void computeLaws(LawsFilter7 lawsFilter, int aperature, int srcRows, int srcCols, 
                 unsigned short *MaskImage, float *lawsImage, double *sourceImage){

	/*
	// hard-wirred to Law's 7 kernels
	*/
	int i, j;
	int lawsLayer;
	int column, row;
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

	for(i = aperature; i < srcRows-aperature; ++i){
	    // get the row array offset for mask and data source. 
	    for(row = -aperature; row <= aperature; ++row){
		maskOffset[row+aperature] = (i+row)*srcCols;
		dataOffset[row+aperature] = maskOffset[row+aperature];
	    }
	    for(j = aperature; j < srcCols-aperature; ++j){
		/*
		// get 7x7 segment and make sure have 100% mask coverage
		*/
		count = 0;
		for(row = -aperature; row <= aperature; ++row){
		    rowNumber = (row+aperature)*kernelSize;
		    for(column = -aperature; column <= aperature; ++column){
			if(MaskImage[maskOffset[row+aperature]+j+column]){
			    myImage[rowNumber+column+aperature] = sourceImage[dataOffset[row+aperature]+j+column];
			    ++count;
			}
		    }
		}
		if(count == fullMask){
		    /*
		    // 100% mask coverage. now do the Law's texture filters
		    */
		    lawsLayer = 0;
		    for(outerKernelNumber = 0; outerKernelNumber < lawsFilter.numberKernels; ++outerKernelNumber){
			/*
			// outer loop pulls the i'th kernel. kernel 0 is the LP kernel
			// the outer loop is the iso-kernel
			*/
			rowFilter = &lawsFilter.lawsKernel[outerKernelNumber][0];
			colFilter = &lawsFilter.lawsKernel[outerKernelNumber][0];
			filterResult1 = lawsConvolution(myImage, rowFilter, colFilter, kernelSize);
			/* lawsLayer 0 is the LP and needs to be used to scale. */
			if(outerKernelNumber){
			    lawsImage[lawsLayer*layerStep + i*srcCols + j] = (float)2.0 * filterResult1;
			}
			else{
			    lawsLL = filterResult1;
			    lawsLL = (float)2.0 * filterResult1;
			    lawsImage[lawsLayer*layerStep + i*srcCols + j] = (float)2.0 * filterResult1;
			}
			++lawsLayer;
			/*
			// now do the inner loop and get the column filters for the other laws kernels
			*/
			for(innerKernelNumber = outerKernelNumber+1;
			                        innerKernelNumber < lawsFilter.numberKernels;
			                        ++innerKernelNumber){
			    colFilter = &lawsFilter.lawsKernel[innerKernelNumber][0];
			    filterResult1 = lawsConvolution(myImage, rowFilter, colFilter, kernelSize);
			    filterResult2 = lawsConvolution(myImage, colFilter, rowFilter, kernelSize);
			    lawsImage[lawsLayer*layerStep + i*srcCols + j] = filterResult1 + filterResult2;
			    ++lawsLayer;
			}
		    }
		}
	    }
	}

	return;

}


int NI_LawsTexture(int num, int rows, int cols, double *src_image, unsigned short *mask, 
		   float *lawsImage, LawsFilter7 lawsFilter){

	int status;
	int number_kernels;
	int kernel_size;
	int filters;
        number_kernels = lawsFilter.numberKernels;
        kernel_size = lawsFilter.kernelLength;
        filters = lawsFilter.numberFilterLayers;
	int aperature = (kernel_size-1)/2;

	computeLaws(lawsFilter, aperature, rows, cols, mask, lawsImage, src_image);

	status = 1;

	return status;

}


int NI_RoiCoOccurence(int samples, int rows, int cols, unsigned short *labelImage,
	              unsigned short *rawImage, int *cocMatrix, int distance, int orientation){ 

	int i, j;
	int offset;
	int d_row;
	int d_col;
	int status;
	int start_row;
	int stop_row;
	int start_col;
	int stop_col;
	int mask;
	int pixel;
	int d_mask_value;
	int d_pixel_value;

	/* built around 8 bit histograms */

	offset = 0;
	if(orientation == 90){
	    start_row = 0;
	    stop_row  = rows;
	    start_col = 0;
	    stop_col  = cols-distance;
	    d_row     = 0;
	    d_col     = distance;
	}
	else if(orientation == 180){
	    start_row = 0;
	    stop_row  = rows-distance;
	    start_col = 0;
	    stop_col  = cols;
	    d_row     = cols*distance;
	    d_col     = 0;
	}
	else if(orientation == 45){
	    start_row = 0;
	    stop_row  = rows-distance;
	    start_col = distance;
	    stop_col  = cols;
	    d_row     = cols*distance;
	    d_col     = -distance;
	}
	else if(orientation == 135){
	    start_row = 0;
	    stop_row  = rows-distance;
	    start_col = 0;
	    stop_col  = cols-distance;
	    d_row     = cols*distance;
	    d_col     = distance;
	}

	for(i = start_row; i < stop_row; ++i){
	    for(j = start_col; j < stop_col; ++j){
		mask = labelImage[offset+j];
		if(mask){
		    /* d rows away from current row */
		    pixel = rawImage[offset+j];
		    d_mask_value = labelImage[offset+d_row+j+d_col];
		    if(d_mask_value){
		        /* over the mask */
		        d_pixel_value = rawImage[offset+d_row+j+d_col];
			/* update the 2D joint histograms */
	                ++cocMatrix[d_pixel_value*256+pixel];
	                ++cocMatrix[d_pixel_value+pixel*256];
		    }
		}
	    }
	    offset += cols;
	}

	status = 1;

	return(status);

}

int NI_GrowRegion2D(int rows, int cols, double *rawimage, unsigned short *label,
                    objStruct *expanded_ROI, objStruct *newgrow_ROI, double low_threshold,
                    double high_threshold, int Label, int N_connectivity){

	int i, j, p, m;
	int offset;
	int offsetM, offsetP;
	int status;
	int T[8], count;
	int LowX;
	int LowY;
	int HighX;
	int HighY;
	double value;
	bool change;

	while(1){
	    change = FALSE;
	    for(i = 1; i < rows-1; ++i){
	        offset  = i * cols;
	        offsetM = offset - cols;
	        offsetP = offset + cols;
	        for(j = 1; j < cols-1; ++j){
	            m = label[offset+j];
		    if(!m){
			/* un-labeled pixel */
	                value = rawimage[offset+j];
	                if((value > low_threshold) && (value < high_threshold)){
			    /* check for N-connectivity */
			    T[0] = label[offset+j+1];
			    T[1] = label[offsetM+j+1];
			    T[2] = label[offsetM+j];
			    T[3] = label[offsetM+j-1];
			    T[4] = label[offset+j-1];
			    T[5] = label[offsetP+j-1];
			    T[6] = label[offsetP+j];
			    T[7] = label[offsetP+j+1];
			    count = 0;
	        	    for(p = 0; p < 8; ++p){
			        if(T[p] == Label){
			            ++count;
			        }
			    }	
			    if(count > N_connectivity){
	            		label[offset+j] = Label;
	    			change = TRUE;
			    } 
			}	
		    }
	        }
	    }
	    if(!change) break;
	}

	/* get new bounding box */
	newgrow_ROI->Left   = expanded_ROI->Left + LowX;
	newgrow_ROI->Right  = newgrow_ROI->Left + (HighX-LowX);
	newgrow_ROI->Bottom = expanded_ROI->Bottom + LowY;
	newgrow_ROI->Top    = expanded_ROI->Bottom + (HighY-LowY);
	newgrow_ROI->Mass   = count;

	status = 1;

	return(status);

}

int NI_GrowRegion3D(int layers, int rows, int cols, double *rawimage, unsigned short *label,
                    objStruct *expanded_ROI, objStruct *newgrow_ROI, double low_threshold,
		    double high_threshold, int Label, int N_connectivity){

	int i, j, k, m, p;
	int offset;
	int ptr;
	int lOffset,  rOffset;
	int lOffsetP, lOffsetN;
	int rOffsetP, rOffsetN;
	int layerSize;
	int status;
	int T[26], count;
	int LowX;
	int LowY;
	int LowZ;
	int HighX;
	int HighY;
	int HighZ;
	float centerX;
	float centerY;
	float centerZ;
	double value;
	bool change;

	layerSize = rows * cols;
	while(1){
	    change = FALSE;
	    for(i = 1; i < layers-1; ++i){
	        lOffset  = i * layerSize;
	        lOffsetP = lOffset+layerSize;
	        lOffsetN = lOffset-layerSize;
	        for(j = 1; j < rows-1; ++j){
		    rOffset = j * cols;
		    rOffsetP = rOffset+cols;
		    rOffsetN = rOffset-cols;
		    for(k = 1; k < cols-1; ++k){
		        m = label[lOffset+rOffset+k];
		        if(!m){
			    /* un-labeled voxel */
	                    value = rawimage[lOffset+rOffset+k];
	                    if((value > low_threshold) && (value < high_threshold)){
			        /* check for N-connectivity */
			        T[0]  = label[lOffset+rOffset+k+1];
			        T[1]  = label[lOffset+rOffsetN+k+1];
			        T[2]  = label[lOffset+rOffsetN+k];
			        T[3]  = label[lOffset+rOffsetN+k-1];
			        T[4]  = label[lOffset+rOffset+k-1];
			        T[5]  = label[lOffset+rOffsetP+k-1];
			        T[6]  = label[lOffset+rOffsetP+k];
			        T[7]  = label[lOffset+rOffsetP+k+1];

			        T[8]  = label[lOffsetN+rOffset+k];
			        T[9]  = label[lOffsetN+rOffset+k+1];
			        T[10] = label[lOffsetN+rOffsetN+k+1];
			        T[11] = label[lOffsetN+rOffsetN+k];
			        T[12] = label[lOffsetN+rOffsetN+k-1];
			        T[13] = label[lOffsetN+rOffset+k-1];
			        T[14] = label[lOffsetN+rOffsetP+k-1];
			        T[15] = label[lOffsetN+rOffsetP+k];
			        T[16] = label[lOffsetN+rOffsetP+k+1];

			        T[17] = label[lOffsetP+rOffset+k];
			        T[18] = label[lOffsetP+rOffset+k+1];
			        T[19] = label[lOffsetP+rOffsetN+k+1];
			        T[20] = label[lOffsetP+rOffsetN+k];
			        T[21] = label[lOffsetP+rOffsetN+k-1];
			        T[22] = label[lOffsetP+rOffset+k-1];
			        T[23] = label[lOffsetP+rOffsetP+k-1];
			        T[24] = label[lOffsetP+rOffsetP+k];
			        T[25] = label[lOffsetP+rOffsetP+k+1];

			        count = 0;
	        	        for(p = 0; p < 26; ++p){
			            if(T[p] == Label){
			                ++count;
				    }
			        }	
			        if(count > N_connectivity){
		        	    label[lOffset+rOffset+k]= Label;
	    			    change = TRUE;
			        } 
			    } 
		        }
		    }
	        }
	    }
	    if(!change) break;
	}

        LowX       = 32767;
	LowY       = 32767;
	LowZ       = 32767;
	HighX      = 0;
	HighY      = 0;
	HighZ      = 0;
	count      = 0;
	centerX    = (float)0.0;
	centerY    = (float)0.0;
	centerZ    = (float)0.0;
	ptr        = 0;
	count      = 0;
	for(i = 0; i < layers; ++i){
	    for(j = 0; j < rows; ++j){
	        for(k = 0; k < cols; ++k){
		    m = label[ptr++];
		    if(m == Label){
		        if(i < LowZ)   LowZ = i;
		        if(j < LowY)   LowY = j;
		        if(k < LowX)   LowX = k;
		        if(i > HighZ) HighZ = i;
		        if(j > HighY) HighY = j;
		        if(k > HighX) HighX = k;
	    	        centerX += (float)k;
	    	        centerY += (float)j;
	    	        centerZ += (float)i;
	    	        ++count;
		    }
		}
	    }
	}

	newgrow_ROI->Left   = expanded_ROI->Left + LowX;
	newgrow_ROI->Right  = newgrow_ROI->Left + (HighX-LowX);
	newgrow_ROI->Bottom = expanded_ROI->Bottom + LowY;
	newgrow_ROI->Top    = newgrow_ROI->Bottom + (HighY-LowY);
	newgrow_ROI->Front  = expanded_ROI->Front + LowZ;
	newgrow_ROI->Back   = expanded_ROI->Front + (HighZ-LowZ);
	newgrow_ROI->Mass   = count;

	status = 1;

	return(status);

}


