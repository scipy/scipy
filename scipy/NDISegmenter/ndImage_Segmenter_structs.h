#ifndef V1_STRUCTSH
#define V1_STRUCTSH

#define bool unsigned char

typedef struct{
    int x;
    int y;
}POINT;

typedef struct{
    int x;
    int y;
    int linkIndex;
    bool haveLink;
}bPOINT;

typedef struct{
    int left;
    int right;
    int top;
    int bottom;
}RECT;

typedef struct{
    char filterName[20];
    float Mean;
    float Variance;
}tTEM;

typedef struct{
    int numberKernels;
    int kernelLength;
    int numberFilterLayers;
    float lawsKernel[6][7];
    char name[7];
}LawsFilter7;

typedef struct{
    // filled in GetObjectStats 
    int L;
    int R;
    int T;
    int B;
    int Label;
    int Area;
    float cX;
    float cY;
    // filled in BuildBoundary
    int   curveClose;
    float cXBoundary;
    float cYBoundary;
    float boundaryLength;
    float minRadius;
    float maxRadius;
    float aveRadius;
    float ratio;
    float compactness;
    // filled in VoxelMeasures
    float voxelMean;
    float voxelVar;
    // filled in TextureMeasures
    float TEM[20];
}objStruct;

typedef struct{
    int numberPoints;
    int curveClose;
    int classify;
    float boundaryLength;
    float minRadius;
    float maxRadius;
    float aveRadius;
    float ratio;
    float compactness;
    float voxelMean;
    float voxelVar;
    RECT rectangle;
    POINT centroid;
    bool isWithin;
    bool closedCurve;
    bool criticalSize;
    int Label;
}boundaryIndex;


typedef struct{
    POINT xy;
}blobBoundary;


//
// prototypes
//
int NI_RegionGrow(int, int, int, int, int, int, int, double *, unsigned short *, int *);   
int NI_TextureMeasures(int, int, int, int, double *, unsigned short *, objStruct objectMetrics[]);
int NI_VoxelMeasures(int, int, int, int, double *, unsigned short *, objStruct objectMetrics[]);
int NI_BuildBoundary(int, int, int, int, unsigned short *, objStruct objectMetrics[]);
int NI_GetObjectStats(int, int, int, unsigned short *, objStruct objectMetrics[]);
int NI_ThinFilter(int, int, int, int, unsigned short *, objStruct objectMetrics[]);
int NI_SobelEdges(int, int, int, double, int, int, int, double, int, double *, unsigned short *, int *);  
int NI_ShenCastanEdges(int, int, int, double, double, int, int, int, double *, unsigned short *, int *);
int NI_CannyEdges(int, int, int, double, double, double, int, int, int, double, int,
	          double *, unsigned short *, int *);

void computeLaws(LawsFilter7, tTEM LawsFeatures[], RECT, int, int, int, int, unsigned char *, float *,
	       	 unsigned short *, float *, double *);
float lawsConvolution(float *, float *, float *, int);
void initLaws(LawsFilter7*);
void getVoxelMeasures(objStruct objectMetrics[], double *, unsigned short *, int, int, int);
void getLawsTexture(LawsFilter7, tTEM LawsFeatures[], objStruct objectMetrics[], double *, unsigned short *, int, int, int);
		      
void morphoFilterBinaryImage(int, int, unsigned short *, int, int);
void buildBinaryImage(int, int, double *, unsigned short *, int, int);
void doRegionGrow(int, int, int, double *, unsigned short *, int, int, int, int);
void buildBoundary(objStruct objectMetrics[], int, unsigned short *, int, int, int);
void getBoundary(unsigned short *, unsigned char *, blobBoundary *, blobBoundary *, 
	         boundaryIndex *, RECT, int, int, int, int, int, int);
void doMorphology(unsigned char *, unsigned char *, unsigned char *, unsigned char *, int olapValuesC[],
       	          int olapValuesO[], unsigned short cmask[11][11], unsigned short omask[11][11],
	          RECT, int, int, int, int);
void getCompactness(unsigned char *, RECT, int, int, float *, float);
void OpenCloseFilter(int olapValues[], int, int, int, int, unsigned char *,  
                     unsigned char *, unsigned short mask[11][11]);
void trackBoundary(unsigned char *, blobBoundary lBoundary[], int, int, blobBoundary, int); 
void getBoundaryMetrics(bPOINT *, float *, float *, float *, float *, float, float, int);
void generateMask(unsigned char *, bPOINT *, int, int, int);
void ThinningFilter(int, int, int, int J_mask[3][30], int K_mask[3][30], unsigned char *, 
	            unsigned char *, unsigned char *, unsigned char *, unsigned char *, unsigned char *);
void initThinFilter(int J_mask[3][30], int K_mask[3][30]);
void Shen_Castan(double, double, int, int, int, int, int, double *, unsigned short *);
void computeISEF(float *, float *, int, int, double);
void ISEF_Horizontal(float *, float *, float *, float *, int, int, double);
void ISEF_Vertical(float *, float *, float *, float *, int, int, double);
void thresholdImage(float *, float *, int, int, int, int);
void computeBandedLaplacian(float *, float *, float *, int, int);
void getZeroCrossings(float *, float *, float *, int, int, int);
float adaptiveGradient(float *, float *, int, int, int, int);
void thresholdEdges(float *, unsigned short *, double, int, int);
void estimateThreshold(float *, float *, float, int, int, float *);
void doSobel(int, int, int, double, int, double *, unsigned short *);
void DGFilters(int, int, int, double, int, float *, float *, double *, double *, float *, float *);
void nonMaxSupress(int, int, float, float, double *, double *, int, float *, float *, float *);
void edgeHysteresis(int, int, double, double, float *, float *);
void edgeThreshold(int, int, double, float *, float *);
int traceEdge(int, int, int, int, double, float *, float *);
float magnitude(float, float);
int ConnectedEdgePoints(int, int, unsigned short *);
void doPreProcess(int, int, int, double *, double, int, int, int);
void filter2D(int, int, int, int, int, float *, double *);
void buildKernel(double, int, int, float *);



#endif
