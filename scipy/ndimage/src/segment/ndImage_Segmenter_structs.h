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

#endif
