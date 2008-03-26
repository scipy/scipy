#ifndef V1_STRUCTSH
#define V1_STRUCTSH

#define bool unsigned char

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
}objStruct;

#endif
