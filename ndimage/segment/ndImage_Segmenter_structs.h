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
    int Left;
    int Right;
    int Top;
    int Bottom;
    int Front;
    int Back;
    int Label;
    int Mass;
    float cX;
    float cY;
    float cZ;
}objStruct;

#endif
