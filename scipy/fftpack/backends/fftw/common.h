#ifndef _SCIPY_FFTW_COMMON_H
#define _SCIPY_FFTW_COMMON_H

#define COPYRFFTW2STD(SRC,DEST,N) { \
  int j,n2=(N)/2; \
  *(DEST) = *(SRC); \
  for (j=1;j<n2;++j) { \
    *(DEST+2*j-1) = *(SRC+j); \
    *(DEST+2*j) = *(SRC+(N)-j); \
  } \
  if (N>1) { \
    *(DEST+2*n2-1) = *(SRC+n2); \
    if ((N)%2) \
      *(DEST+2*n2) = *(SRC+(N)-n2); \
  } \
}
#define COPYINVRFFTW2STD(SRC,DEST,N) { \
  int j,n2=(N)/2; \
  *(DEST) = *(SRC); \
  for (j=1;j<n2;++j) { \
    *(DEST+j) = *(SRC+2*j-1); \
    *(DEST+(N)-j) = *(SRC+2*j); \
  } \
  if (N>1) {\
    *(DEST+n2) = *(SRC+2*n2-1); \
    if ((N)%2) \
      *(DEST+(N)-n2) = *(SRC+2*n2); \
  } \
}

#endif
