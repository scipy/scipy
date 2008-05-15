#ifndef _SCIPY_DJBFFT_COMMON_H
#define _SCIPY_DJBFFT_COMMON_H

#include <cycliccache.h>

#define	COPYSTD2DJB(SRC,DEST,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+1) = *(SRC+n2); \
  for (j=(N)/2-1,k=2;j>0;--j,k+=2) { \
    *(DEST+k) = *(SRC+n2+j); \
    *(DEST+k+1) = *(SRC+j); \
  } \
}

#define	COPYINVDJB2STD(SRC,DEST,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+n2) = *(SRC+1); \
  for (j=(N)/2-1,k=2;j>0;--j,k+=2) { \
    *(DEST+n2+j) = *(SRC+k); \
    *(DEST+j) = *(SRC+k+1); \
  } \
}

#define	COPYINVDJB2STD2(SRC,DEST,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+(N)-1) = *(SRC+(N)-1); \
  for (j=1,k=1;j<n2;++j,k+=2) { \
    *(DEST+n2+j-1) = *(SRC+k); \
    *(DEST+j) = *(SRC+k+1); \
  } \
}

#define COPYDJB2STD(SRC,DEST,FRQ,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+N-1) = *(SRC+1); \
  for (k=2;k<N-1;k+=2) { \
    j = FRQ[k]; \
    if (j>n2) { \
      j = 2*(N-j); \
      *(DEST+j-1) = *(SRC+k); \
      *(DEST+j) = -*(SRC+k+1); \
    } else { \
      j *= 2; \
      *(DEST+j-1) = *(SRC+k); \
      *(DEST+j) = *(SRC+k+1); \
    } \
  } \
}
#define COPYINVSTD2DJB(SRC,DEST,NORMALIZE,FRQ,N) { \
  int n2 = (N)/2,k,j; \
  if (NORMALIZE) { \
    *(DEST) = *(SRC); \
    *(DEST+1) = *(SRC+N-1); \
  } else { \
    *(DEST) = (*(SRC))*0.5; \
    *(DEST+1) = (*(SRC+N-1))*0.5; \
  } \
  for (k=2;k<N-1;k+=2) { \
    j = FRQ[k]; \
    if (j>n2) { \
      j = 2*(N-j); \
      *(DEST+k) = *(SRC+j-1); \
      *(DEST+k+1) = -*(SRC+j); \
    } else { \
      j *= 2; \
      *(DEST+k) = *(SRC+j-1); \
      *(DEST+k+1) = *(SRC+j); \
    } \
  } \
}
namespace fft {

class DJBFFTCacheId : public CacheId {
        public:
                DJBFFTCacheId(int n) : CacheId(n) {};
};

};

#endif
