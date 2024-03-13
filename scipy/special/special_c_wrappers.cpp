extern "C" {
#include "special_c_wrappers.h"
}

#include "_special.h"

extern "C" npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp) {
    return hyp2f1_complex(a, b, c, zp);
}


    


