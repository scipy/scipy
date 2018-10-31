#include <Python.h>
#include <math.h>
#include <numpy/npy_math.h>
#include "misc.h"

double gammasgn(double x)
{
    double fx;

    if (npy_isnan(x)) {
      return x;
    }
    if (x > 0) {
        return 1.0;
    }
    else {
        fx = floor(x);
        if (x - fx == 0.0) {
            return 0.0;
        }
        else if ((int)fx % 2) {
            return -1.0;
        }
        else {
            return 1.0;
        }
    }
}
