// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_POWER_METHOD_H_
#define IPX_POWER_METHOD_H_

#include <cmath>
#include "ipx_internal.h"
#include "utils.h"

namespace ipx {

// Power method for estimating the maximum eigenvalue of a linear operator f.
// @func: function object that is called by func(v,fv) to evaluate fv=f(v).
// @v: vector of dimension of the linear operator. On return holds an
//     approximate eigenvector corresponding to the maximum eigenvalue of f.
// Returns an estimate for the maximum eigenvalue of f.

template <typename F>
double PowerMethod(F func, Vector& v) {
    const Int maxiter = 100;
    const double tol = 1e-3;
    const Int dim = v.size();
    Vector fv(dim);

    // Construct starting vector and normalize.
    for (Int i = 0; i < dim; i++)
        v[i] = 1.0 + 1.0/(i+1);
    v /= Twonorm(v);

    // Run power method
    double lambda = 0.0;
    Int iter = 0;
    while (iter++ < maxiter) {
        func(v, fv);
        double lambda_old = lambda;
        lambda = Twonorm(fv);
        v = fv/lambda;
        if (std::abs(lambda-lambda_old) <= tol*lambda)
            break;
    }
    return lambda;
}

}  // namespace ipx

#endif // IPX_POWER_METHOD_H_
