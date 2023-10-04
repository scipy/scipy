#ifndef STIRLING_H
#define STIRLING_H

/* c implementation of Stirling numbers of the second kind */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

//#include "_extra_special.h"
#include "_lambertw.h"


/*
 *
 *     Stirling numbers of the second kind
 *
 *
 *
 * SYNOPSIS: Stirling numbers of the second kind count the
 *  number of ways to make a partition of n distinct elements
 *  into k non-empty subsets.
 *
 * DESCRIPTION: n is the number of distinct elements of the set
 *  to be partitioned and k is the number of non-empty subsets.
 *  The values for n < 0 or k < 0 are interpreted as 0. If you
 *
 * ACCURACY: The method returns a double type.
 *
 * NOTE: this file is *NOT* part of the originaldistribution and was
 *  added by Lucas Roberts
 */

// Dynamic programming approach

double stirling2_dp(int n, int k){
    if ((n == 0 && k == 0) || (n==1 && k==1)) {
        return 1.;
    }
    if (k <= 0 || k > n || n < 0){
        return 0.;
    }
    int arraySize = k <= n - k + 1 ? k : n - k + 1;
    double *curr = (double *) malloc(arraySize * sizeof(double));
    for (int i = 0; i < arraySize; i++){
        curr[i] = 1.;
    }
    if (k <= n - k + 1) {
        for (int i = 1; i < n - k + 1; i++){
            for (int j = 1; j < k; j++){
                curr[j] = (j + 1) * curr[j] + curr[j - 1];
                // supported in c99: https://devdocs.io/c/numeric/math/isinf
                if (isinf(curr[j])){
                    free(curr);
                    return INFINITY; // numeric overflow
                }
            }
        }
    } else {
        for (int i = 1; i < k; i++){
            for (int j = 1; j < n - k + 1; j++){
                curr[j] = (i + 1) * curr[j - 1] + curr[j];
                // supported in c99: https://devdocs.io/c/numeric/math/isinf
                if (isinf(curr[j])){
                    free(curr);
                    return INFINITY; // numeric overflow
                }
            }
        }
    }
    double output = curr[arraySize - 1];
    free(curr);
    return output;
}


// TODO: replace this standin binomial coefficients
// https://github.com/scipy/scipy/pull/19471
int binom(int n, int k){
  assert(n >= 1);
  assert(k >= 1);
  int64_t res  = 1;
  for (int i = n - k + 1; i <= n; ++i)
      res *= i;
  for (int i = 2; i <= k; ++i)
      res /= i;
  return res;
}

// second order Temme approximation
double stirling2_temme(int n, int k){
  if ((n == k && n >= 0) || (n > 0 && k==1)){
      return 1.;
  }
  if (k <= 0 || k > n || n < 0){
      return 0.;
  }
  double mu = k / (double)n;
  double delta = 1. / mu * pow(exp(-1), (-1/mu));
  // note: lambert returns complex value, we only want the real part
  double x0 = 1. / mu;
  // matching k=0, tolerance=1e-8 from _lambertw.py
  std::complex<double> lwv = scipy::special::lambertw(-delta, 0, 1e-8);
  x0 += lwv.real();
  double t0 = (n - k) / (double)k;
  double F = sqrt(t0/((1 + t0)*(x0 - t0)));
  double A = -n * log(x0) + k * log(pow(exp(1), x0) - 1) - k * t0 + (n - k) * log(t0);
  // write F1 as numerator and denominator and apply Horner rule to num
  double x0t0 = x0*t0;
  double t03 = t0*t0*t0;
  double num = (((2*x0 +1)*x0 + 3)*x0 + ((-6*x0 + 8*x0t0 - 5)*x0t0 -6*t03)*x0t0);
  num += (2*t0*t0 + 4)*t03 - 2*x0*x0*x0;
  double denom = (24*F*(1 + t0) * (1 + t0)*(x0 - t0)*(x0 - t0)*(x0 - t0)*(x0-t0));
  double F1 = num / denom;
  double val = exp(A) * pow(k, n - k) * (F - F1/k);
  val *= binom(n, k);
  return val;
}


/*
 *  This is the main entrypoint from stirling2 which handles dispatch to each
 *  of the (private) functions here in the file that implement specific
 *  ways of approximating the stirling numbers of the second kind
 *
 */

double stirling2_inexact(int n, int k) {
    //stubbed out dispatch, TODO: implement
    return 1.;
}


#endif
