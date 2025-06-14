#include "propack/_common.h"

void scompute_mu(float* restrict mu, const int j, const float delta, const float eta, int* restrict indices)
{
    // If exists find peaks higher than delta and check for the skirts of
    // the peaks that are higher than eta on both sides for all peaks.

    if (delta < eta) { return; }  // Malformed input
    int i = 0, ip = 0, k = 0, s = 0, tmp = 0;
    indices[0] = 0;
    while (i < j) {
        // Find next k > i where |mu[k]| > delta
        k = -1;
        for (int kk = i; kk < j; ++kk) {
            if (fabs(mu[kk]) > delta) { k = kk; break; }
        }
        // No peaks, quit while loop
        if (k == -1) { break; }
        // Peak is at k, check right side of the peak for a drop below eta
        s = k;
        for (int ss = k; ss >= i; ss--) {
            if (fabs(mu[ss]) < eta) { s = ss+1; break; }
        }
        // 10
        indices[ip++] = s;
        // Left side of the peak
        i = j;
        for (int kk = k; kk < j; kk++) {
            if (fabs(mu[kk]) < eta) { i = kk; break; }
        }
        // 20
        indices[ip++] = i-1;
    }
    //30
    // Insert a sentinel marking the end
    indices[ip++] = j;
}


void sset_mu(const int k, float* restrict mu, int* const restrict indices, const float val)
{
    int i = 0, p, q;
    while ((indices[i] < k) && (indices[i] > 0))
    {
        p = indices[i];
        q = indices[i+1];
        for (int j = p; j <= q; j++)
        {
            mu[j] = val;
        }
        i += 2;
    }
}


void
supdate_mu(
    float* mumax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha,
    float* restrict beta, const float anorm, const float eps1)
{
    float d = 0.0f;

    if (j == 0)
    {
        d = eps1 * (hypotf(alpha[j], beta[j]) + alpha[0]) + eps1 * anorm;
        mu[0] = eps1 / beta[0];
        *mumax = fabsf(mu[0]);
    } else {
        mu[0] = alpha[0] * nu[0] - alpha[j] * mu[0];
        d = eps1 * (hypotf(alpha[j], beta[j]) + alpha[0]) + eps1 * anorm;
        mu[0] = (mu[0] + copysignf(d, mu[0])) / beta[j];
        *mumax = fabsf(mu[0]);
        for (int k = 1; k < j - 1; k++)
        {
            mu[k] = alpha[k]*nu[k] +beta[k-1]*nu[k-1]-alpha[j]*mu[k];
            d = eps1 * (hypotf(alpha[j], beta[j]) + hypotf(alpha[k], beta[k-1])) + eps1 * anorm;
            mu[k] = (mu[k] + copysignf(d, mu[k])) / beta[j];
            *mumax = fmax(*mumax, fabsf(mu[k]));
        }
        mu[j-1] = beta[j-1] * nu[j-1];
        d = eps1 * (hypotf(alpha[j], beta[j]) + hypotf(alpha[j], beta[j-1])) + eps1 * anorm;
        mu[j-1] = (mu[j-1] + copysignf(d, mu[j-1])) / beta[j];
        *mumax = fmaxf(*mumax, fabsf(mu[j-1]));
    }
    mu[j] = 1.0f;
}


void supdate_nu(
    float* numax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha,
    float* restrict beta, const float anorm, const float eps1)
{
    float d = 0.0f;

    if (j > 0)
    {
        *numax = 0.0f;
        for (int k = 0; k < j-1; k++)
        {
            nu[k] += beta[k] * mu[k+1] + alpha[k] * mu[k] - beta[j-1] * nu[k];
            d = eps1 * (hypotf(alpha[k], beta[k]) + hypotf(alpha[j], beta[j-1])) + eps1 * anorm;
            nu[k] = (nu[k] + copysignf(d, nu[k])) / (alpha[j]);
            *numax = fmaxf(*numax, fabsf(nu[k]));
        }
        nu[j-1] = 1.0f;
    }
}


void dcompute_mu(double* restrict mu, const int j, const double delta, const double eta, int* restrict indices)
{
    // If exists find peaks higher than delta and check for the skirts of
    // the peaks that are higher than eta on both sides for all peaks.

    if (delta < eta) { return; }  // Malformed input
    int i = 0, ip = 0, k = 0, s = 0, tmp = 0;
    indices[0] = 0;
    while (i < j) {
        // Find next k > i where |mu[k]| > delta
        k = -1;
        for (int kk = i; kk < j; ++kk) {
            if (fabs(mu[kk]) > delta) { k = kk; break; }
        }
        // No peaks, quit while loop
        if (k == -1) { break; }
        // Peak is at k, check right side of the peak for a drop below eta
        s = k;
        for (int ss = k; ss >= i; ss--) {
            if (fabs(mu[ss]) < eta) { s = ss+1; break; }
        }
        // 10
        indices[ip++] = s;
        // Left side of the peak
        i = j;
        for (int kk = k; kk < j; kk++) {
            if (fabs(mu[kk]) < eta) { i = kk; break; }
        }
        // 20
        indices[ip++] = i-1;
    }
    //30
    // Insert a sentinel marking the end
    indices[ip++] = j;
}


void dset_mu(const int k, double* restrict mu, int* const restrict indices, const double val)
{
    int i = 0, p, q;
    while ((indices[i] < k) && (indices[i] > 0))
    {
        p = indices[i];
        q = indices[i+1];
        for (int j = p; j <= q; j++)
        {
            mu[j] = val;
        }
        i += 2;
    }
}


void dupdate_mu(
    double* mumax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha,
    double* restrict beta, const double anorm, const double eps1)
{
    double d = 0.0;

    if (j == 0)
    {
        d = eps1 * (hypot(alpha[j], beta[j]) + alpha[0]) + eps1 * anorm;
        mu[0] = eps1 / beta[0];
        *mumax = fabs(mu[0]);
    } else {
        mu[0] = alpha[0] * nu[0] - alpha[j] * mu[0];
        d = eps1 * (hypot(alpha[j], beta[j]) + alpha[0]) + eps1 * anorm;
        mu[0] = (mu[0] + copysign(d, mu[0])) / beta[j];
        *mumax = fabs(mu[0]);
        for (int k = 1; k < j - 1; k++)
        {
            mu[k] = alpha[k]*nu[k] +beta[k-1]*nu[k-1]-alpha[j]*mu[k];
            d = eps1 * (hypot(alpha[j], beta[j]) + hypot(alpha[k], beta[k-1])) + eps1 * anorm;
            mu[k] = (mu[k] + copysign(d, mu[k])) / beta[j];
            *mumax = fmax(*mumax, fabs(mu[k]));
        }
        mu[j-1] = beta[j-1] * nu[j-1];
        d = eps1 * (hypot(alpha[j], beta[j]) + hypot(alpha[j], beta[j-1])) + eps1 * anorm;
        mu[j-1] = (mu[j-1] + copysign(d, mu[j-1])) / beta[j];
        *mumax = fmax(*mumax, fabs(mu[j-1]));
    }
    mu[j] = 1.0;
}


void dupdate_nu(
    double* numax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha,
    double* restrict beta, const double anorm, const double eps1)
{
    double d = 0.0;

    if (j > 0)
    {
        *numax = 0.0;
        for (int k = 0; k < j-1; k++)
        {
            nu[k] += beta[k] * mu[k+1] + alpha[k] * mu[k] - beta[j-1] * nu[k];
            d = eps1 * (hypot(alpha[k], beta[k]) + hypot(alpha[j], beta[j-1])) + eps1 * anorm;
            nu[k] = (nu[k] + copysign(d, nu[k])) / (alpha[j]);
            *numax = fmax(*numax, fabs(nu[k]));
        }
        nu[j-1] = 1.0;
    }
}


