#include "common.h"
#include <math.h>
#include <float.h>


// ==================================================================================
// Random 64bit integer generator for random floating point number purposes; based on
// https://prng.di.unimi.it/xoshiro256plus.c

static uint64_t rol64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t xoshiro256p(uint64_t* s) {
    uint64_t const result = s[0] + s[3];
    uint64_t const t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rol64(s[3], 45);

    return result;
}

double random_double(uint64_t* state) {
    union {
        uint64_t u;
        double d;
    } x;

    uint64_t r = xoshiro256p(state);
    r >>= 11;  // top 53 bits
    x.u = ((uint64_t)1023 << 52) | r;
    return x.d - 1.0;
}

float random_float(uint64_t* state) {
    union {
        uint32_t u;
        float f;
    } x;

    uint32_t r = (uint32_t)(xoshiro256p(state) >> 32);  // top 32 bits, waste half of the bits
    r >>= 9;  // keep top 23 bits
    x.u = (127 << 23) | r;

    return x.f - 1.0f;
}

// ==================================================================================


void sbsvdstep(const int jobu, const int jobv, int m, int n, int k, float sigma, float* D, float* E, float* U, int ldu, float* V, int ldv)
{
    int int1 = 1;
    float c, s, r;
    // Perform an implicit LQ SVD sweep with shift sigma.
    if (k < 2) { return; }

    // Compute the initial rotation based on B*B^T - sigma^2
    float x = D[0]*D[0] - sigma*sigma;
    float y = E[0]*D[0];

    // Chase the bulge down the lower bidiagonal with Givens rotations.
    // Below "y" is the bulge and "x" is the element used to eliminate it.
    for (int i = 0; i < k-1; i++) {
        if (i > 0)
        {
            slartg_(&x, &y, &c, &s, &E[i-1]);
        } else {
            slartg_(&x, &y, &c, &s, &r);
        }
        x = c*D[i] + s*E[i];
        E[i] = -s*D[i] + c*E[i];
        D[i] = x;
        y = s*D[i+1];
        D[i+1] = c*D[i+1];

        if ((jobu) && (m > 0))
        {
            srot_(&m, &U[i*ldu], &int1, &U[(i+1)*ldu], &int1, &c, &s);
        }

        slartg_(&x, &y, &c, &s, &D[i]);
        x = c*E[i] + s*D[i+1];
        D[i+1] = -s*E[i] + c*D[i+1];
        E[i] = x;
        y = s*E[i+1];
        E[i+1] = c*E[i+1];

        if ((jobv) && (n > 0))
        {
            srot_(&n, &V[i*ldv], &int1, &V[(i+1)*ldv], &int1, &c, &s);
        }
    }

    slartg_(&x, &y, &c, &s, &E[k-2]);
    x = c*D[k-1] + s*E[k-1];
    E[k-1] = -s*D[k-1] + c*E[k-1];
    D[k-1] = x;

    if ((jobu) && (m > 0))
    {
        srot_(&m, &U[(k-1)*ldu], &int1, &U[k*ldu], &int1, &c, &s);
    }

    return;
}


void sbdqr(const int ignorelast, const int jobq, const int n, float* restrict D, float* restrict E,
           float* c1, float* c2, float* restrict Qt, int ldq)
{
    float flt1 = 1.0f, flt0 = 0.0f, cs, sn, r;
    if (n < 1) { return; }

    if (jobq)
    {
        // Reset Qt to the identity matrix.
        int nplus1 = n + 1;
        slaset_("A", &nplus1, &nplus1, &flt0, &flt1, Qt, &ldq);
    }
    for (int i = 0; i < n-1; i++)
    {
        slartg_(&D[i], &E[i], &cs, &sn, &r);
        D[i] = r;
        E[i]   = sn*D[i+1];
        D[i+1] = cs*D[i+1];
        if (jobq)
        {
            // Apply the Givens rotation to Qt.
            for (int j = 0; j <= i; j++)
            {
                Qt[i+1 + j*ldq] = -sn*Qt[i + j*ldq];
                Qt[i   + j*ldq] =  cs*Qt[i + j*ldq];
            }
            Qt[i   + (i+1)*ldq] = sn;
            Qt[i+1 + (i+1)*ldq] = cs;
        }
    }
    if (!ignorelast)
    {
        slartg_(&D[n-1], &E[n-1], &cs, &sn, &r);
        D[n-1] = r;
        E[n-1] = 0.0f;
        *c1 = sn;
        *c2 = cs;
        if (jobq)
        {
            // Apply the last Givens rotation to Qt.
            for (int j = 0; j < n; j++)
            {
                Qt[n   + j*ldq] = -sn*Qt[n-1 + j*ldq];
                Qt[n-1 + j*ldq] =  cs*Qt[n-1 + j*ldq];
            }
            Qt[n-1 + n*ldq] = sn;
            Qt[n   + n*ldq] = cs;
        }
    }

    return;
}


void srefinebounds(const int n, const int k, float* restrict theta, float* restrict bound, const float tol, const float eps34)
{
    float gap = 0.0f;
    if (k < 2) { return; };
    for (int i = 0; i < k; i++)
    {
        for (int pm1 = -1; pm1 <= 1; pm1 += 2)
        {
            if (((pm1 == 1) && (i < k-1)) || ((pm1 == -1) && (i > 0)))
            {
                if (fabsf(theta[i] - theta[i+pm1]) < eps34*(theta[i]))
                {
                    if ((bound[i] > tol) && (bound[i + pm1] > tol))
                    {
                        bound[i + pm1] = hypotf(bound[i], bound[i + pm1]);
                        bound[i] = 0.0f;
                    }
                }
            }
        }
    }
    for (int i = 0; i < k; i++)
    {
        if ((i < k-1) || (k == n))
        {
            // We cannot compute a reliable value for the gap of the last
            // Ritz value unless we know it is an approximation to the
            // smallest singular value (k.eq.n). In this case we can take the
            // distance to the next bigger one as the gap, which can really
            // save us from getting stuck on matrices with a single isolated tiny
            // singular value.
            if (i == 0)
            {
                gap = fabsf(theta[i] - theta[i+1]) - fmaxf(bound[i], bound[i+1]);
            } else if (i == n-1) {
                gap = fabsf(theta[i-1] - theta[i]) - fmaxf(bound[i-1], bound[i]);
            } else {
                gap = fabsf(theta[i] - theta[i+1]) - fmaxf(bound[i], bound[i+1]);
                gap = fminf(gap, fabsf(theta[i-1] - theta[i]) - fmaxf(bound[i-1], bound[i]));
            }
            if (gap > bound[i])
            {
                bound[i] = bound[i] * (bound[i] / gap);
            }
        }
    }
}


void scompute_int(float* restrict mu, const int j, const float delta, const float eta, int* restrict indices)
{

    if (delta < eta) { return; } // Malformed input

    int interval_count = 0;
    int current_pos = 0;

    // Process the array from left to right
    while (current_pos <= j) {
        // Find the next peak (value > delta)
        int peak_idx = -1;
        for (int k = current_pos; k <= j; k++) {
            if (fabsf(mu[k]) > delta) {
                peak_idx = k;
                break;
            }
        }

        // If no peak found, we're done
        if (peak_idx == -1) {
            break;
        }

        // Find the left edge of the interval
        // Go backwards from the peak to find where values drop below eta
        int left_edge = peak_idx;
        for (int k = peak_idx; k >= current_pos; k--)
        {
            if (fabsf(mu[k]) >= eta)
            {
                left_edge = k;
            } else {
                break;
            }
        }

        // Find the right edge of the interval
        // Go forwards from the peak to find where values drop below eta
        int right_edge = peak_idx;
        for (int k = peak_idx; k <= j; k++) {
            if (fabsf(mu[k]) >= eta)
            {
                right_edge = k;
            } else {
                break;
            }
        }

        // Step 4: Store the interval [left_edge, right_edge]
        indices[interval_count * 2] = left_edge;
        indices[interval_count * 2 + 1] = right_edge;
        interval_count++;

        // Step 5: Move past this interval for the next search
        current_pos = right_edge + 1;
    }

    // Add a sentinel interval to mark the end of intervals
    indices[interval_count * 2]     = j + 1;
    indices[interval_count * 2 + 1] = j + 1;

}


void sset_mu(const int k, float* restrict mu, int* const restrict indices, const float val)
{
    int i = 0, p, q;
    while (indices[i] <= k)
    {
        if ((i > 1) && (indices[i+1] == 0)) { break; }
        p = indices[i];
        q = indices[i+1];
        for (int j = p; j <= q; j++) { mu[j] = val; }
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
        for (int k = 1; k < j; k++)
        {
            mu[k] = alpha[k]*nu[k] + beta[k-1]*nu[k-1] - alpha[j]*mu[k];
            d = eps1 * (hypotf(alpha[j], beta[j]) + hypotf(alpha[k], beta[k-1])) + eps1 * anorm;
            mu[k] = (mu[k] + copysignf(d, mu[k])) / beta[j];
            *mumax = fmaxf(*mumax, fabsf(mu[k]));
        }
        mu[j] = beta[j-1] * nu[j-1];
        d = eps1 * (hypotf(alpha[j], beta[j]) + hypotf(alpha[j], beta[j-1])) + eps1 * anorm;
        mu[j] = (mu[j] + copysignf(d, mu[j])) / beta[j];
        *mumax = fmaxf(*mumax, fabsf(mu[j]));
    }
    mu[j+1] = 1.0f;
}


void supdate_nu(
    float* numax, float* restrict mu, float* restrict nu, const int j, float* restrict alpha,
    float* restrict beta, const float anorm, const float eps1)
{
    float d = 0.0f;

    if (j > 0)
    {
        *numax = 0.0f;
        for (int k = 0; k < j; k++)
        {
            nu[k] = beta[k] * mu[k+1] + alpha[k] * mu[k] - beta[j-1] * nu[k];
            d = eps1 * (hypotf(alpha[k], beta[k]) + hypotf(alpha[j], beta[j-1])) + eps1 * anorm;
            nu[k] = (nu[k] + copysignf(d, nu[k])) / alpha[j];
            *numax = fmaxf(*numax, fabsf(nu[k]));
        }
        nu[j] = 1.0f;
    }
}


void dbsvdstep(const int jobu, const int jobv, int m, int n, int k, double sigma, double* D, double* E, double* U, int ldu, double* V, int ldv)
{
    int int1 = 1;
    double c, s, r;
    // Perform an implicit LQ SVD sweep with shift sigma.
    if (k < 2) { return; }

    // Compute the initial rotation based on B*B^T - sigma^2
    double x = D[0]*D[0] - sigma*sigma;
    double y = E[0]*D[0];

    // Chase the bulge down the lower bidiagonal with Givens rotations.
    // Below "y" is the bulge and "x" is the element used to eliminate it.
    for (int i = 0; i < k-1; i++) {
        if (i > 0)
        {
            dlartg_(&x, &y, &c, &s, &E[i-1]);
        } else {
            dlartg_(&x, &y, &c, &s, &r);
        }
        x = c*D[i] + s*E[i];
        E[i] = -s*D[i] + c*E[i];
        D[i] = x;
        y = s*D[i+1];
        D[i+1] = c*D[i+1];

        if ((jobu) && (m > 0))
        {
            drot_(&m, &U[i*ldu], &int1, &U[(i+1)*ldu], &int1, &c, &s);
        }

        dlartg_(&x, &y, &c, &s, &D[i]);
        x = c*E[i] + s*D[i+1];
        D[i+1] = -s*E[i] + c*D[i+1];
        E[i] = x;
        y = s*E[i+1];
        E[i+1] = c*E[i+1];

        if ((jobv) && (n > 0))
        {
            drot_(&n, &V[i*ldv], &int1, &V[(i+1)*ldv], &int1, &c, &s);
        }
    }

    dlartg_(&x, &y, &c, &s, &E[k-2]);
    x = c*D[k-1] + s*E[k-1];
    E[k-1] = -s*D[k-1] + c*E[k-1];
    D[k-1] = x;

    if ((jobu) && (m > 0))
    {
        drot_(&m, &U[(k-1)*ldu], &int1, &U[k*ldu], &int1, &c, &s);
    }

    return;
}


void dbdqr(const int ignorelast, const int jobq, const int n, double* restrict D, double* restrict E,
           double* c1, double* c2, double* restrict Qt, int ldq)
{
    double dbl1 = 1.0, dbl0 = 0.0, cs, sn, r;
    if (n < 1) { return; }

    if (jobq)
    {
        // Reset Qt to the identity matrix.
        int nplus1 = n + 1;
        dlaset_("A", &nplus1, &nplus1, &dbl0, &dbl1, Qt, &ldq);
    }
    for (int i = 0; i < n-1; i++)
    {
        dlartg_(&D[i], &E[i], &cs, &sn, &r);
        D[i] = r;
        E[i]   = sn*D[i+1];
        D[i+1] = cs*D[i+1];
        if (jobq)
        {
            // Apply the Givens rotation to Qt.
            for (int j = 0; j <= i; j++)
            {
                Qt[i+1 + j*ldq] = -sn*Qt[i + j*ldq];
                Qt[i   + j*ldq] =  cs*Qt[i + j*ldq];
            }
            Qt[i   + (i+1)*ldq] = sn;
            Qt[i+1 + (i+1)*ldq] = cs;
        }
    }
    if (!ignorelast)
    {
        dlartg_(&D[n-1], &E[n-1], &cs, &sn, &r);
        D[n-1] = r;
        E[n-1] = 0.0;
        *c1 = sn;
        *c2 = cs;
        if (jobq)
        {
            // Apply the last Givens rotation to Qt.
            for (int j = 0; j < n; j++)
            {
                Qt[n   + j*ldq] = -sn*Qt[n-1 + j*ldq];
                Qt[n-1 + j*ldq] =  cs*Qt[n-1 + j*ldq];
            }
            Qt[n-1 + n*ldq] = sn;
            Qt[n   + n*ldq] = cs;
        }
    }

    return;
}


void drefinebounds(const int n, const int k, double* restrict theta, double* restrict bound, const double tol, const double eps34)
{
    double gap = 0.0;
    if (k < 2) { return; };
    for (int i = 0; i < k; i++)
    {
        for (int pm1 = -1; pm1 <= 1; pm1 += 2)
        {
            if (((pm1 == 1) && (i < k-1)) || ((pm1 == -1) && (i > 0)))
            {
                if (fabs(theta[i] - theta[i+pm1]) < eps34*(theta[i]))
                {
                    if ((bound[i] > tol) && (bound[i + pm1] > tol))
                    {
                        bound[i + pm1] = hypot(bound[i], bound[i + pm1]);
                        bound[i] = 0.0;
                    }
                }
            }
        }
    }
    for (int i = 0; i < k; i++)
    {
        if ((i < k-1) || (k == n))
        {
            // We cannot compute a reliable value for the gap of the last
            // Ritz value unless we know it is an approximation to the
            // smallest singular value (k.eq.n). In this case we can take the
            // distance to the next bigger one as the gap, which can really
            // save us from getting stuck on matrices with a single isolated tiny
            // singular value.
            if (i == 0)
            {
                gap = fabs(theta[i] - theta[i+1]) - fmax(bound[i], bound[i+1]);
            } else if (i == n-1) {
                gap = fabs(theta[i-1] - theta[i]) - fmax(bound[i-1], bound[i]);
            } else {
                gap = fabs(theta[i] - theta[i+1]) - fmax(bound[i], bound[i+1]);
                gap = fmin(gap, fabs(theta[i-1] - theta[i]) - fmax(bound[i-1], bound[i]));
            }
            if (gap > bound[i])
            {
                bound[i] = bound[i] * (bound[i] / gap);
            }
        }
    }
}


void dcompute_int(double* restrict mu, const int j, const double delta, const double eta, int* restrict indices)
{

    if (delta < eta) { return; } // Malformed input

    int interval_count = 0;
    int current_pos = 0;

    // Process the array from left to right
    while (current_pos <= j) {
        // Find the next peak (value > delta)
        int peak_idx = -1;
        for (int k = current_pos; k <= j; k++) {
            if (fabs(mu[k]) > delta) {
                peak_idx = k;
                break;
            }
        }

        // If no peak found, we're done
        if (peak_idx == -1) {
            break;
        }

        // Find the left edge of the interval
        // Go backwards from the peak to find where values drop below eta
        int left_edge = peak_idx;
        for (int k = peak_idx; k >= current_pos; k--)
        {
            if (fabs(mu[k]) >= eta)
            {
                left_edge = k;
            } else {
                break;
            }
        }

        // Find the right edge of the interval
        // Go forwards from the peak to find where values drop below eta
        int right_edge = peak_idx;
        for (int k = peak_idx; k <= j; k++) {
            if (fabs(mu[k]) >= eta)
            {
                right_edge = k;
            } else {
                break;
            }
        }

        // Step 4: Store the interval [left_edge, right_edge]
        indices[interval_count * 2] = left_edge;
        indices[interval_count * 2 + 1] = right_edge;
        interval_count++;

        // Step 5: Move past this interval for the next search
        current_pos = right_edge + 1;
    }

    // Add a sentinel interval to mark the end of intervals
    indices[interval_count * 2]     = j + 1;
    indices[interval_count * 2 + 1] = j + 1;

}


void dset_mu(const int k, double* restrict mu, int* const restrict indices, const double val)
{
    int i = 0, p, q;
    while (indices[i] <= k)
    {
        if ((i > 1) && (indices[i+1] == 0)) { break; }
        p = indices[i];
        q = indices[i+1];
        for (int j = p; j <= q; j++) { mu[j] = val; }
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
        for (int k = 1; k < j; k++)
        {
            mu[k] = alpha[k]*nu[k] + beta[k-1]*nu[k-1] - alpha[j]*mu[k];
            d = eps1 * (hypot(alpha[j], beta[j]) + hypot(alpha[k], beta[k-1])) + eps1 * anorm;
            mu[k] = (mu[k] + copysign(d, mu[k])) / beta[j];
            *mumax = fmax(*mumax, fabs(mu[k]));
        }
        mu[j] = beta[j-1] * nu[j-1];
        d = eps1 * (hypot(alpha[j], beta[j]) + hypot(alpha[j], beta[j-1])) + eps1 * anorm;
        mu[j] = (mu[j] + copysign(d, mu[j])) / beta[j];
        *mumax = fmax(*mumax, fabs(mu[j]));
    }
    mu[j+1] = 1.0;
}


void dupdate_nu(
    double* numax, double* restrict mu, double* restrict nu, const int j, double* restrict alpha,
    double* restrict beta, const double anorm, const double eps1)
{
    double d = 0.0;

    if (j > 0)
    {
        *numax = 0.0;
        for (int k = 0; k < j; k++)
        {
            nu[k] = beta[k] * mu[k+1] + alpha[k] * mu[k] - beta[j-1] * nu[k];
            d = eps1 * (hypot(alpha[k], beta[k]) + hypot(alpha[j], beta[j-1])) + eps1 * anorm;
            nu[k] = (nu[k] + copysign(d, nu[k])) / alpha[j];
            *numax = fmax(*numax, fabs(nu[k]));
        }
        nu[j] = 1.0;
    }
}
