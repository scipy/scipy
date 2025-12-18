#include "_common_array_utils.h"

static float snorm1(float*, float*, const Py_ssize_t);
static double dnorm1(double*, double*, const Py_ssize_t);
static float cnorm1(SCIPY_C*, float*, const Py_ssize_t);
static double znorm1(SCIPY_Z*, double*, const Py_ssize_t);
static float snorm1est(float*, int);
static double dnorm1est(double*, int);
static float cnorm1est(SCIPY_C*, int);
static double znorm1est(SCIPY_Z*, int);


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************/

 static float
snorm1(float* A, float* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    float temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = fabsf(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += fabsf(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

static double
dnorm1(double* A, double* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    double temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = fabs(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += fabs(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

static float
cnorm1(SCIPY_C* A, float* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    float temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = cabsf(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += cabsf(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

static double
znorm1(SCIPY_Z* A, double* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    double temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = cabs(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += cabs(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

static float
snorm1est(float* A, int n)
{
    int from = 2*n, to = n, kase = 0, tempint, int1 = 1;
    int isave[3];
    float est, dbl1 = 1.0, dbl0 = 0.0;
    char* opA;
    float* work_arr = PyMem_RawMalloc(3*n*sizeof(float));
    if (!work_arr) { return -100; }
    int* iwork_arr = PyMem_RawMalloc(n*sizeof(int));
    if (!iwork_arr) { PyMem_RawFree(work_arr);return -101; }
    // 1-norm estimator by reverse communication
    // dlacon( n, v, x, isgn, est, kase )
    slacn2_(&n, work_arr, &work_arr[n], iwork_arr, &est, &kase, isave);

    while (kase)
    {
        opA = (kase == 1 ? "T" : "N");
        tempint = from;
        from = to;
        to = tempint;
        sgemv_(opA, &n, &n, &dbl1, A, &n, &work_arr[from], &int1, &dbl0, &work_arr[to], &int1);
        slacn2_(&n, work_arr, &work_arr[to], iwork_arr, &est, &kase, isave);
    }

    PyMem_RawFree(work_arr);
    PyMem_RawFree(iwork_arr);

    return est;
}


static double
dnorm1est(double* A, int n)
{
    int from = 2*n, to = n, kase = 0, tempint, int1 = 1;
    int isave[3];
    double est, dbl1 = 1.0, dbl0 = 0.0;
    char* opA;
    double* work_arr = PyMem_RawMalloc(3*n*sizeof(double));
    if (!work_arr) { return -100; }
    int* iwork_arr = PyMem_RawMalloc(n*sizeof(int));
    if (!iwork_arr) { PyMem_RawFree(work_arr);return -101; }
    // 1-norm estimator by reverse communication
    // dlacon( n, v, x, isgn, est, kase )
    dlacn2_(&n, work_arr, &work_arr[n], iwork_arr, &est, &kase, isave);

    while (kase)
    {
        opA = (kase == 1 ? "T" : "N");
        tempint = from;
        from = to;
        to = tempint;
        dgemv_(opA, &n, &n, &dbl1, A, &n, &work_arr[from], &int1, &dbl0, &work_arr[to], &int1);
        dlacn2_(&n, work_arr, &work_arr[to], iwork_arr, &est, &kase, isave);
    }

    PyMem_RawFree(work_arr);
    PyMem_RawFree(iwork_arr);

    return est;
}


static float
cnorm1est(SCIPY_C* A, int n)
{
    int from = 2*n, to = n, kase = 0, tempint, int1 = 1;
    int isave[3];
    float est;
    SCIPY_C dbl1 = CPLX_C(1.0, 0.0), dbl0 = CPLX_C(0.0, 0.0);
    char* opA;
    SCIPY_C* work_arr = PyMem_RawMalloc(3*n*sizeof(SCIPY_C));
    if (!work_arr) { return -100; }
    // clacon( n, v, x, est, kase )
    clacn2_(&n, work_arr, &work_arr[n], &est, &kase, isave);

    while (kase)
    {
        opA = (kase == 1 ? "C" : "N");
        tempint = from;
        from = to;
        to = tempint;
        cgemv_(opA, &n, &n, &dbl1, A, &n, &work_arr[from], &int1, &dbl0, &work_arr[to], &int1);
        clacn2_(&n, work_arr, &work_arr[to], &est, &kase, isave);
    }

    PyMem_RawFree(work_arr);

    return est;
}


static double
znorm1est(SCIPY_Z* A, int n)
{
    int from = 2*n, to = n, kase = 0, tempint, int1 = 1;
    int isave[3];
    double est;
    SCIPY_Z dbl1 = CPLX_Z(1.0, 0.0), dbl0 = CPLX_Z(0.0, 0.0);
    char* opA;
    SCIPY_Z* work_arr = PyMem_RawMalloc(3*n*sizeof(SCIPY_Z));
    if (!work_arr) { return -100; }
    // zlacon( n, v, x, est, kase )
    zlacn2_(&n, work_arr, &work_arr[n], &est, &kase, isave);

    while (kase)
    {
        opA = (kase == 1 ? "C" : "N");
        tempint = from;
        from = to;
        to = tempint;
        zgemv_(opA, &n, &n, &dbl1, A, &n, &work_arr[from], &int1, &dbl0, &work_arr[to], &int1);
        zlacn2_(&n, work_arr, &work_arr[to], &est, &kase, isave);
    }

    PyMem_RawFree(work_arr);

    return est;
}

/*******************************************************************************
 *******************************************************************************
 *******************************************************************************/


static void
pick_pade_structure_s(float* Am, const Py_ssize_t size_n, int* m, int* s)
{
    Py_ssize_t i, j;
    Py_ssize_t dims[2];
    int lm = 0, int1 = 1, n = (int)size_n;
    float normA, dbl1 = 1.0, dbl0 = 0.0;
    float d4, d6, d8, d10, eta0, eta1, eta2, eta3, eta4, two_pow_s, temp, test;
    float theta[5];
    float coeff[5];
    // work_arr is two n cols that will be used multiplying absA, alternating.
    float* work_arr = PyMem_RawMalloc(2*n*sizeof(float));
    if (!work_arr) { *m = -1; return; }
    float* absA = PyMem_RawMalloc(n*n*sizeof(float));
    if (!absA) { *m = -2; PyMem_RawFree(work_arr); return; }

    dims[0] = n;
    dims[1] = n;
    theta[0] = 1.495585217958292e-002;
    theta[1] = 2.539398330063230e-001;
    theta[2] = 9.504178996162932e-001;
    theta[3] = 2.097847961257068e+000;
    theta[4] = 4.250000000000000e+000;
    coeff[0] = 6.0081481933593750e-03;
    coeff[1] = 5.9956512451171875e+02;
    coeff[2] = 2.6750196800000000e+08;
    coeff[3] = 3.5252480874905600e+14;
    coeff[4] = 6.7502722515508048e+27;

    // Initialize the first n part of work_arr
    for (i = 0; i < n; i++) { work_arr[i] = 1.0; }

    // absA = np.abs(Am[0, :n, :n])
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            absA[j + i*n] = fabsf(Am[j + i*n]);
        }
    }

    // First spin = normest(|A|, m=1)
    sgemv_("N", &n, &n, &dbl1, absA, &n, work_arr, &int1, &dbl0, &work_arr[n], &int1);
    normA = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > normA) { normA = work_arr[n+i]; } }


    sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[0*n*n], &n, &Am[0*n*n], &n, &dbl0, &Am[1*n*n], &n);
    sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[1*n*n], &n, &Am[1*n*n], &n, &dbl0, &Am[2*n*n], &n);
    sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n*n], &n, &Am[1*n*n], &n, &dbl0, &Am[3*n*n], &n);
    d4 = powf(snorm1(&Am[2*n*n], work_arr, n), 0.25);
    d6 = powf(snorm1(&Am[3*n*n], work_arr, n), 1.0/6.0);
    eta0 = fmaxf(d4, d6);
    eta1 = eta0;

    // m = 3
    // -------
    // 1-norm of A**7
    // Alternating matvecs
    // absA * work_arr[n:] = work_arr[:n]
    // absA * work_arr[:n] = work_arr[n:]
    for (i = 0; i < 3; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[0])/6);
    lm = (lm < 0 ? 0 : lm);
    if ((eta0 < theta[0]) && lm == 0)
    {
        *m = 3;
        goto FREE_EXIT;
    }

    // m = 5
    // -------
    // 1-norm of A**11
    for (i = 0; i < 2; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[1])/10);
    lm = (lm < 0 ? 0 : lm);

    if ((eta1 < theta[1]) && lm == 0)
    {
        *m = 5;
        goto FREE_EXIT;
    }

    // m = 7
    // -------
    if (n < 400)
    {
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &dbl0, &Am[4*n*n], &n);
        d8 = powf(snorm1(&Am[4*n*n], work_arr, n), 0.125);
    } else {
        test = snorm1est(&Am[0], 8);
        // If memory error in s1normest
        if (test <= -100.0) { *m = -3; goto FREE_EXIT; }
        d8 = powf(test, 0.125);
    }

    eta2 = fmaxf(d6, d8);

    // 1-norm of A**15
    for (i = 0; i < 2; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[2])/14);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[2]) && lm == 0)
    {
        *m = 7;
        goto FREE_EXIT;
    }

    // m = 9
    // -------
    // 1-norm of A**19
    for (i = 0; i < 2; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[3])/18);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[3]) && lm == 0)
    {
        if (n >= 400) {
            sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &dbl0, &Am[4*n*n], &n);
        }
        *m = 9;
        goto FREE_EXIT;
    }

    // m = 13
    // -------
    // Scale-square
    if (n < 400)
    {
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[3*n*n], &n, &Am[2*n*n], &n, &dbl0, &Am[4*n*n], &n);
        d10 = powf(snorm1(&Am[4*n*n], work_arr, n), 0.1);
    } else {
        test = snorm1est(&Am[0], 10);
        // If memory error in s1normest
        if (test <= -100.0) { *m = -4; goto FREE_EXIT; }
        d10 = powf(test, 0.1);
    }

    eta3 = fmaxf(d8, d10);
    eta4 = fminf(eta2, eta3);
    *s = (int)ceilf(log2f(eta4 / theta[4]));
    *s = (*s < 0 ? 0 : *s);

    if (*s != 0)
    {
        two_pow_s = powf(2.0, -(*s));
        for (i = 0; i < n*n; i++) { absA[i] *= two_pow_s; }

        // work_arr has spun 19 times already
        two_pow_s = powf(2.0, -19.0*(*s));
        for (i = 0; i < n; i++) { work_arr[n + i] *= two_pow_s; }
        normA *= powf(2.0, -(*s));
    }

    // 1-norm of A**27
    for (i = 0; i < 4; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[4])/26.0);
    *s += (lm < 0 ? 0 : lm );
    *m = 13;

    if (*s != 0)
    {
        float s0 = powf(2.0, -*s);
        float s1 = powf(4.0, -*s);
        float s2 = powf(16.0, -*s);
        float s3 = powf(64.0, -*s);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[        n*i + j] *= s0;
                Am[  n*n + n*i + j] *= s1;
                Am[2*n*n + n*i + j] *= s2;
                Am[3*n*n + n*i + j] *= s3;
            }
        }
    }

FREE_EXIT:
    PyMem_RawFree(absA);
    PyMem_RawFree(work_arr);
    return;
}


static void
pick_pade_structure_d(double* Am, const Py_ssize_t size_n, int* m, int* s)
{
    Py_ssize_t i, j;
    Py_ssize_t dims[2];
    int lm = 0, int1 = 1, n = (int)size_n;
    double normA, dbl1 = 1.0, dbl0 = 0.0;
    double d4, d6, d8, d10, eta0, eta1, eta2, eta3, eta4, two_pow_s, temp, test;
    double theta[5];
    double coeff[5];
    // work_arr is two n cols that will be used multiplying absA, alternating.
    double* work_arr = PyMem_RawMalloc(2*n*sizeof(double));
    if (!work_arr) { *m = -1; return; }
    double* absA = PyMem_RawMalloc(n*n*sizeof(double));
    if (!absA) { *m = -2; PyMem_RawFree(work_arr); return; }

    dims[0] = n;
    dims[1] = n;
    theta[0] = 1.495585217958292e-002;
    theta[1] = 2.539398330063230e-001;
    theta[2] = 9.504178996162932e-001;
    theta[3] = 2.097847961257068e+000;
    theta[4] = 4.250000000000000e+000;
    coeff[0] = 1.1191048088221578e-11;
    coeff[1] = 1.1167770708198077e-06;
    coeff[2] = 4.9826124310493469e-01;
    coeff[3] = 6.5662862500000000e+05;
    coeff[4] = 1.2573361865339437e+19;

    // Initialize the first n part of work_arr
    for (i = 0; i < n; i++) { work_arr[i] = 1.0; }

    // absA = np.abs(Am[0, :n, :n])
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            absA[j + i*n] = fabs(Am[j + i*n]);
        }
    }

    // First spin = normest(|A|, m=1)
    dgemv_("N", &n, &n, &dbl1, absA, &n, work_arr, &int1, &dbl0, &work_arr[n], &int1);
    normA = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > normA) { normA = work_arr[n+i]; } }


    dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[0*n*n], &n, &Am[0*n*n], &n, &dbl0, &Am[1*n*n], &n);
    dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[1*n*n], &n, &Am[1*n*n], &n, &dbl0, &Am[2*n*n], &n);
    dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n*n], &n, &Am[1*n*n], &n, &dbl0, &Am[3*n*n], &n);
    d4 = pow(dnorm1(&Am[2*n*n], work_arr, n), 0.25);
    d6 = pow(dnorm1(&Am[3*n*n], work_arr, n), 1.0/6.0);
    eta0 = fmax(d4, d6);
    eta1 = eta0;

    // m = 3
    // -------
    // 1-norm of A**7
    // Alternating matvecs
    // absA * work_arr[n:] = work_arr[:n]
    // absA * work_arr[:n] = work_arr[n:]
    for (i = 0; i < 3; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[0])/6);
    lm = (lm < 0 ? 0 : lm);
    if ((eta0 < theta[0]) && lm == 0)
    {
        *m = 3;
        goto FREE_EXIT;
    }

    // m = 5
    // -------
    // 1-norm of A**11
    for (i = 0; i < 2; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[1])/10);
    lm = (lm < 0 ? 0 : lm);

    if ((eta1 < theta[1]) && lm == 0)
    {
        *m = 5;
        goto FREE_EXIT;
    }

    // m = 7
    // -------
    if (n < 400)
    {
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &dbl0, &Am[4*n*n], &n);
        d8 = pow(dnorm1(&Am[4*n*n], work_arr, n), 0.125);
    } else {
        test = dnorm1est(&Am[0], 8);
        // If memory error in d1normest
        if (test <= -100.0) { *m = -3; goto FREE_EXIT; }
        d8 = pow(test, 0.125);
    }

    eta2 = fmax(d6, d8);

    // 1-norm of A**15
    for (i = 0; i < 2; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[2])/14);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[2]) && lm == 0)
    {
        *m = 7;
        goto FREE_EXIT;
    }

    // m = 9
    // -------
    // 1-norm of A**19
    for (i = 0; i < 2; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[3])/18);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[3]) && lm == 0)
    {
        if (n >= 400) {
            dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &dbl0, &Am[4*n*n], &n);
        }
        *m = 9;
        goto FREE_EXIT;
    }

    // m = 13
    // -------
    // Scale-square
    if (n < 400)
    {
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[3*n*n], &n, &Am[2*n*n], &n, &dbl0, &Am[4*n*n], &n);
        d10 = pow(dnorm1(&Am[4*n*n], work_arr, n), 0.1);
    } else {
        test = dnorm1est(&Am[0], 10);
        // If memory error in d1normest
        if (test <= -100.0) { *m = -4; goto FREE_EXIT; }
        d10 = pow(test, 0.1);
    }

    eta3 = fmax(d8, d10);
    eta4 = fmin(eta2, eta3);
    *s = (int)ceil(log2(eta4 / theta[4]));
    *s = (*s < 0 ? 0 : *s);

    if (*s != 0)
    {
        two_pow_s = pow(2.0, -(*s));
        for (i = 0; i < n*n; i++) { absA[i] *= two_pow_s; }

        // work_arr has spun 19 times already
        two_pow_s = pow(2.0, -19.0*(*s));
        for (i = 0; i < n; i++) { work_arr[n + i] *= two_pow_s; }
        normA *= pow(2.0, -(*s));
    }

    // 1-norm of A**27
    for (i = 0; i < 4; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[4])/26.0);
    *s += (lm < 0 ? 0 : lm );
    *m = 13;

    if (*s != 0)
    {
        double s0 = pow(2.0, -*s);
        double s1 = pow(4.0, -*s);
        double s2 = pow(16.0, -*s);
        double s3 = pow(64.0, -*s);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[        n*i + j] *= s0;
                Am[  n*n + n*i + j] *= s1;
                Am[2*n*n + n*i + j] *= s2;
                Am[3*n*n + n*i + j] *= s3;
            }
        }
    }

FREE_EXIT:
    PyMem_RawFree(absA);
    PyMem_RawFree(work_arr);
    return;
}


static void
pick_pade_structure_c(SCIPY_C* Am, const Py_ssize_t size_n, int* m, int* s)
{
    Py_ssize_t i, j;
    Py_ssize_t dims[2];
    int lm = 0, int1 = 1, n = (int)size_n;
    float normA;
    float dbl1 = 1.0, dbl0 = 0.0;
    SCIPY_C cdbl1 = CPLX_C(1.0, 0.0), cdbl0 = CPLX_C(0.0, 0.0);
    float d4, d6, d8, d10, eta0, eta1, eta2, eta3, eta4, two_pow_s, temp, test;
    float theta[5];
    float coeff[5];
    // work_arr is two n cols that will be used multiplying absA, alternating.
    float* work_arr = PyMem_RawMalloc(2*n*sizeof(float));
    if (!work_arr) { *m = -1; return; }
    float* absA = PyMem_RawMalloc(n*n*sizeof(float));
    if (!absA) { *m = -2; PyMem_RawFree(work_arr); return; }

    dims[0] = n;
    dims[1] = n;
    theta[0] = 1.495585217958292e-002;
    theta[1] = 2.539398330063230e-001;
    theta[2] = 9.504178996162932e-001;
    theta[3] = 2.097847961257068e+000;
    theta[4] = 4.250000000000000e+000;
    coeff[0] = 6.0081481933593750e-03;
    coeff[1] = 5.9956512451171875e+02;
    coeff[2] = 2.6750196800000000e+08;
    coeff[3] = 3.5252480874905600e+14;
    coeff[4] = 6.7502722515508048e+27;

    // Initialize the first n part of work_arr
    for (i = 0; i < n; i++) { work_arr[i] = 1.0; }

    // absA = np.abs(Am[0, :n, :n])
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            absA[j + i*n] = cabsf(Am[j + i*n]);
        }
    }

    // First spin = normest(|A|, m=1)
    sgemv_("N", &n, &n, &dbl1, absA, &n, work_arr, &int1, &dbl0, &work_arr[n], &int1);
    normA = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > normA) { normA = work_arr[n+i]; } }


    cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[0*n*n], &n, &Am[0*n*n], &n, &cdbl0, &Am[1*n*n], &n);
    cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[1*n*n], &n, &Am[1*n*n], &n, &cdbl0, &Am[2*n*n], &n);
    cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n*n], &n, &Am[1*n*n], &n, &cdbl0, &Am[3*n*n], &n);
    d4 = powf(cnorm1(&Am[2*n*n], work_arr, n), 0.25);
    d6 = powf(cnorm1(&Am[3*n*n], work_arr, n), 1.0/6.0);
    eta0 = fmaxf(d4, d6);
    eta1 = eta0;

    // m = 3
    // -------
    // 1-norm of A**7
    // Alternating matvecs
    // absA * work_arr[n:] = work_arr[:n]
    // absA * work_arr[:n] = work_arr[n:]
    for (i = 0; i < 3; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[0])/6);
    lm = (lm < 0 ? 0 : lm);
    if ((eta0 < theta[0]) && lm == 0)
    {
        *m = 3;
        goto FREE_EXIT;
    }

    // m = 5
    // -------
    // 1-norm of A**11
    for (i = 0; i < 2; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[1])/10);
    lm = (lm < 0 ? 0 : lm);

    if ((eta1 < theta[1]) && lm == 0)
    {
        *m = 5;
        goto FREE_EXIT;
    }

    // m = 7
    // -------
    if (n < 400)
    {
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &cdbl0, &Am[4*n*n], &n);
        d8 = powf(cnorm1(&Am[4*n*n], work_arr, n), 0.125);
    } else {
        test = cnorm1est(&Am[0], 8);
        // If memory error in c1normest
        if (test <= -100.0) { *m = -3; goto FREE_EXIT; }
        d8 = powf(test, 0.125);
    }

    eta2 = fmaxf(d6, d8);

    // 1-norm of A**15
    for (i = 0; i < 2; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[2])/14);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[2]) && lm == 0)
    {
        *m = 7;
        goto FREE_EXIT;
    }

    // m = 9
    // -------
    // 1-norm of A**19
    for (i = 0; i < 2; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[3])/18);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[3]) && lm == 0)
    {
        if (n >= 400) {
            cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &cdbl0, &Am[4*n*n], &n);
        }
        *m = 9;
        goto FREE_EXIT;
    }

    // m = 13
    // -------
    // Scale-square
    if (n < 400)
    {
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[3*n*n], &n, &Am[2*n*n], &n, &cdbl0, &Am[4*n*n], &n);
        d10 = powf(cnorm1(&Am[4*n*n], work_arr, n), 0.1);
    } else {
        test = cnorm1est(&Am[0], 10);
        // If memory error in c1normest
        if (test <= -100.0) { *m = -4; goto FREE_EXIT; }
        d10 = powf(test, 0.1);
    }

    eta3 = fmaxf(d8, d10);
    eta4 = fminf(eta2, eta3);
    *s = (int)ceilf(log2f(eta4 / theta[4]));
    *s = (*s < 0 ? 0 : *s);

    if (*s != 0)
    {
        two_pow_s = powf(2.0, -(*s));
        for (i = 0; i < n*n; i++) { absA[i] *= two_pow_s; }

        // work_arr has spun 19 times already
        two_pow_s = powf(2.0, -19.0*(*s));
        for (i = 0; i < n; i++) { work_arr[n + i] *= two_pow_s; }
        normA *= powf(2.0, -(*s));
    }

    // 1-norm of A**27
    for (i = 0; i < 4; i++)
    {
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        sgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceilf(log2f(temp/normA/coeff[4])/26.0);
    *s += (lm < 0 ? 0 : lm );
    *m = 13;

    if (*s != 0)
    {
        float s0 = powf(2.0, -*s);
        float s1 = powf(4.0, -*s);
        float s2 = powf(16.0, -*s);
        float s3 = powf(64.0, -*s);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                Am[        n*i + j] = _FCmulcr(Am[        n*i + j], s0);
                Am[  n*n + n*i + j] = _FCmulcr(Am[  n*n + n*i + j], s1);
                Am[2*n*n + n*i + j] = _FCmulcr(Am[2*n*n + n*i + j], s2);
                Am[3*n*n + n*i + j] = _FCmulcr(Am[3*n*n + n*i + j], s3);
#else
                Am[        n*i + j] *= s0;
                Am[  n*n + n*i + j] *= s1;
                Am[2*n*n + n*i + j] *= s2;
                Am[3*n*n + n*i + j] *= s3;
#endif
            }
        }
    }

FREE_EXIT:
    PyMem_RawFree(absA);
    PyMem_RawFree(work_arr);
    return;
}


static void
pick_pade_structure_z(SCIPY_Z* Am, const Py_ssize_t size_n, int* m, int* s)
{
    Py_ssize_t i, j;
    Py_ssize_t dims[2];
    int lm = 0, int1 = 1, n = (int)size_n;
    double normA;
    double dbl1 = 1.0, dbl0 = 0.0;
    SCIPY_Z cdbl1 = CPLX_Z(1.0, 0.0), cdbl0 = CPLX_Z(0.0, 0.0);
    double d4, d6, d8, d10, eta0, eta1, eta2, eta3, eta4, two_pow_s, temp, test;
    double theta[5];
    double coeff[5];
    // work_arr is two n cols that will be used multiplying absA, alternating.
    double* work_arr = PyMem_RawMalloc(2*n*sizeof(double));
    if (!work_arr) { *m = -1; return; }
    double* absA = PyMem_RawMalloc(n*n*sizeof(double));
    if (!absA) { *m = -2; PyMem_RawFree(work_arr); return; }

    dims[0] = n;
    dims[1] = n;
    theta[0] = 1.495585217958292e-002;
    theta[1] = 2.539398330063230e-001;
    theta[2] = 9.504178996162932e-001;
    theta[3] = 2.097847961257068e+000;
    theta[4] = 4.250000000000000e+000;
    coeff[0] = 1.1191048088221578e-11;
    coeff[1] = 1.1167770708198077e-06;
    coeff[2] = 4.9826124310493469e-01;
    coeff[3] = 6.5662862500000000e+05;
    coeff[4] = 1.2573361865339437e+19;

    // Initialize the first n part of work_arr
    for (i = 0; i < n; i++) { work_arr[i] = 1.0; }

    // absA = np.abs(Am[0, :n, :n])
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            absA[j + i*n] = cabs(Am[j + i*n]);
        }
    }

    // First spin = normest(|A|, m=1)
    dgemv_("N", &n, &n, &dbl1, absA, &n, work_arr, &int1, &dbl0, &work_arr[n], &int1);
    normA = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > normA) { normA = work_arr[n+i]; } }


    zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[0*n*n], &n, &Am[0*n*n], &n, &cdbl0, &Am[1*n*n], &n);
    zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[1*n*n], &n, &Am[1*n*n], &n, &cdbl0, &Am[2*n*n], &n);
    zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n*n], &n, &Am[1*n*n], &n, &cdbl0, &Am[3*n*n], &n);
    d4 = pow(znorm1(&Am[2*n*n], work_arr, n), 0.25);
    d6 = pow(znorm1(&Am[3*n*n], work_arr, n), 1.0/6.0);
    eta0 = fmax(d4, d6);
    eta1 = eta0;

    // m = 3
    // -------
    // 1-norm of A**7
    // Alternating matvecs
    // absA * work_arr[n:] = work_arr[:n]
    // absA * work_arr[:n] = work_arr[n:]
    for (i = 0; i < 3; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[0])/6);
    lm = (lm < 0 ? 0 : lm);
    if ((eta0 < theta[0]) && lm == 0)
    {
        *m = 3;
        goto FREE_EXIT;
    }

    // m = 5
    // -------
    // 1-norm of A**11
    for (i = 0; i < 2; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[1])/10);
    lm = (lm < 0 ? 0 : lm);

    if ((eta1 < theta[1]) && lm == 0)
    {
        *m = 5;
        goto FREE_EXIT;
    }

    // m = 7
    // -------
    if (n < 400)
    {
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &cdbl0, &Am[4*n*n], &n);
        d8 = pow(znorm1(&Am[4*n*n], work_arr, n), 0.125);
    } else {
        test = znorm1est(&Am[0], 8);
        // If memory error in z1normest
        if (test <= -100.0) { *m = -3; goto FREE_EXIT; }
        d8 = pow(test, 0.125);
    }

    eta2 = fmax(d6, d8);

    // 1-norm of A**15
    for (i = 0; i < 2; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[2])/14);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[2]) && lm == 0)
    {
        *m = 7;
        goto FREE_EXIT;
    }

    // m = 9
    // -------
    // 1-norm of A**19
    for (i = 0; i < 2; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[3])/18);
    lm = (lm < 0 ? 0 : lm);

    if ((eta2 < theta[3]) && lm == 0)
    {
        // Add the deferred matmul for m=9
        if (n >= 400)
        {
            zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n*n], &n, &Am[2*n*n], &n, &cdbl0, &Am[4*n*n], &n);
        }
        *m = 9;
        goto FREE_EXIT;
    }

    // m = 13
    // -------
    // Scale-square
    if (n < 400)
    {
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[3*n*n], &n, &Am[2*n*n], &n, &cdbl0, &Am[4*n*n], &n);
        d10 = pow(znorm1(&Am[4*n*n], work_arr, n), 0.1);
    } else {
        test = znorm1est(&Am[0], 10);
        // If memory error in z1normest
        if (test <= -100.0) { *m = -4; goto FREE_EXIT; }
        d10 = pow(test, 0.1);
    }

    eta3 = fmax(d8, d10);
    eta4 = fmin(eta2, eta3);
    *s = (int)ceil(log2(eta4 / theta[4]));
    *s = (*s < 0 ? 0 : *s);

    if (*s != 0)
    {
        two_pow_s = pow(2.0, -(*s));
        for (i = 0; i < n*n; i++) { absA[i] *= two_pow_s; }

        // work_arr has spun 19 times already
        two_pow_s = pow(2.0, -19.0*(*s));
        for (i = 0; i < n; i++) { work_arr[n + i] *= two_pow_s; }
        normA *= pow(2.0, -(*s));
    }

    // 1-norm of A**27
    for (i = 0; i < 4; i++)
    {
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[n], &int1, &dbl0, &work_arr[0], &int1);
        dgemv_("N", &n, &n, &dbl1, absA, &n, &work_arr[0], &int1, &dbl0, &work_arr[n], &int1);
    }
    temp = 0.0;
    for (i = 0; i < n; i++)  { if (work_arr[n+i] > temp) { temp = work_arr[n+i]; } }

    lm = (int)ceil(log2(temp/normA/coeff[4])/26.0);
    *s += (lm < 0 ? 0 : lm );
    *m = 13;

    if (*s != 0)
    {
        double s0 = pow(2.0, -*s);
        double s1 = pow(4.0, -*s);
        double s2 = pow(16.0, -*s);
        double s3 = pow(64.0, -*s);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                Am[        n*i + j] = _Cmulcr(Am[        n*i + j], s0);
                Am[  n*n + n*i + j] = _Cmulcr(Am[  n*n + n*i + j], s1);
                Am[2*n*n + n*i + j] = _Cmulcr(Am[2*n*n + n*i + j], s2);
                Am[3*n*n + n*i + j] = _Cmulcr(Am[3*n*n + n*i + j], s3);
#else
                Am[        n*i + j] *= s0;
                Am[  n*n + n*i + j] *= s1;
                Am[2*n*n + n*i + j] *= s2;
                Am[3*n*n + n*i + j] *= s3;
#endif
            }
        }
    }

FREE_EXIT:
    PyMem_RawFree(absA);
    PyMem_RawFree(work_arr);
    return;
}


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************/

static void
pade_UV_calc_s(float* Am, const Py_ssize_t size_n, const int m, int* info)
{
    // U, V computation.
    // We utilize the unused powers of Am as scratch memory depending on m value.
    // Final U is stored Am[3], final V is stored Am[1] for all m.
    // b[m] = 1.0 for all m.

    int int1 = 1, n = (int)size_n;
    int n2 = n*n, i, j;
    float b[13];
    float two = 2.0, dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0;

    int* ipiv = PyMem_RawMalloc(n*sizeof(int));
    if (!ipiv) { *info = -11; return; }

    if (m == 3)
    {
        // Am[2], Am[3] and Am[4] are free
        b[0] = 120.0;
        b[1] = 60.0;
        b[2] = 12.0;

        // U = Am[0] @ Am[1] + b[1]*Am[0]
        scopy_(&n2, Am, &int1, &Am[3*n2], &int1);
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[n2], &n, Am, &n, &b[1], &Am[3*n2], &n);

        // V = b[2]*Am[1] + b[0]*I_n
        sscal_(&n2, &b[2], &Am[n2], &int1);
        for (i = 0; i < n; i++) { Am[n2 + n*i + i] += b[0]; }

    } else if (m == 5) {
        // Am[3] and Am[4] are free
        b[0] = 30240.0;
        b[1] = 15120.0;
        b[2] = 3360.0;
        b[3] = 420.0;
        b[4] = 30.0;

        // U = Am[0] @ (b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        scopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[4*n2 + n*i + j] = Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
            }
        }
        // Am[0] @ (1.0*Am[2] + 420*Am[1] + 15120*I)
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, Am, &n, &dbl0, &Am[3*n2], &n);

        // V
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[n2 + n*i + j] = b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
            }
        }


    } else if (m == 7) {
        // Am[4] is free
        b[0] = 17297280.0;
        b[1] = 8648640.0;
        b[2] = 1995840.0;
        b[3] = 277200.0;
        b[4] = 25200.0;
        b[5] = 1512.0;
        b[6] = 56.0;

        // U = Am[0] @ (b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        scopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[4*n2 + n*i + j] = Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
            }
        }
        // Defer the Am[0] multiplication since we still need Am[3] for V

        // V
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[n2 + n*i + j] = b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
            }
        }

        // Now overwrite Am[3] with U
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, Am, &n, &dbl0, &Am[3*n2], &n);

    } else if (m == 9) {

        b[0] = 17643225600.0;
        b[1] = 8821612800.0;
        b[2] = 2075673600.0;
        b[3] = 302702400.0;
        b[4] = 30270240.0;
        b[5] = 2162160.0;
        b[6] = 110880.0;
        b[7] = 3960.0;
        b[8] = 90.0;

        // U = Am[0] @ (b[9]*Am[4] + b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[8]*Am[4] + b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // We form U and V together since all powers are used simultaneously.
        // U in Am[4] before matmul and V in Am[1]
        float temp1, temp2;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                temp1 = Am[4*n2 + n*i + j];
                temp2 = Am[n2 + n*i + j];
                Am[4*n2 + n*i + j] = temp1 + b[7]*Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*temp2;
                Am[n2 + n*i + j] = b[8]*temp1 + b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*temp2;
                if (i == j)
                {
                    Am[4*n2 + n*i + j] += b[1];
                    Am[n2 + n*i + j] += b[0];
                }
            }
        }
        // Now overwrite Am[3] with U
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, Am, &n, &dbl0, &Am[3*n2], &n);

    } else if (m == 13) {
        b[0]  = 64764752532480000.0;
        b[1]  = 32382376266240000.0;
        b[2]  = 7771770303897600.0;
        b[3]  = 1187353796428800.0;
        b[4]  = 129060195264000.0;
        b[5]  = 10559470521600.0;
        b[6]  = 670442572800.0;
        b[7]  = 33522128640.0;
        b[8]  = 1323241920.0;
        b[9]  = 40840800.0;
        b[10] = 960960.0;
        b[11] = 16380.0;
        b[12] = 182.0;

        // U = Am[0] @ (Am[3] @ (b[13]*Am[3] + b[11]*Am[2] + b[9]*Am[1]) +
        //              b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = Am[3] @ (b[12]*Am[3] + b[10]*Am[2] + b[8]*Am[1]) +
        //              b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // In simpler notation, both are of the form.
        //          U = K @ (L @ M + N)
        //          V =      P @ Q + R
        // So we need A[0], A[3] all the way. We can scratch A[1], A[2] and A[4]
        // but it will not be enough and hence the costly malloc.

        // K = Am[0], L = Am[3], M = work, N = Am[2]
        // P = Am[3], Q = Am[4], R = Am[1]

        float* work = PyMem_RawMalloc(n2*sizeof(float));
        if (!work) { *info = -12; return; }
        float temp1, temp2, temp3;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                temp1 = Am[1*n2 + n*i + j];
                temp2 = Am[2*n2 + n*i + j];
                temp3 = Am[3*n2 + n*i + j];

                Am[2*n2 + n*i + j] =  b[7]*temp3 +  b[5]*temp2 + b[3]*temp1;
                Am[4*n2 + n*i + j] = b[12]*temp3 + b[10]*temp2 + b[8]*temp1;
                Am[1*n2 + n*i + j] =  b[6]*temp3 +  b[4]*temp2 + b[2]*temp1;
                work[n*i + j]      =       temp3 + b[11]*temp2 + b[9]*temp1;
                if (i == j)
                {
                    Am[2*n2 + n*i + j] += b[1];
                    Am[1*n2 + n*i + j] += b[0];
                }
            }
        }

        // V = P @ Q + R
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, &Am[3*n2], &n, &dbl1, &Am[1*n2], &n);

        // U = K @ (L @ M + N)
        sgemm_("N", "N", &n, &n, &n, &dbl1, work, &n, &Am[3*n2], &n, &dbl1, &Am[2*n2], &n);
        sgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n2], &n, Am, &n, &dbl1, &Am[3*n2], &n);

        PyMem_RawFree(work);
    }

    // inv(V-U) (V+U) = inv(V-U) (V-U+2U) = I + 2 inv(V-U) U
    // subtract U from V
    saxpy_(&n2, &dblm1, &Am[3*n2], &int1, &Am[n2], &int1);

    // Solve for (V-U) X = U
    swap_cf_s(&Am[3*n2], &Am[2*n2], n, n, n);

    sgetrf_(&n, &n, &Am[n2], &n, ipiv, info);
    sgetrs_("T", &n, &n, &Am[n2], &n, ipiv, &Am[2*n2], &n, info);
    PyMem_RawFree(ipiv);
    sscal_(&n2, &two, &Am[2*n2], &int1);
    for (i = 0; i < n; i++) { Am[2*n2 + n*i + i] += 1.0; }

    swap_cf_s(&Am[2*n2], Am, n, n, n);

    return;
}


static void
pade_UV_calc_d(double* Am, const Py_ssize_t size_n, const int m, int* info)
{
    // U, V computation.
    // We utilize the unused powers of Am as scratch memory depending on m value.
    // Final U is stored Am[3], final V is stored Am[1] for all m.
    // b[m] = 1.0 for all m.

    int int1 = 1, n = (int)size_n;
    int n2 = n*n, i, j;
    double b[13];
    double two = 2.0, dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0;

    int* ipiv = PyMem_RawMalloc(n*sizeof(int));
    if (!ipiv) { *info = -11; return; }

    if (m == 3)
    {
        // Am[2], Am[3] and Am[4] are free
        b[0] = 120.0;
        b[1] = 60.0;
        b[2] = 12.0;

        // U = Am[0] @ Am[1] + 60.*Am[0]
        // V = 12.*Am[1] + 120*I_n

        // U = Am[0] @ Am[1] + b[1]*Am[0]
        dcopy_(&n2, Am, &int1, &Am[3*n2], &int1);
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[n2], &n, Am, &n, &b[1], &Am[3*n2], &n);

        // V = b[2]*Am[1] + b[0]*I_n
        dscal_(&n2, &b[2], &Am[n2], &int1);
        for (i = 0; i < n; i++) { Am[n2 + n*i + i] += b[0]; }

    } else if (m == 5) {
        // Am[3] and Am[4] are free
        b[0] = 30240.0;
        b[1] = 15120.0;
        b[2] = 3360.0;
        b[3] = 420.0;
        b[4] = 30.0;

        // U = Am[0] @ (b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        dcopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[4*n2 + n*i + j] = Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
            }
        }
        // Am[0] @ (1.0*Am[2] + 420*Am[1] + 15120*I)
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, Am, &n, &dbl0, &Am[3*n2], &n);

        // V
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[n2 + n*i + j] = b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
            }
        }


    } else if (m == 7) {
        // Am[4] is free
        b[0] = 17297280.0;
        b[1] = 8648640.0;
        b[2] = 1995840.0;
        b[3] = 277200.0;
        b[4] = 25200.0;
        b[5] = 1512.0;
        b[6] = 56.0;

        // U = Am[0] @ (b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        dcopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[4*n2 + n*i + j] = Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
            }
        }
        // Defer the Am[0] multiplication since we still need Am[3] for V

        // V
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                Am[n2 + n*i + j] = b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
            }
        }

        // Now overwrite Am[3] with U
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, Am, &n, &dbl0, &Am[3*n2], &n);

    } else if (m == 9) {

        b[0] = 17643225600.0;
        b[1] = 8821612800.0;
        b[2] = 2075673600.0;
        b[3] = 302702400.0;
        b[4] = 30270240.0;
        b[5] = 2162160.0;
        b[6] = 110880.0;
        b[7] = 3960.0;
        b[8] = 90.0;

        // U = Am[0] @ (b[9]*Am[4] + b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[8]*Am[4] + b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // We form U and V together since all powers are used simultaneously.
        // U in Am[4] before matmul and V in Am[1]
        double temp1, temp2;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                temp1 = Am[4*n2 + n*i + j];
                temp2 = Am[n2 + n*i + j];
                Am[4*n2 + n*i + j] = temp1 + b[7]*Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*temp2;
                Am[n2 + n*i + j] = b[8]*temp1 + b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*temp2;
                if (i == j)
                {
                    Am[4*n2 + n*i + j] += b[1];
                    Am[n2 + n*i + j] += b[0];
                }
            }
        }
        // Now overwrite Am[3] with U
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, Am, &n, &dbl0, &Am[3*n2], &n);

    } else if (m == 13) {
        b[0]  = 64764752532480000.0;
        b[1]  = 32382376266240000.0;
        b[2]  = 7771770303897600.0;
        b[3]  = 1187353796428800.0;
        b[4]  = 129060195264000.0;
        b[5]  = 10559470521600.0;
        b[6]  = 670442572800.0;
        b[7]  = 33522128640.0;
        b[8]  = 1323241920.0;
        b[9]  = 40840800.0;
        b[10] = 960960.0;
        b[11] = 16380.0;
        b[12] = 182.0;

        // U = Am[0] @ (Am[3] @ (b[13]*Am[3] + b[11]*Am[2] + b[9]*Am[1]) +
        //              b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = Am[3] @ (b[12]*Am[3] + b[10]*Am[2] + b[8]*Am[1]) +
        //              b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // In simpler notation, both are of the form.
        //          U = K @ (L @ M + N)
        //          V =      P @ Q + R
        // So we need A[0], A[3] all the way. We can scratch A[1], A[2] and A[4]
        // but it will not be enough and hence the costly malloc.

        // K = Am[0], L = Am[3], M = work, N = Am[2]
        // P = Am[3], Q = Am[4], R = Am[1]

        double* work = PyMem_RawMalloc(n2*sizeof(double));
        if (!work) { *info = -12; return; }
        double temp1, temp2, temp3;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                temp1 = Am[1*n2 + n*i + j];
                temp2 = Am[2*n2 + n*i + j];
                temp3 = Am[3*n2 + n*i + j];

                Am[2*n2 + n*i + j] =  b[7]*temp3 +  b[5]*temp2 + b[3]*temp1;
                Am[4*n2 + n*i + j] = b[12]*temp3 + b[10]*temp2 + b[8]*temp1;
                Am[1*n2 + n*i + j] =  b[6]*temp3 +  b[4]*temp2 + b[2]*temp1;
                work[n*i + j]      =       temp3 + b[11]*temp2 + b[9]*temp1;
                if (i == j)
                {
                    Am[2*n2 + n*i + j] += b[1];
                    Am[1*n2 + n*i + j] += b[0];
                }
            }
        }

        // V = P @ Q + R
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[4*n2], &n, &Am[3*n2], &n, &dbl1, &Am[n2], &n);

        // U = K @ (L @ M + N)
        dgemm_("N", "N", &n, &n, &n, &dbl1, work, &n, &Am[3*n2], &n, &dbl1, &Am[2*n2], &n);
        dgemm_("N", "N", &n, &n, &n, &dbl1, &Am[2*n2], &n, Am, &n, &dbl1, &Am[3*n2], &n);

        PyMem_RawFree(work);
    }

    // inv(V-U) (V+U) = inv(V-U) (V-U+2U) = I + 2 inv(V-U) U
    // subtract U from V
    daxpy_(&n2, &dblm1, &Am[3*n2], &int1, &Am[n2], &int1);

    // Solve for (V-U) X = U
    swap_cf_d(&Am[3*n2], &Am[2*n2], n, n, n);

    dgetrf_(&n, &n, &Am[n2], &n, ipiv, info);
    dgetrs_("T", &n, &n, &Am[n2], &n, ipiv, &Am[2*n2], &n, info);
    PyMem_RawFree(ipiv);
    dscal_(&n2, &two, &Am[2*n2], &int1);
    for (i = 0; i < n; i++) { Am[2*n2 + n*i + i] += 1.0; }

    swap_cf_d(&Am[2*n2], Am, n, n, n);
    return;
}


static void
pade_UV_calc_c(SCIPY_C* Am, const Py_ssize_t size_n, const int m, int* info)
{
    // U, V computation.
    // We utilize the unused powers of Am as scratch memory depending on m value.
    // Final U is stored Am[3], final V is stored Am[1] for all m.
    // b[m] = 1.0 for all m.

    int int1 = 1, n = (int)size_n;
    int n2 = n*n, i, j;
    float b[13];
    float two = 2.0;
    SCIPY_C temp1, temp2, temp3, cb1;
    SCIPY_C cdbl1 = CPLX_C(1.0, 0.0), cdbl0 = CPLX_C(0.0, 0.0), cdblm1 = CPLX_C(-1.0, 0.0);

#if defined(_MSC_VER)
    SCIPY_C inter1, inter2, inter3, inter4;
#endif

    int* ipiv = PyMem_RawMalloc(n*sizeof(int));
    if (!ipiv) { *info = -11; return; }

    if (m == 3)
    {
        // Am[2], Am[3] and Am[4] are free
        b[0] = 120.0;
        b[1] = 60.0;
        b[2] = 12.0;
        cb1 = CPLX_C(b[1], 0.0);

        // U = Am[0] @ Am[1] + 60.*Am[0]
        // V = 12.*Am[1] + 120*I_n

        // U = Am[0] @ Am[1] + b[1]*Am[0]
        ccopy_(&n2, Am, &int1, &Am[3*n2], &int1);
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[n2], &n, Am, &n, &cb1, &Am[3*n2], &n);

        // V = b[2]*Am[1] + b[0]*I_n
        csscal_(&n2, &b[2], &Am[n2], &int1);
#if defined(_MSC_VER)
            for (i = 0; i < n; i++) {
                Am[n2 + n*i + i] = CPLX_C(crealf(Am[n2 + n*i + i]) + b[0], cimagf(Am[n2 + n*i + i]));
            }
#else
            for (i = 0; i < n; i++) { Am[n2 + n*i + i] += b[0]; }
#endif

    } else if (m == 5) {
        // Am[3] and Am[4] are free
        b[0] = 30240.0;
        b[1] = 15120.0;
        b[2] = 3360.0;
        b[3] = 420.0;
        b[4] = 30.0;

        // U = Am[0] @ (b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        ccopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _FCmulcr(Am[n2 + n*i + j], b[3]);
                inter2 = Am[2*n2 + n*i + j];
                Am[4*n2 + n*i + j] =  CPLX_C(crealf(inter1) + crealf(inter2), cimagf(inter1)+ cimagf(inter2));
                if (i == j) { Am[4*n2 + n*i + j] = CPLX_C(crealf(Am[4*n2 + n*i + j]) + b[1], cimagf(Am[4*n2 + n*i + j])); }
#else
                Am[4*n2 + n*i + j] = Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
#endif
            }
        }

        // Am[0] @ (1.0*Am[2] + 420*Am[1] + 15120*I)
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

        // V
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _FCmulcr(Am[n2 + n*i + j], b[2]);
                inter2 = _FCmulcr(Am[2*n2 + n*i + j], b[4]);
                Am[n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2), cimagf(inter1)+ cimagf(inter2));
                if (i == j) { Am[n2 + n*i + j] = CPLX_C(crealf(Am[n2 + n*i + j]) + b[0], cimagf(Am[n2 + n*i + j])); }
#else
                Am[n2 + n*i + j] = b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
#endif
            }
        }

    } else if (m == 7) {
        // Am[4] is free
        b[0] = 17297280.0;
        b[1] = 8648640.0;
        b[2] = 1995840.0;
        b[3] = 277200.0;
        b[4] = 25200.0;
        b[5] = 1512.0;
        b[6] = 56.0;

        // U = Am[0] @ (b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        ccopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _FCmulcr(Am[n2 + n*i + j], b[3]);
                inter2 = _FCmulcr(Am[2*n2 + n*i + j], b[5]);
                inter3 = Am[3*n2 + n*i + j];
                Am[4*n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3),
                                            cimagf(inter1) + cimagf(inter2) + cimagf(inter3));
                if (i == j) { Am[4*n2 + n*i + j]  = CPLX_C(crealf(Am[4*n2 + n*i + j]) + b[1], cimagf(Am[4*n2 + n*i + j])); }
#else
                Am[4*n2 + n*i + j] = Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
#endif
            }
        }

        // Defer the Am[0] multiplication since we still need Am[3] for V
        // V

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _FCmulcr(Am[n2 + n*i + j], b[2]);
                inter2 = _FCmulcr(Am[2*n2 + n*i + j], b[4]);
                inter3 = _FCmulcr(Am[3*n2 + n*i + j], b[6]);
                Am[n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3),
                                          cimagf(inter1) + cimagf(inter2) + cimagf(inter3));
                if (i == j) { Am[n2 + n*i + j]  = CPLX_C(crealf(Am[n2 + n*i + j]) + b[0], cimagf(Am[n2 + n*i + j])); }
#else
                Am[n2 + n*i + j] = b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
#endif
            }
        }

        // Now overwrite Am[3] with U
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

    } else if (m == 9) {

        b[0] = 17643225600.0;
        b[1] = 8821612800.0;
        b[2] = 2075673600.0;
        b[3] = 302702400.0;
        b[4] = 30270240.0;
        b[5] = 2162160.0;
        b[6] = 110880.0;
        b[7] = 3960.0;
        b[8] = 90.0;

        // U = Am[0] @ (b[9]*Am[4] + b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[8]*Am[4] + b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // We form U and V together since all powers are used simultaneously.
        // U in Am[4] before matmul and V in Am[1]

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                temp1 = Am[4*n2 + n*i + j];
                temp2 = Am[n2 + n*i + j];

                inter2 = _FCmulcr(temp2, b[3]);
                inter3 = _FCmulcr(Am[2*n2 + n*i + j], b[5]);
                inter4 = _FCmulcr(Am[3*n2 + n*i + j], b[7]);

                Am[4*n2 + n*i + j] = CPLX_C(crealf(temp1) + crealf(inter2) + crealf(inter3) + crealf(inter4),
                                            cimagf(temp1) + cimagf(inter2) + cimagf(inter3) + cimagf(inter4));

                inter1 = _FCmulcr(temp1, b[8]);
                inter2 = _FCmulcr(temp2, b[2]);
                inter3 = _FCmulcr(Am[2*n2 + n*i + j], b[4]);
                inter4 = _FCmulcr(Am[3*n2 + n*i + j], b[6]);

                Am[n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3) + crealf(inter4),
                                          cimagf(inter1) + cimagf(inter2) + cimagf(inter3) + cimagf(inter4));
                if (i == j)
                {
                    Am[4*n2 + n*i + j] = CPLX_C(crealf(Am[4*n2 + n*i + j]) + b[1], cimagf(Am[4*n2 + n*i + j]));
                    Am[n2 + n*i + j] = CPLX_C(crealf(Am[n2 + n*i + j]) + b[0], cimagf(Am[n2 + n*i + j]));
                }
#else
                temp1 = Am[4*n2 + n*i + j];
                temp2 = Am[n2 + n*i + j];
                Am[4*n2 + n*i + j] = temp1 + b[7]*Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*temp2;
                Am[n2 + n*i + j] = b[8]*temp1 + b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*temp2;
                if (i == j)
                {
                    Am[4*n2 + n*i + j] += b[1];
                    Am[n2 + n*i + j] += b[0];
                }
#endif
            }
        }

        // Now overwrite Am[3] with U
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

    } else if (m == 13) {
        b[0]  = 64764752532480000.0;
        b[1]  = 32382376266240000.0;
        b[2]  = 7771770303897600.0;
        b[3]  = 1187353796428800.0;
        b[4]  = 129060195264000.0;
        b[5]  = 10559470521600.0;
        b[6]  = 670442572800.0;
        b[7]  = 33522128640.0;
        b[8]  = 1323241920.0;
        b[9]  = 40840800.0;
        b[10] = 960960.0;
        b[11] = 16380.0;
        b[12] = 182.0;

        // U = Am[0] @ (Am[3] @ (b[13]*Am[3] + b[11]*Am[2] + b[9]*Am[1]) +
        //              b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = Am[3] @ (b[12]*Am[3] + b[10]*Am[2] + b[8]*Am[1]) +
        //              b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // In simpler notation, both are of the form.
        //          U = K @ (L @ M + N)
        //          V =      P @ Q + R
        // So we need A[0], A[3] all the way. We can scratch A[1], A[2] and A[4]
        // but it will not be enough and hence the costly malloc.

        // K = Am[0], L = Am[3], M = work, N = Am[2]
        // P = Am[3], Q = Am[4], R = Am[1]

        SCIPY_C* work = PyMem_RawMalloc(n2*sizeof(SCIPY_C));
        if (!work) { *info = -12; return; }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                temp1 = Am[1*n2 + n*i + j];
                temp2 = Am[2*n2 + n*i + j];
                temp3 = Am[3*n2 + n*i + j];

                inter1 = _FCmulcr(temp1, b[3]);
                inter2 = _FCmulcr(temp2, b[5]);
                inter3 = _FCmulcr(temp3, b[7]);

                Am[2*n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3),
                                            cimagf(inter1) + cimagf(inter2) + cimagf(inter3));

                inter1 = _FCmulcr(temp1, b[8]);
                inter2 = _FCmulcr(temp2, b[10]);
                inter3 = _FCmulcr(temp3, b[12]);

                Am[4*n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3),
                                            cimagf(inter1) + cimagf(inter2) + cimagf(inter3));

                inter1 = _FCmulcr(temp1, b[2]);
                inter2 = _FCmulcr(temp2, b[4]);
                inter3 = _FCmulcr(temp3, b[6]);

                Am[1*n2 + n*i + j] = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3),
                                            cimagf(inter1) + cimagf(inter2) + cimagf(inter3));

                inter1 = _FCmulcr(temp1, b[9]);
                inter2 = _FCmulcr(temp2, b[11]);
                inter3 = temp3;

                work[n*i + j]      = CPLX_C(crealf(inter1) + crealf(inter2) + crealf(inter3),
                                            cimagf(inter1) + cimagf(inter2) + cimagf(inter3));
                if (i == j)
                {
                    Am[2*n2 + n*i + j] = CPLX_C(crealf(Am[2*n2 + n*i + j]) + b[1], cimagf(Am[2*n2 + n*i + j]));
                    Am[1*n2 + n*i + j] = CPLX_C(crealf(Am[1*n2 + n*i + j]) + b[0], cimagf(Am[1*n2 + n*i + j]));
                }
#else
                temp1 = Am[1*n2 + n*i + j];
                temp2 = Am[2*n2 + n*i + j];
                temp3 = Am[3*n2 + n*i + j];

                Am[2*n2 + n*i + j] =  b[7]*temp3 +  b[5]*temp2 + b[3]*temp1;
                Am[4*n2 + n*i + j] = b[12]*temp3 + b[10]*temp2 + b[8]*temp1;
                Am[1*n2 + n*i + j] =  b[6]*temp3 +  b[4]*temp2 + b[2]*temp1;
                work[n*i + j]      =       temp3 + b[11]*temp2 + b[9]*temp1;
                if (i == j)
                {
                    Am[2*n2 + n*i + j] += b[1];
                    Am[1*n2 + n*i + j] += b[0];
                }
#endif
            }
        }

        // V = P @ Q + R
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, &Am[3*n2], &n, &cdbl1, &Am[1*n2], &n);

        // U = K @ (L @ M + N)
        cgemm_("N", "N", &n, &n, &n, &cdbl1, work, &n, &Am[3*n2], &n, &cdbl1, &Am[2*n2], &n);
        cgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

        PyMem_RawFree(work);
    }

    // inv(V-U) (V+U) = inv(V-U) (V-U+2U) = I + 2 inv(V-U) U
    // subtract U from V
    caxpy_(&n2, &cdblm1, &Am[3*n2], &int1, &Am[n2], &int1);

    // Solve for (V-U) X = U
    swap_cf_c(&Am[3*n2], &Am[2*n2], n, n, n);

    cgetrf_(&n, &n, &Am[n2], &n, ipiv, info);
    cgetrs_("T", &n, &n, &Am[n2], &n, ipiv, &Am[2*n2], &n, info);
    PyMem_RawFree(ipiv);
    csscal_(&n2, &two, &Am[2*n2], &int1);

#if defined(_MSC_VER)
    for (i = 0; i < n; i++) { Am[2*n2 + n*i + i] = CPLX_C(crealf(Am[2*n2 + n*i + i]) + 1.0f, cimagf(Am[2*n2 + n*i + i])); }
#else
    for (i = 0; i < n; i++) { Am[2*n2 + n*i + i] += cdbl1; }
#endif

    swap_cf_c(&Am[2*n2], Am, n, n, n);

    return;
}


static void
pade_UV_calc_z(SCIPY_Z* Am, const Py_ssize_t size_n, const int m, int* info)
{
    // U, V computation.
    // We utilize the unused powers of Am as scratch memory depending on m value.
    // Final U is stored Am[3], final V is stored Am[1] for all m.
    // b[m] = 1.0 for all m.

    int int1 = 1, n = (int)size_n;
    int n2 = n*n, i, j;
    double b[13];
    double two = 2.0;
    SCIPY_Z temp1, temp2, temp3, cb1;
    SCIPY_Z cdbl1 = CPLX_Z(1.0, 0.0), cdbl0 = CPLX_Z(0.0, 0.0), cdblm1 = CPLX_Z(-1.0, 0.0);

#if defined(_MSC_VER)
    SCIPY_Z inter1, inter2, inter3, inter4;
#endif

    int* ipiv = PyMem_RawMalloc(n*sizeof(int));
    if (!ipiv) { *info = -11; return; }

    if (m == 3)
    {
        // Am[2], Am[3] and Am[4] are free
        b[0] = 120.0;
        b[1] = 60.0;
        b[2] = 12.0;
        cb1 = CPLX_Z(b[1], 0.0);

        // U = Am[0] @ Am[1] + b[1]*Am[0]
        zcopy_(&n2, Am, &int1, &Am[3*n2], &int1);
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[n2], &n, Am, &n, &cb1, &Am[3*n2], &n);

        // V = b[2]*Am[1] + b[0]*I_n
        zdscal_(&n2, &b[2], &Am[n2], &int1);
#if defined(_MSC_VER)
        for (i = 0; i < n; i++) { Am[n2 + n*i + i] = CPLX_Z(creal(Am[n2 + n*i + i]) + b[0], cimag(Am[n2 + n*i + i])); }
#else
        for (i = 0; i < n; i++) { Am[n2 + n*i + i] += b[0]; }
#endif

    } else if (m == 5) {
        // Am[3] and Am[4] are free
        b[0] = 30240.0;
        b[1] = 15120.0;
        b[2] = 3360.0;
        b[3] = 420.0;
        b[4] = 30.0;

        // U = Am[0] @ (b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        zcopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _Cmulcr(Am[n2 + n*i + j], b[3]);
                inter2 = Am[2*n2 + n*i + j];
                Am[4*n2 + n*i + j] =  CPLX_Z(creal(inter1) + creal(inter2), cimag(inter1)+ cimag(inter2));
                if (i == j) { Am[4*n2 + n*i + j] = CPLX_Z(creal(Am[4*n2 + n*i + j]) + b[1], cimag(Am[4*n2 + n*i + j])); }
#else
                Am[4*n2 + n*i + j] = Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
#endif
            }
        }
        // Am[0] @ (1.0*Am[2] + 420*Am[1] + 15120*I)
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

        // V
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _Cmulcr(Am[n2 + n*i + j], b[2]);
                inter2 = _Cmulcr(Am[2*n2 + n*i + j], b[4]);
                Am[n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2), cimag(inter1)+ cimag(inter2));
                if (i == j) { Am[n2 + n*i + j] = CPLX_Z(creal(Am[n2 + n*i + j]) + b[0], cimag(Am[n2 + n*i + j])); }
#else
                Am[n2 + n*i + j] = b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
#endif
            }
        }

    } else if (m == 7) {
        // Am[4] is free
        b[0] = 17297280.0;
        b[1] = 8648640.0;
        b[2] = 1995840.0;
        b[3] = 277200.0;
        b[4] = 25200.0;
        b[5] = 1512.0;
        b[6] = 56.0;

        // U = Am[0] @ (b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // Am[4] = Am[1]
        zcopy_(&n2, &Am[n2], &int1, &Am[4*n2], &int1);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _Cmulcr(Am[n2 + n*i + j], b[3]);
                inter2 = _Cmulcr(Am[2*n2 + n*i + j], b[5]);
                inter3 = Am[3*n2 + n*i + j];
                Am[4*n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3),
                                            cimag(inter1) + cimag(inter2) + cimag(inter3));
                if (i == j) { Am[4*n2 + n*i + j]  = CPLX_Z(creal(Am[4*n2 + n*i + j]) + b[1], cimag(Am[4*n2 + n*i + j])); }
#else
                Am[4*n2 + n*i + j] = Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*Am[n2 + n*i + j];
                if (i == j) { Am[4*n2 + n*i + j] += b[1]; }
#endif
            }
        }

        // Defer the Am[0] multiplication since we still need Am[3] for V
        // V

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                inter1 = _Cmulcr(Am[n2 + n*i + j], b[2]);
                inter2 = _Cmulcr(Am[2*n2 + n*i + j], b[4]);
                inter3 = _Cmulcr(Am[3*n2 + n*i + j], b[6]);
                Am[n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3),
                                          cimag(inter1) + cimag(inter2) + cimag(inter3));
                if (i == j) { Am[n2 + n*i + j]  = CPLX_Z(creal(Am[n2 + n*i + j]) + b[0], cimag(Am[n2 + n*i + j])); }
#else
                Am[n2 + n*i + j] = b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*Am[n2 + n*i + j];
                if (i == j) { Am[n2 + n*i + j] += b[0]; }
#endif
            }
        }

        // Now overwrite Am[3] with U
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

    } else if (m == 9) {

        b[0] = 17643225600.0;
        b[1] = 8821612800.0;
        b[2] = 2075673600.0;
        b[3] = 302702400.0;
        b[4] = 30270240.0;
        b[5] = 2162160.0;
        b[6] = 110880.0;
        b[7] = 3960.0;
        b[8] = 90.0;

        // U = Am[0] @ (b[9]*Am[4] + b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = b[8]*Am[4] + b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // We form U and V together since all powers are used simultaneously.
        // U in Am[4] before matmul and V in Am[1]

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                temp1 = Am[4*n2 + n*i + j];
                temp2 = Am[n2 + n*i + j];

                inter2 = _Cmulcr(temp2, b[3]);
                inter3 = _Cmulcr(Am[2*n2 + n*i + j], b[5]);
                inter4 = _Cmulcr(Am[3*n2 + n*i + j], b[7]);

                Am[4*n2 + n*i + j] = CPLX_Z(creal(temp1) + creal(inter2) + creal(inter3) + creal(inter4),
                                            cimag(temp1) + cimag(inter2) + cimag(inter3) + cimag(inter4));

                inter1 = _Cmulcr(temp1, b[8]);
                inter2 = _Cmulcr(temp2, b[2]);
                inter3 = _Cmulcr(Am[2*n2 + n*i + j], b[4]);
                inter4 = _Cmulcr(Am[3*n2 + n*i + j], b[6]);

                Am[n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3) + creal(inter4),
                                          cimag(inter1) + cimag(inter2) + cimag(inter3) + cimag(inter4));
                if (i == j)
                {
                    Am[4*n2 + n*i + j] = CPLX_Z(creal(Am[4*n2 + n*i + j]) + b[1], cimag(Am[4*n2 + n*i + j]));
                    Am[n2 + n*i + j] = CPLX_Z(creal(Am[n2 + n*i + j]) + b[0], cimag(Am[n2 + n*i + j]));
                }
#else
                temp1 = Am[4*n2 + n*i + j];
                temp2 = Am[n2 + n*i + j];
                Am[4*n2 + n*i + j] = temp1 + b[7]*Am[3*n2 + n*i + j] + b[5]*Am[2*n2 + n*i + j] + b[3]*temp2;
                Am[n2 + n*i + j] = b[8]*temp1 + b[6]*Am[3*n2 + n*i + j] + b[4]*Am[2*n2 + n*i + j] + b[2]*temp2;
                if (i == j)
                {
                    Am[4*n2 + n*i + j] += b[1];
                    Am[n2 + n*i + j] += b[0];
                }
#endif
            }
        }

        // Now overwrite Am[3] with U
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

    } else if (m == 13) {
        b[0]  = 64764752532480000.0;
        b[1]  = 32382376266240000.0;
        b[2]  = 7771770303897600.0;
        b[3]  = 1187353796428800.0;
        b[4]  = 129060195264000.0;
        b[5]  = 10559470521600.0;
        b[6]  = 670442572800.0;
        b[7]  = 33522128640.0;
        b[8]  = 1323241920.0;
        b[9]  = 40840800.0;
        b[10] = 960960.0;
        b[11] = 16380.0;
        b[12] = 182.0;

        // U = Am[0] @ (Am[3] @ (b[13]*Am[3] + b[11]*Am[2] + b[9]*Am[1]) +
        //              b[7]*Am[3] + b[5]*Am[2] + b[3]*Am[1] + b[1]*I_n)
        // V = Am[3] @ (b[12]*Am[3] + b[10]*Am[2] + b[8]*Am[1]) +
        //              b[6]*Am[3] + b[4]*Am[2] + b[2]*Am[1] + b[0]*I_n

        // In simpler notation, both are of the form.
        //          U = K @ (L @ M + N)
        //          V =      P @ Q + R
        // So we need A[0], A[3] all the way. We can scratch A[1], A[2] and A[4]
        // but it will not be enough and hence the costly malloc.

        // K = Am[0], L = Am[3], M = work, N = Am[2]
        // P = Am[3], Q = Am[4], R = Am[1]

        SCIPY_Z* work = PyMem_RawMalloc(n2*sizeof(SCIPY_Z));
        if (!work) { *info = -12; return; }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
#if defined(_MSC_VER)
                temp1 = Am[1*n2 + n*i + j];
                temp2 = Am[2*n2 + n*i + j];
                temp3 = Am[3*n2 + n*i + j];

                inter1 = _Cmulcr(temp1, b[3]);
                inter2 = _Cmulcr(temp2, b[5]);
                inter3 = _Cmulcr(temp3, b[7]);

                Am[2*n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3),
                                            cimag(inter1) + cimag(inter2) + cimag(inter3));

                inter1 = _Cmulcr(temp1, b[8]);
                inter2 = _Cmulcr(temp2, b[10]);
                inter3 = _Cmulcr(temp3, b[12]);

                Am[4*n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3),
                                            cimag(inter1) + cimag(inter2) + cimag(inter3));

                inter1 = _Cmulcr(temp1, b[2]);
                inter2 = _Cmulcr(temp2, b[4]);
                inter3 = _Cmulcr(temp3, b[6]);

                Am[1*n2 + n*i + j] = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3),
                                            cimag(inter1) + cimag(inter2) + cimag(inter3));

                inter1 = _Cmulcr(temp1, b[9]);
                inter2 = _Cmulcr(temp2, b[11]);
                inter3 = temp3;

                work[n*i + j]      = CPLX_Z(creal(inter1) + creal(inter2) + creal(inter3),
                                            cimag(inter1) + cimag(inter2) + cimag(inter3));
                if (i == j)
                {
                    Am[2*n2 + n*i + j] = CPLX_Z(creal(Am[2*n2 + n*i + j]) + b[1], cimag(Am[2*n2 + n*i + j]));
                    Am[1*n2 + n*i + j] = CPLX_Z(creal(Am[1*n2 + n*i + j]) + b[0], cimag(Am[1*n2 + n*i + j]));
                }
#else
                temp1 = Am[1*n2 + n*i + j];
                temp2 = Am[2*n2 + n*i + j];
                temp3 = Am[3*n2 + n*i + j];

                Am[2*n2 + n*i + j] =  b[7]*temp3 +  b[5]*temp2 + b[3]*temp1;
                Am[4*n2 + n*i + j] = b[12]*temp3 + b[10]*temp2 + b[8]*temp1;
                Am[1*n2 + n*i + j] =  b[6]*temp3 +  b[4]*temp2 + b[2]*temp1;
                work[n*i + j]      =       temp3 + b[11]*temp2 + b[9]*temp1;
                if (i == j)
                {
                    Am[2*n2 + n*i + j] += b[1];
                    Am[1*n2 + n*i + j] += b[0];
                }
#endif
            }
        }

        // V = P @ Q + R
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[4*n2], &n, &Am[3*n2], &n, &cdbl1, &Am[1*n2], &n);

        // U = K @ (L @ M + N)
        zgemm_("N", "N", &n, &n, &n, &cdbl1, work, &n, &Am[3*n2], &n, &cdbl1, &Am[2*n2], &n);
        zgemm_("N", "N", &n, &n, &n, &cdbl1, &Am[2*n2], &n, Am, &n, &cdbl0, &Am[3*n2], &n);

        PyMem_RawFree(work);
    }

    // inv(V-U) (V+U) = inv(V-U) (V-U+2U) = I + 2 inv(V-U) U
    // subtract U from V
    zaxpy_(&n2, &cdblm1, &Am[3*n2], &int1, &Am[n2], &int1);
    // Solve for (V-U) X = U
    swap_cf_z(&Am[3*n2], &Am[2*n2], n, n, n);

    zgetrf_(&n, &n, &Am[n2], &n, ipiv, info);
    zgetrs_("T", &n, &n, &Am[n2], &n, ipiv, &Am[2*n2], &n, info);
    PyMem_RawFree(ipiv);
    zdscal_(&n2, &two, &Am[2*n2], &int1);

#if defined(_MSC_VER)
    for (i = 0; i < n; i++) { Am[2*n2 + n*i + i] = CPLX_Z(creal(Am[2*n2 + n*i + i]) + 1.0, cimag(Am[2*n2 + n*i + i])); }
#else
    for (i = 0; i < n; i++) { Am[2*n2 + n*i + i] += cdbl1; }
#endif

    swap_cf_z(&Am[2*n2], &Am[0], n, n, n);

    return;
}


void
matrix_exponential_s(PyArrayObject* a, float* restrict result, int* info)
{
    int m = 0, s = 0, is_lower = 0;
    npy_intp lband = 0, uband = 0;
    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    float* restrict a_data = (float*)PyArray_DATA(a);
    int ndim = PyArray_NDIM(a);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(a);      // Array shape
    int n = (int)shape[ndim - 1];                // Slice size
    npy_intp* restrict strides = PyArray_STRIDES(a);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    float* restrict Am = malloc(sizeof(float)*(n*n*5 + 2*n));
    if (Am == NULL) { *info = -100; return; }
    float* restrict Am1 = &Am[n*n];
    // These two arrays are only used for the triangular case for scaling/squaring
    float* restrict diag_aw = &Am[5*n*n];
    float* restrict sd = &Am[5*n*n + n];

    /*====================================================================
    |                    MAIN nxn SLICE LOOP                             |
    ====================================================================*/
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // See sqrtm for explanation of the loop and offset calculations
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        float* restrict slice_ptr = (float*)(a_data + (offset/sizeof(float)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Am1[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(float)) + (j*strides[ndim - 1]/sizeof(float)));
            }
        }
        bandwidth_s(Am1, n, n, &lband, &uband);
        if ((lband == 0) && (uband == 0)) {
            for (npy_intp i = 0; i < n; i++) {
                result[idx*n*n + i*n + i] = expf(Am1[i*n + i]);
            }
            continue;
        }
        swap_cf_s(Am1, Am, n, n, n);


        is_lower = (uband == 0);
        if ((lband == 0) || (uband == 0)) {
            if (is_lower) {
                // Column-major: iterate columns, grab diagonal and the subdiagonal elements
                for (npy_intp i = 0; i < n - 1; i++) {
                    diag_aw[i] = Am[i * n + i];         // A[i, i]
                    sd[i] = Am[i * n + (i + 1)];        // A[i+1, i]
                }
                diag_aw[n - 1] = Am[(n - 1) * n + (n - 1)];  // last diagonal
            } else {
                // Column-major: iterate columns, grab superdiagonal and diagonal elements
                diag_aw[0] = Am[0];  // first diagonal A[0,0]
                for (npy_intp i = 0; i < n - 1; i++) {
                    sd[i] = Am[(i + 1) * n + i];        // A[i, i+1]
                    diag_aw[i + 1] = Am[(i + 1) * n + (i + 1)];  // A[i+1, i+1]
                }
            }
        }

        pick_pade_structure_s(Am, n, &m, &s);
        if (m < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }
        pade_UV_calc_s(Am, n, m, info);
        if (*info < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }

        float* temp1 = Am;
        float* temp2 = Am1;
        if (s > 0) {
            // Squaring needed.
            if ((lband == 0) || (uband == 0)) {
                // Triangular case - use Fragment 2.1 of Al-Mohy and Higham (2009)
                for (int iter = s - 1; iter >= 0; iter--) {
                    sgemm_("N", "N", &n, &n, &n, &(float){1.0f}, temp1, &n, temp1, &n, &(float){0.0f}, temp2, &n);
                    float* swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;

                    float scale = exp2f(-(float)iter);  // 2^(-iter)

                    if (is_lower) {

                        for (npy_intp i = 0; i < n - 1; i++) {
                            float d_i = diag_aw[i] * scale;
                            float d_ip1 = diag_aw[i + 1] * scale;
                            float diff = d_ip1 - d_i;
                            float exp_sinch = (diff == 0.0f) ? expf(d_i) : (expf(d_ip1) - expf(d_i)) / diff;
                            temp1[i * n + i] = expf(d_i);                       // diagonal
                            temp1[i * n + (i + 1)] = exp_sinch * sd[i] * scale; // subdiagonal
                        }
                        temp1[(n-1) * n + (n-1)] = expf(diag_aw[n-1] * scale);  // last diagonal

                    } else {

                        temp1[0] = expf(diag_aw[0] * scale);                    // first diagonal
                        for (npy_intp i = 0; i < n - 1; i++) {
                            float d_i = diag_aw[i] * scale;
                            float d_ip1 = diag_aw[i + 1] * scale;
                            float diff = d_ip1 - d_i;
                            float exp_sinch = (diff == 0.0f) ? expf(d_i) : (expf(d_ip1) - expf(d_i)) / diff;
                            temp1[(i + 1) * n + i] = exp_sinch * sd[i] * scale; // superdiagonal
                            temp1[(i + 1) * n + (i + 1)] = expf(d_ip1);         // diagonal
                        }

                    }
                }
            } else {
                // General dense case, compute A**(2**s) by repeated squaring.
                for (int i = 0; i < s; i++) {
                    sgemm_("N", "N", &n, &n, &n, &(float){1.0f}, temp1, &n, temp1, &n, &(float){0.0f}, temp2, &n);
                    // Swap pointers
                    float *swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;
                }
            }
        }

        // Copy back the result to output array
        swap_cf_s(temp1, &result[idx*n*n], n, n, n);
    }

    free(Am);
}


void
matrix_exponential_d(PyArrayObject* a, double* restrict result, int* info)
{
    int m = 0, s = 0, is_lower = 0;
    npy_intp lband = 0, uband = 0;
    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    double* restrict a_data = (double*)PyArray_DATA(a);
    int ndim = PyArray_NDIM(a);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(a);      // Array shape
    int n = (int)shape[ndim - 1];                // Slice size
    npy_intp* restrict strides = PyArray_STRIDES(a);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    double* restrict Am = malloc(sizeof(double)*(n*n*5 + 2*n));
    if (Am == NULL) { *info = -100; return; }
    double* restrict Am1 = &Am[n*n];
    // These two arrays are only used for the triangular case for scaling/squaring
    double* restrict diag_aw = &Am[5*n*n];
    double* restrict sd = &Am[5*n*n + n];

    /*====================================================================
    |                    MAIN nxn SLICE LOOP                             |
    ====================================================================*/
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // See sqrtm for explanation of the loop and offset calculations
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        double* restrict slice_ptr = (double*)(a_data + (offset/sizeof(double)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Am1[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(double)) + (j*strides[ndim - 1]/sizeof(double)));
            }
        }
        bandwidth_d(Am1, n, n, &lband, &uband);
        if ((lband == 0) && (uband == 0)) {
            for (npy_intp i = 0; i < n; i++) {
                result[idx*n*n + i*n + i] = exp(Am1[i*n + i]);
            }
            continue;
        }
        swap_cf_d(Am1, Am, n, n, n);


        is_lower = (uband == 0);
        if ((lband == 0) || (uband == 0)) {
            if (is_lower) {
                // Column-major: iterate columns, grab diagonal and the subdiagonal elements
                for (npy_intp i = 0; i < n - 1; i++) {
                    diag_aw[i] = Am[i * n + i];         // A[i, i]
                    sd[i] = Am[i * n + (i + 1)];        // A[i+1, i]
                }
                diag_aw[n - 1] = Am[(n - 1) * n + (n - 1)];  // last diagonal
            } else {
                // Column-major: iterate columns, grab superdiagonal and diagonal elements
                diag_aw[0] = Am[0];  // first diagonal A[0,0]
                for (npy_intp i = 0; i < n - 1; i++) {
                    sd[i] = Am[(i + 1) * n + i];        // A[i, i+1]
                    diag_aw[i + 1] = Am[(i + 1) * n + (i + 1)];  // A[i+1, i+1]
                }
            }
        }

        pick_pade_structure_d(Am, n, &m, &s);
        if (m < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }
        pade_UV_calc_d(Am, n, m, info);
        if (*info < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }

        double* temp1 = Am;
        double* temp2 = Am1;
        if (s > 0) {
            // Squaring needed.
            if ((lband == 0) || (uband == 0)) {
                // Triangular case - use Fragment 2.1 of Al-Mohy and Higham (2009)
                for (int iter = s - 1; iter >= 0; iter--) {
                    dgemm_("N", "N", &n, &n, &n, &(double){1.0}, temp1, &n, temp1, &n, &(double){0.0}, temp2, &n);
                    double* swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;

                    double scale = exp2(-(double)iter);  // 2^(-iter)

                    if (is_lower) {

                        for (npy_intp i = 0; i < n - 1; i++) {
                            double d_i = diag_aw[i] * scale;
                            double d_ip1 = diag_aw[i + 1] * scale;
                            double diff = d_ip1 - d_i;
                            double exp_sinch = (diff == 0.0) ? exp(d_i) : (exp(d_ip1) - exp(d_i)) / diff;
                            temp1[i * n + i] = exp(d_i);                        // diagonal
                            temp1[i * n + (i + 1)] = exp_sinch * sd[i] * scale; // subdiagonal
                        }
                        temp1[(n-1) * n + (n-1)] = exp(diag_aw[n-1] * scale);   // last diagonal

                    } else {

                        temp1[0] = exp(diag_aw[0] * scale);                     // first diagonal
                        for (npy_intp i = 0; i < n - 1; i++) {
                            double d_i = diag_aw[i] * scale;
                            double d_ip1 = diag_aw[i + 1] * scale;
                            double diff = d_ip1 - d_i;
                            double exp_sinch = (diff == 0.0) ? exp(d_i) : (exp(d_ip1) - exp(d_i)) / diff;
                            temp1[(i + 1) * n + i] = exp_sinch * sd[i] * scale; // superdiagonal
                            temp1[(i + 1) * n + (i + 1)] = exp(d_ip1);          // diagonal
                        }

                    }
                }
            } else {
                // General dense case, compute A**(2**s) by repeated squaring.
                for (int i = 0; i < s; i++) {
                    dgemm_("N", "N", &n, &n, &n, &(double){1.0}, temp1, &n, temp1, &n, &(double){0.0}, temp2, &n);
                    // Swap pointers
                    double *swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;
                }
            }
        }

        // Copy back the result to output array
        swap_cf_d(temp1, &result[idx*n*n], n, n, n);
    }

    free(Am);
}


void
matrix_exponential_c(PyArrayObject* a, SCIPY_C* restrict result, int* info)
{
    int m = 0, s = 0, is_lower = 0;
    npy_intp lband = 0, uband = 0;
    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    SCIPY_C* restrict a_data = (SCIPY_C*)PyArray_DATA(a);
    int ndim = PyArray_NDIM(a);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(a);      // Array shape
    int n = (int)shape[ndim - 1];                // Slice size
    npy_intp* restrict strides = PyArray_STRIDES(a);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    SCIPY_C* restrict Am = malloc(sizeof(SCIPY_C)*(n*n*5 + 2*n));
    if (Am == NULL) { *info = -100; return; }
    SCIPY_C* restrict Am1 = &Am[n*n];
    // These two arrays are only used for the triangular case for scaling/squaring
    SCIPY_C* restrict diag_aw = &Am[5*n*n];
    SCIPY_C* restrict sd = &Am[5*n*n + n];

    /*====================================================================
    |                    MAIN nxn SLICE LOOP                             |
    ====================================================================*/
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // See sqrtm for explanation of the loop and offset calculations
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        SCIPY_C* restrict slice_ptr = (SCIPY_C*)(a_data + (offset/sizeof(SCIPY_C)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Am1[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(SCIPY_C)) + (j*strides[ndim - 1]/sizeof(SCIPY_C)));
            }
        }
        bandwidth_c(Am1, n, n, &lband, &uband);
        if ((lband == 0) && (uband == 0)) {
            for (npy_intp i = 0; i < n; i++) {
                result[idx*n*n + i*n + i] = cexpf(Am1[i*n + i]);
            }
            continue;
        }
        swap_cf_c(Am1, Am, n, n, n);


        is_lower = (uband == 0);
        if ((lband == 0) || (uband == 0)) {
            if (is_lower) {
                // Column-major: iterate columns, grab diagonal and the subdiagonal elements
                for (npy_intp i = 0; i < n - 1; i++) {
                    diag_aw[i] = Am[i * n + i];         // A[i, i]
                    sd[i] = Am[i * n + (i + 1)];        // A[i+1, i]
                }
                diag_aw[n - 1] = Am[(n - 1) * n + (n - 1)];  // last diagonal
            } else {
                // Column-major: iterate columns, grab superdiagonal and diagonal elements
                diag_aw[0] = Am[0];  // first diagonal A[0,0]
                for (npy_intp i = 0; i < n - 1; i++) {
                    sd[i] = Am[(i + 1) * n + i];        // A[i, i+1]
                    diag_aw[i + 1] = Am[(i + 1) * n + (i + 1)];  // A[i+1, i+1]
                }
            }
        }

        pick_pade_structure_c(Am, n, &m, &s);
        if (m < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }
        pade_UV_calc_c(Am, n, m, info);
        if (*info < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }

        SCIPY_C* temp1 = Am;
        SCIPY_C* temp2 = Am1;
        if (s > 0) {
            // Squaring needed.
            if ((lband == 0) || (uband == 0)) {
                // Triangular case - use Fragment 2.1 of Al-Mohy and Higham (2009)
                for (int iter = s - 1; iter >= 0; iter--) {
#if defined(_MSC_VER)
                    SCIPY_C c_one = CPLX_C(1.0f, 0.0f);
                    SCIPY_C c_zero = CPLX_C(0.0f, 0.0f);
                    cgemm_("N", "N", &n, &n, &n, &c_one, temp1, &n, temp1, &n, &c_zero, temp2, &n);
#else
                    cgemm_("N", "N", &n, &n, &n, &(SCIPY_C){CPLX_C(1.0f, 0.0f)}, temp1, &n, temp1, &n, &(SCIPY_C){CPLX_C(0.0f, 0.0f)}, temp2, &n);
#endif
                    SCIPY_C* swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;

                    float scale = exp2f(-(float)iter);  // 2^(-iter)

                    if (is_lower) {

                        for (npy_intp i = 0; i < n - 1; i++) {
#if defined(_MSC_VER)
                            SCIPY_C d_i = _FCmulcr(diag_aw[i], scale);
                            SCIPY_C d_ip1 = _FCmulcr(diag_aw[i + 1], scale);
                            SCIPY_C diff = CPLX_C(crealf(d_ip1) - crealf(d_i), cimagf(d_ip1) - cimagf(d_i));
                            SCIPY_C exp_d_i = cexpf(d_i);
                            SCIPY_C exp_d_ip1 = cexpf(d_ip1);
                            SCIPY_C exp_diff = CPLX_C(crealf(exp_d_ip1) - crealf(exp_d_i), cimagf(exp_d_ip1) - cimagf(exp_d_i));
                            float diff_abs2 = crealf(diff)*crealf(diff) + cimagf(diff)*cimagf(diff);
                            SCIPY_C diff_inv = CPLX_C(crealf(diff)/diff_abs2, -cimagf(diff)/diff_abs2);
                            SCIPY_C exp_sinch = (cabsf(diff) == 0.0f) ? exp_d_i : _FCmulcc(exp_diff, diff_inv);
                            temp1[i * n + i] = exp_d_i;                                            // diagonal
                            temp1[i * n + (i + 1)] = _FCmulcr(_FCmulcc(exp_sinch, sd[i]), scale);  // subdiagonal
#else
                            SCIPY_C d_i = diag_aw[i] * scale;
                            SCIPY_C d_ip1 = diag_aw[i + 1] * scale;
                            SCIPY_C diff = d_ip1 - d_i;
                            SCIPY_C exp_sinch = (cabsf(diff) == 0.0f) ? cexpf(d_i) : (cexpf(d_ip1) - cexpf(d_i)) / diff;
                            temp1[i * n + i] = cexpf(d_i);                      // diagonal
                            temp1[i * n + (i + 1)] = exp_sinch * sd[i] * scale; // subdiagonal
#endif
                        }
#if defined(_MSC_VER)
                        temp1[(n-1) * n + (n-1)] = cexpf(_FCmulcr(diag_aw[n-1], scale)); // last diagonal
#else
                        temp1[(n-1) * n + (n-1)] = cexpf(diag_aw[n-1] * scale); // last diagonal
#endif

                    } else {

#if defined(_MSC_VER)
                        temp1[0] = cexpf(_FCmulcr(diag_aw[0], scale));          // first diagonal
#else
                        temp1[0] = cexpf(diag_aw[0] * scale);                   // first diagonal
#endif
                        for (npy_intp i = 0; i < n - 1; i++) {
#if defined(_MSC_VER)
                            SCIPY_C d_i = _FCmulcr(diag_aw[i], scale);
                            SCIPY_C d_ip1 = _FCmulcr(diag_aw[i + 1], scale);
                            SCIPY_C diff = CPLX_C(crealf(d_ip1) - crealf(d_i), cimagf(d_ip1) - cimagf(d_i));
                            SCIPY_C exp_d_i = cexpf(d_i);
                            SCIPY_C exp_d_ip1 = cexpf(d_ip1);
                            SCIPY_C exp_diff = CPLX_C(crealf(exp_d_ip1) - crealf(exp_d_i), cimagf(exp_d_ip1) - cimagf(exp_d_i));
                            float diff_abs2 = crealf(diff)*crealf(diff) + cimagf(diff)*cimagf(diff);
                            SCIPY_C diff_inv = CPLX_C(crealf(diff)/diff_abs2, -cimagf(diff)/diff_abs2);
                            SCIPY_C exp_sinch = (cabsf(diff) == 0.0f) ? exp_d_i : _FCmulcc(exp_diff, diff_inv);
                            temp1[(i + 1) * n + i] = _FCmulcr(_FCmulcc(exp_sinch, sd[i]), scale); // superdiagonal
                            temp1[(i + 1) * n + (i + 1)] = exp_d_ip1;                             // diagonal
#else
                            SCIPY_C d_i = diag_aw[i] * scale;
                            SCIPY_C d_ip1 = diag_aw[i + 1] * scale;
                            SCIPY_C diff = d_ip1 - d_i;
                            SCIPY_C exp_sinch = (cabsf(diff) == 0.0f) ? cexpf(d_i) : (cexpf(d_ip1) - cexpf(d_i)) / diff;
                            temp1[(i + 1) * n + i] = exp_sinch * sd[i] * scale; // superdiagonal
                            temp1[(i + 1) * n + (i + 1)] = cexpf(d_ip1);        // diagonal
#endif
                        }

                    }
                }
            } else {
                // General dense case, compute A**(2**s) by repeated squaring.
                for (int i = 0; i < s; i++) {
#if defined(_MSC_VER)
                    SCIPY_C c_one = CPLX_C(1.0f, 0.0f);
                    SCIPY_C c_zero = CPLX_C(0.0f, 0.0f);
                    cgemm_("N", "N", &n, &n, &n, &c_one, temp1, &n, temp1, &n, &c_zero, temp2, &n);
#else
                    cgemm_("N", "N", &n, &n, &n, &(SCIPY_C){CPLX_C(1.0f, 0.0f)}, temp1, &n, temp1, &n, &(SCIPY_C){CPLX_C(0.0f, 0.0f)}, temp2, &n);
#endif
                    // Swap pointers
                    SCIPY_C *swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;
                }
            }
        }

        // Copy back the result to output array
        swap_cf_c(temp1, &result[idx*n*n], n, n, n);
    }

    free(Am);
}


void
matrix_exponential_z(PyArrayObject* a, SCIPY_Z* restrict result, int* info)
{
    int m = 0, s = 0, is_lower = 0;
    npy_intp lband = 0, uband = 0;
    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    SCIPY_Z* restrict a_data = (SCIPY_Z*)PyArray_DATA(a);
    int ndim = PyArray_NDIM(a);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(a);      // Array shape
    int n = (int)shape[ndim - 1];                // Slice size
    npy_intp* restrict strides = PyArray_STRIDES(a);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    SCIPY_Z* restrict Am = malloc(sizeof(SCIPY_Z)*(n*n*5 + 2*n));
    if (Am == NULL) { *info = -100; return; }
    SCIPY_Z* restrict Am1 = &Am[n*n];
    // These two arrays are only used for the triangular case for scaling/squaring
    SCIPY_Z* restrict diag_aw = &Am[5*n*n];
    SCIPY_Z* restrict sd = &Am[5*n*n + n];

    /*====================================================================
    |                    MAIN nxn SLICE LOOP                             |
    ====================================================================*/
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        // See sqrtm for explanation of the loop and offset calculations
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        SCIPY_Z* restrict slice_ptr = (SCIPY_Z*)(a_data + (offset/sizeof(SCIPY_Z)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Am1[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(SCIPY_Z)) + (j*strides[ndim - 1]/sizeof(SCIPY_Z)));
            }
        }
        bandwidth_z(Am1, n, n, &lband, &uband);
        if ((lband == 0) && (uband == 0)) {
            for (npy_intp i = 0; i < n; i++) {
                result[idx*n*n + i*n + i] = cexp(Am1[i*n + i]);
            }
            continue;
        }
        swap_cf_z(Am1, Am, n, n, n);


        is_lower = (uband == 0);
        if ((lband == 0) || (uband == 0)) {
            if (is_lower) {
                // Column-major: iterate columns, grab diagonal and the subdiagonal elements
                for (npy_intp i = 0; i < n - 1; i++) {
                    diag_aw[i] = Am[i * n + i];         // A[i, i]
                    sd[i] = Am[i * n + (i + 1)];        // A[i+1, i]
                }
                diag_aw[n - 1] = Am[(n - 1) * n + (n - 1)];  // last diagonal
            } else {
                // Column-major: iterate columns, grab superdiagonal and diagonal elements
                diag_aw[0] = Am[0];  // first diagonal A[0,0]
                for (npy_intp i = 0; i < n - 1; i++) {
                    sd[i] = Am[(i + 1) * n + i];        // A[i, i+1]
                    diag_aw[i + 1] = Am[(i + 1) * n + (i + 1)];  // A[i+1, i+1]
                }
            }
        }

        pick_pade_structure_z(Am, n, &m, &s);
        if (m < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }
        pade_UV_calc_z(Am, n, m, info);
        if (*info < 0) {
            free(Am);
            *info = -100 + *info;
            return;
        }

        SCIPY_Z* temp1 = Am;
        SCIPY_Z* temp2 = Am1;
        if (s > 0) {
            // Squaring needed.
            if ((lband == 0) || (uband == 0)) {
                // Triangular case - use Fragment 2.1 of Al-Mohy and Higham (2009)
                for (int iter = s - 1; iter >= 0; iter--) {
#if defined(_MSC_VER)
                    SCIPY_Z z_one = CPLX_Z(1.0, 0.0);
                    SCIPY_Z z_zero = CPLX_Z(0.0, 0.0);
                    zgemm_("N", "N", &n, &n, &n, &z_one, temp1, &n, temp1, &n, &z_zero, temp2, &n);
#else
                    zgemm_("N", "N", &n, &n, &n, &(SCIPY_Z){CPLX_Z(1.0, 0.0)}, temp1, &n, temp1, &n, &(SCIPY_Z){CPLX_Z(0.0, 0.0)}, temp2, &n);
#endif
                    SCIPY_Z* swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;

                    double scale = exp2(-(double)iter);  // 2^(-iter)

                    if (is_lower) {

                        for (npy_intp i = 0; i < n - 1; i++) {
#if defined(_MSC_VER)
                            SCIPY_Z d_i = _Cmulcr(diag_aw[i], scale);
                            SCIPY_Z d_ip1 = _Cmulcr(diag_aw[i + 1], scale);
                            SCIPY_Z diff = CPLX_Z(creal(d_ip1) - creal(d_i), cimag(d_ip1) - cimag(d_i));
                            SCIPY_Z exp_d_i = cexp(d_i);
                            SCIPY_Z exp_d_ip1 = cexp(d_ip1);
                            SCIPY_Z exp_diff = CPLX_Z(creal(exp_d_ip1) - creal(exp_d_i), cimag(exp_d_ip1) - cimag(exp_d_i));
                            double diff_abs2 = creal(diff)*creal(diff) + cimag(diff)*cimag(diff);
                            SCIPY_Z diff_inv = CPLX_Z(creal(diff)/diff_abs2, -cimag(diff)/diff_abs2);
                            SCIPY_Z exp_sinch = (cabs(diff) == 0.0) ? exp_d_i : _Cmulcc(exp_diff, diff_inv);
                            temp1[i * n + i] = exp_d_i;                                          // diagonal
                            temp1[i * n + (i + 1)] = _Cmulcr(_Cmulcc(exp_sinch, sd[i]), scale); // subdiagonal
#else
                            SCIPY_Z d_i = diag_aw[i] * scale;
                            SCIPY_Z d_ip1 = diag_aw[i + 1] * scale;
                            SCIPY_Z diff = d_ip1 - d_i;
                            SCIPY_Z exp_sinch = (cabs(diff) == 0.0) ? cexp(d_i) : (cexp(d_ip1) - cexp(d_i)) / diff;
                            temp1[i * n + i] = cexp(d_i);                       // diagonal
                            temp1[i * n + (i + 1)] = exp_sinch * sd[i] * scale; // subdiagonal
#endif
                        }
#if defined(_MSC_VER)
                        temp1[(n-1) * n + (n-1)] = cexp(_Cmulcr(diag_aw[n-1], scale)); // last diagonal
#else
                        temp1[(n-1) * n + (n-1)] = cexp(diag_aw[n-1] * scale); // last diagonal
#endif

                    } else {

#if defined(_MSC_VER)
                        temp1[0] = cexp(_Cmulcr(diag_aw[0], scale));           // first diagonal
#else
                        temp1[0] = cexp(diag_aw[0] * scale);                    // first diagonal
#endif
                        for (npy_intp i = 0; i < n - 1; i++) {
#if defined(_MSC_VER)
                            SCIPY_Z d_i = _Cmulcr(diag_aw[i], scale);
                            SCIPY_Z d_ip1 = _Cmulcr(diag_aw[i + 1], scale);
                            SCIPY_Z diff = CPLX_Z(creal(d_ip1) - creal(d_i), cimag(d_ip1) - cimag(d_i));
                            SCIPY_Z exp_d_i = cexp(d_i);
                            SCIPY_Z exp_d_ip1 = cexp(d_ip1);
                            SCIPY_Z exp_diff = CPLX_Z(creal(exp_d_ip1) - creal(exp_d_i), cimag(exp_d_ip1) - cimag(exp_d_i));
                            double diff_abs2 = creal(diff)*creal(diff) + cimag(diff)*cimag(diff);
                            SCIPY_Z diff_inv = CPLX_Z(creal(diff)/diff_abs2, -cimag(diff)/diff_abs2);
                            SCIPY_Z exp_sinch = (cabs(diff) == 0.0) ? exp_d_i : _Cmulcc(exp_diff, diff_inv);
                            temp1[(i + 1) * n + i] = _Cmulcr(_Cmulcc(exp_sinch, sd[i]), scale); // superdiagonal
                            temp1[(i + 1) * n + (i + 1)] = exp_d_ip1;                           // diagonal
#else
                            SCIPY_Z d_i = diag_aw[i] * scale;
                            SCIPY_Z d_ip1 = diag_aw[i + 1] * scale;
                            SCIPY_Z diff = d_ip1 - d_i;
                            SCIPY_Z exp_sinch = (cabs(diff) == 0.0) ? cexp(d_i) : (cexp(d_ip1) - cexp(d_i)) / diff;
                            temp1[(i + 1) * n + i] = exp_sinch * sd[i] * scale; // superdiagonal
                            temp1[(i + 1) * n + (i + 1)] = cexp(d_ip1);         // diagonal
#endif
                        }

                    }
                }
            } else {
                // General dense case, compute A**(2**s) by repeated squaring.
                for (int i = 0; i < s; i++) {
#if defined(_MSC_VER)
                    SCIPY_Z z_one = CPLX_Z(1.0, 0.0);
                    SCIPY_Z z_zero = CPLX_Z(0.0, 0.0);
                    zgemm_("N", "N", &n, &n, &n, &z_one, temp1, &n, temp1, &n, &z_zero, temp2, &n);
#else
                    zgemm_("N", "N", &n, &n, &n, &(SCIPY_Z){CPLX_Z(1.0, 0.0)}, temp1, &n, temp1, &n, &(SCIPY_Z){CPLX_Z(0.0, 0.0)}, temp2, &n);
#endif
                    // Swap pointers
                    SCIPY_Z *swap = temp1;
                    temp1 = temp2;
                    temp2 = swap;
                }
            }
        }

        // Copy back the result to output array
        swap_cf_z(temp1, &result[idx*n*n], n, n, n);
    }

    free(Am);
}
