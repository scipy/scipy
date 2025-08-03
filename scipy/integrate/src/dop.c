#include "dop.h"

static void
dp86co(const int n, dopri_fcn* fcn, double* x, double* y, double *xend, double hmax,
       double h, double* rtol, double* atol, const int itol, dopri_solout* solout, const int iout,
       int* ierr, const int nmax, const double uround, const int meth, const int nstiff,
       const double safe, const double beta, const double fac1, const double fac2,
       double* work, int* iwork, int* nrd, double* rpar, int* ipar, int* nfcn,
       int* nstep, int* naccept, int* nreject);

static double
hinit853(int n, dopri_fcn* fcn, double* x, double* y, double posneg,
         double* f0, double* f1, double* y1, int iord, double hmax, double* atol,
         double* rtol, int itol, double* rpar, int* ipar);

static void
dopcor(const int n, dopri_fcn* fcn, double* x, double* y, double *xend,
       double hmax, double h, double* rtol, double* atol, const int itol,
       dopri_solout* solout, const int iout, int* ierr, const int nmax,
       const double uround, const int meth, int nstiff, const double safe,
       const double beta, const double fac1, const double fac2, double* work,
       int* iwork, int* nrd, double* rpar, int* ipar, int* nfcn, int* nstep,
       int* naccept, int* nreject);

static double
hinit(int n, dopri_fcn* fcn, double* x, double* y, double posneg, double* f0,
      double* f1, double* y1, int iord, double hmax, double* atol, double* rtol,
      int itol, double* rpar, int* ipar);

void
dopri853(const int n, dopri_fcn* fcn, double* x, double* y, double *xend, double* rtol,
         double* atol, const int itol, dopri_solout* solout, const int iout, double* work,
         int* iwork, double* rpar, int* ipar, int* ierr)
{
    int i, meth, nfcn, nmax, nrdens, nstep, nstiff, naccept, nreject, arret;
    double beta, h, hmax, fac1, fac2, safe, uround;
    nfcn = 0;
    nstep = 0;
    naccept = 0;
    nreject = 0;
    arret = 0;
    *ierr = 0;

    // -------- nmax, the maximal number of steps -----
    if (iwork[0] == 0)
    {
        nmax = 100000;
    } else {
        nmax = iwork[0];
        if (nmax <= 0)
        {
            arret = 11;
        }
    }

    // -------- meth, coefficients of the method
    if (iwork[1] == 0)
    {
        meth = 1;
    } else {
        meth = iwork[1];
        if ((meth <= 0) || (meth >= 4))
        {
            arret = 12;
        }
    }
    // -------- nstiff, parameter for stiffness detection
    nstiff = iwork[3];
    if (nstiff == 0) { nstiff = 1000; }
    if (nstiff < 0) { nstiff = nmax + 10; }

    // -------- nrdens, number of dense output components
    nrdens = iwork[4];
    if ((nrdens < 0) || (nrdens > n))
    {
        arret = 13;
    } else {
        if (nrdens == n) { for (i = 0; i < nrdens; i++) { iwork[i + 20] = i; } }
    }

    // -------- uround, smallest number satisfying 1.d0+uround>1.d0
    if (work[0] == 0.0)
    {
        uround = 2.220446049250313e-16;
    } else {
        uround = work[0];
        if ((uround <= 1e-35) || (uround >= 1.0))
        {
            arret = 14;
        }
    }

    // -------  safety factor -------------
    if (work[1] == 0.0)
    {
        safe = 0.9;
    } else {
        safe = work[1];
        if ((safe <= 1e-4) || (safe >= 1.0))
        {
            arret = 15;
        }
    }

    // -------  fac1,fac2, parameters for step size selection
    fac1 = (work[2] == 0.0 ? 1.0 / 3.0 : work[2]);
    fac2 = (work[3] == 0.0 ? 6.0 : work[3]);

    // --------- beta for step control stabilization -----------
    if (work[4] <= 0.0)
    {
        beta = 0.0;
    } else {
        beta = work[4];
        if (beta > 0.2)
        {
            arret = 16;
        }
    }

    // -------- maximal step size
    if (work[5] == 0.0)
    {
        hmax = *xend - *x;
    } else {
        hmax = work[5];
    }

    // -------- initial step size
    h = work[6];

    // -------- when a fail has occurred, we return with wrong input error code
    if (arret > 0) { *ierr = -arret; return; }

    // -------- call to core integrator ------------
    dp86co(n, fcn, x, y, xend, hmax, h, rtol, atol, itol, solout, iout, ierr, nmax,
           uround, meth, nstiff, safe, beta, fac1, fac2, work, iwork, &nrdens, rpar,
           ipar, &nfcn, &nstep, &naccept, &nreject);

    work[6] = h;
    iwork[16] = nfcn;
    iwork[17] = nstep;
    iwork[18] = naccept;
    iwork[19] = nreject;

    return;
}

static void
dp86co(const int n, dopri_fcn* fcn, double *x, double* y, double *xend, double hmax,
       double h, double* rtol, double* atol, const int itol, dopri_solout* solout, const int iout,
       int* ierr, const int nmax, const double uround, const int meth, const int nstiff,
       const double safe, const double beta, const double fac1, const double fac2,
       double* work, int* iwork, int* nrd, double* rpar, int* ipar, int* nfcn,
       int* nstep, int* naccept, int* nreject)
{

    double atoli, rtoli, bspl, facold, expo1, facc1, facc2, posneg, hlamb, hnew, fac;
    double err, err2, fac11, deno, hout, xold, xph, ydiff, sk, erri, stnum, stden;
    int irtrn, iasti, iord, i, j, last = 0, nonsti = 0, reject = 0;
    (void)meth;
    // Pointers into work
    double* restrict k1   = &work[20];
    double* restrict k2   = &work[20 + n];
    double* restrict k3   = &work[20 + n*2];
    double* restrict k4   = &work[20 + n*3];
    double* restrict k5   = &work[20 + n*4];
    double* restrict k6   = &work[20 + n*5];
    double* restrict k7   = &work[20 + n*6];
    double* restrict k8   = &work[20 + n*7];
    double* restrict k9   = &work[20 + n*8];
    double* restrict k10  = &work[20 + n*9];
    double* restrict y1   = &work[20 + n*10];
    double* restrict cont = &work[20 + n*11];
    int* icomp   = &iwork[20];

    // Integration parameters
    const double c2    =  0.526001519587677318785587544488e-01;
    const double c3    =  0.789002279381515978178381316732e-01;
    const double c4    =  0.118350341907227396726757197510e+00;
    const double c5    =  0.281649658092772603273242802490e+00;
    const double c6    =  0.333333333333333333333333333333e+00;
    const double c7    =  0.25e+00;
    const double c8    =  0.307692307692307692307692307692e+00;
    const double c9    =  0.651282051282051282051282051282e+00;
    const double c10   =  0.6e+00;
    const double c11   =  0.857142857142857142857142857142e+00;
    const double c14   =  0.1e+00;
    const double c15   =  0.2e+00;
    const double c16   =  0.777777777777777777777777777778e+00;
    const double b1    =  5.42937341165687622380535766363e-2;
    const double b6    =  4.45031289275240888144113950566e0;
    const double b7    =  1.89151789931450038304281599044e0;
    const double b8    = -5.8012039600105847814672114227e0;
    const double b9    =  3.1116436695781989440891606237e-1;
    const double b10   = -1.52160949662516078556178806805e-1;
    const double b11   =  2.01365400804030348374776537501e-1;
    const double b12   =  4.47106157277725905176885569043e-2;
    const double bhh1  =  0.244094488188976377952755905512e+00;
    const double bhh2  =  0.733846688281611857341361741547e+00;
    const double bhh3  =  0.220588235294117647058823529412e-01;
    const double er1   =  0.1312004499419488073250102996e-01;
    const double er6   = -0.1225156446376204440720569753e+01;
    const double er7   = -0.4957589496572501915214079952e+00;
    const double er8   =  0.1664377182454986536961530415e+01;
    const double er9   = -0.3503288487499736816886487290e+00;
    const double er10  =  0.3341791187130174790297318841e+00;
    const double er11  =  0.8192320648511571246570742613e-01;
    const double er12  = -0.2235530786388629525884427845e-01;
    const double a21   =  5.26001519587677318785587544488e-2;
    const double a31   =  1.97250569845378994544595329183e-2;
    const double a32   =  5.91751709536136983633785987549e-2;
    const double a41   =  2.95875854768068491816892993775e-2;
    const double a43   =  8.87627564304205475450678981324e-2;
    const double a51   =  2.41365134159266685502369798665e-1;
    const double a53   = -8.84549479328286085344864962717e-1;
    const double a54   =  9.24834003261792003115737966543e-1;
    const double a61   =  3.7037037037037037037037037037e-2;
    const double a64   =  1.70828608729473871279604482173e-1;
    const double a65   =  1.25467687566822425016691814123e-1;
    const double a71   =  3.7109375e-2;
    const double a74   =  1.70252211019544039314978060272e-1;
    const double a75   =  6.02165389804559606850219397283e-2;
    const double a76   = -1.7578125e-2;
    const double a81   =  3.70920001185047927108779319836e-2;
    const double a84   =  1.70383925712239993810214054705e-1;
    const double a85   =  1.07262030446373284651809199168e-1;
    const double a86   = -1.53194377486244017527936158236e-2;
    const double a87   =  8.27378916381402288758473766002e-3;
    const double a91   =  6.24110958716075717114429577812e-1;
    const double a94   = -3.36089262944694129406857109825e0;
    const double a95   = -8.68219346841726006818189891453e-1;
    const double a96   =  2.75920996994467083049415600797e1;
    const double a97   =  2.01540675504778934086186788979e1;
    const double a98   = -4.34898841810699588477366255144e1;
    const double a101  =  4.77662536438264365890433908527e-1;
    const double a104  = -2.48811461997166764192642586468e0;
    const double a105  = -5.90290826836842996371446475743e-1;
    const double a106  =  2.12300514481811942347288949897e1;
    const double a107  =  1.52792336328824235832596922938e1;
    const double a108  = -3.32882109689848629194453265587e1;
    const double a109  = -2.03312017085086261358222928593e-2;
    const double a111  = -9.3714243008598732571704021658e-1;
    const double a114  =  5.18637242884406370830023853209e0;
    const double a115  =  1.09143734899672957818500254654e0;
    const double a116  = -8.14978701074692612513997267357e0;
    const double a117  = -1.85200656599969598641566180701e1;
    const double a118  =  2.27394870993505042818970056734e1;
    const double a119  =  2.49360555267965238987089396762e0;
    const double a1110 = -3.0467644718982195003823669022e0;
    const double a121  =  2.27331014751653820792359768449e0;
    const double a124  = -1.05344954667372501984066689879e1;
    const double a125  = -2.00087205822486249909675718444e0;
    const double a126  = -1.79589318631187989172765950534e1;
    const double a127  =  2.79488845294199600508499808837e1;
    const double a128  = -2.85899827713502369474065508674e0;
    const double a129  = -8.87285693353062954433549289258e0;
    const double a1210 =  1.23605671757943030647266201528e1;
    const double a1211 =  6.43392746015763530355970484046e-1;
    const double a141  =  5.61675022830479523392909219681e-2;
    const double a147  =  2.53500210216624811088794765333e-1;
    const double a148  = -2.46239037470802489917441475441e-1;
    const double a149  = -1.24191423263816360469010140626e-1;
    const double a1410 =  1.5329179827876569731206322685e-1;
    const double a1411 =  8.20105229563468988491666602057e-3;
    const double a1412 =  7.56789766054569976138603589584e-3;
    const double a1413 = -8.298e-3;
    const double a151  =  3.18346481635021405060768473261e-2;
    const double a156  =  2.83009096723667755288322961402e-2;
    const double a157  =  5.35419883074385676223797384372e-2;
    const double a158  = -5.49237485713909884646569340306e-2;
    const double a1511 = -1.08347328697249322858509316994e-4;
    const double a1512 =  3.82571090835658412954920192323e-4;
    const double a1513 = -3.40465008687404560802977114492e-4;
    const double a1514 =  1.41312443674632500278074618366e-1;
    const double a161  = -4.28896301583791923408573538692e-1;
    const double a166  = -4.69762141536116384314449447206e0;
    const double a167  =  7.68342119606259904184240953878e0;
    const double a168  =  4.06898981839711007970213554331e0;
    const double a169  =  3.56727187455281109270669543021e-1;
    const double a1613 = -1.39902416515901462129418009734e-3;
    const double a1614 =  2.9475147891527723389556272149e0;
    const double a1615 = -9.15095847217987001081870187138e0;
    const double d41   = -0.84289382761090128651353491142e+01;
    const double d46   =  0.56671495351937776962531783590e+00;
    const double d47   = -0.30689499459498916912797304727e+01;
    const double d48   =  0.23846676565120698287728149680e+01;
    const double d49   =  0.21170345824450282767155149946e+01;
    const double d410  = -0.87139158377797299206789907490e+00;
    const double d411  =  0.22404374302607882758541771650e+01;
    const double d412  =  0.63157877876946881815570249290e+00;
    const double d413  = -0.88990336451333310820698117400e-01;
    const double d414  =  0.18148505520854727256656404962e+02;
    const double d415  = -0.91946323924783554000451984436e+01;
    const double d416  = -0.44360363875948939664310572000e+01;
    const double d51   =  0.10427508642579134603413151009e+02;
    const double d56   =  0.24228349177525818288430175319e+03;
    const double d57   =  0.16520045171727028198505394887e+03;
    const double d58   = -0.37454675472269020279518312152e+03;
    const double d59   = -0.22113666853125306036270938578e+02;
    const double d510  =  0.77334326684722638389603898808e+01;
    const double d511  = -0.30674084731089398182061213626e+02;
    const double d512  = -0.93321305264302278729567221706e+01;
    const double d513  =  0.15697238121770843886131091075e+02;
    const double d514  = -0.31139403219565177677282850411e+02;
    const double d515  = -0.93529243588444783865713862664e+01;
    const double d516  =  0.35816841486394083752465898540e+02;
    const double d61   =  0.19985053242002433820987653617e+02;
    const double d66   = -0.38703730874935176555105901742e+03;
    const double d67   = -0.18917813819516756882830838328e+03;
    const double d68   =  0.52780815920542364900561016686e+03;
    const double d69   = -0.11573902539959630126141871134e+02;
    const double d610  =  0.68812326946963000169666922661e+01;
    const double d611  = -0.10006050966910838403183860980e+01;
    const double d612  =  0.77771377980534432092869265740e+00;
    const double d613  = -0.27782057523535084065932004339e+01;
    const double d614  = -0.60196695231264120758267380846e+02;
    const double d615  =  0.84320405506677161018159903784e+02;
    const double d616  =  0.11992291136182789328035130030e+02;
    const double d71   = -0.25693933462703749003312586129e+02;
    const double d76   = -0.15418974869023643374053993627e+03;
    const double d77   = -0.23152937917604549567536039109e+03;
    const double d78   =  0.35763911791061412378285349910e+03;
    const double d79   =  0.93405324183624310003907691704e+02;
    const double d710  = -0.37458323136451633156875139351e+02;
    const double d711  =  0.10409964950896230045147246184e+03;
    const double d712  =  0.29840293426660503123344363579e+02;
    const double d713  = -0.43533456590011143754432175058e+02;
    const double d714  =  0.96324553959188282948394950600e+02;
    const double d715  = -0.39177261675615439165231486172e+02;
    const double d716  = -0.14972683625798562581422125276e+03;

    // Initializations
    facold = 1.0e-4;
    expo1 = 1.25e-1 - beta * 0.2;
    facc1 = 1.0 / fac1;
    facc2 = 1.0 / fac2;
    posneg = copysign(1.0, *xend - *x);
    // Initial preparations
    atoli = atol[0];
    rtoli = rtol[0];
    hlamb = 0.0;
    iasti = 0;

    fcn(n, *x, y, k1, rpar, ipar);
    hmax = fabs(hmax);
    iord = 8;
    if (h == 0.0) {
        h = hinit853(n, fcn, x, y, posneg, k1, k2, k3, iord, hmax, atol, rtol, itol, rpar, ipar);
    }

    *nfcn += 2;
    xold = *x;

    if (iout >= 1) {
        irtrn = 1;
        hout = 1.0;
        solout(*naccept + 1, xold, *x, y, n, cont, icomp, *nrd, rpar, ipar, &irtrn);
        if (irtrn < 0) { *ierr = -3; return; }
    }

    // basic integration step
    while (1) {
        if (*nstep > nmax) { *ierr = -2; return; }
        if (0.1 * fabs(h) <= fabs(*x) * uround) { *ierr = -3; return; }
        if ((*x + (1.01 * h) - *xend) * posneg > 0.0) {
            h = *xend - *x;
            last = 1;
        }

        *nstep += 1;

        // the twelve stages
        if (irtrn >= 2) { fcn(n, *x, y, k1, rpar, ipar); }

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*a21*k1[i]; }  // 22
        fcn(n, *x + c2*h, y1, k2, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a31*k1[i] + a32*k2[i]); }  // 23
        fcn(n, *x + c3*h, y1, k3, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a41*k1[i] + a43*k3[i]); }  // 24
        fcn(n, *x + c4*h, y1, k4, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a51*k1[i] + a53*k3[i] + a54*k4[i]); }  // 25
        fcn(n, *x + c5*h, y1, k5, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a61*k1[i] + a64*k4[i] + a65*k5[i]); }  // 26
        fcn(n, *x + c6*h, y1, k6, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a71*k1[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]); }  // 27
        fcn(n, *x + c7*h, y1, k7, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a81*k1[i] + a84*k4[i] + a85*k5[i] + a86*k6[i] + a87*k7[i]); }  // 28
        fcn(n, *x + c8*h, y1, k8, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a91*k1[i] + a94*k4[i] + a95*k5[i] + a96*k6[i] + a97*k7[i] + a98*k8[i]); }  // 29
        fcn(n, *x + c9*h, y1, k9, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a101*k1[i] + a104*k4[i] + a105*k5[i] + a106*k6[i] + a107*k7[i] + a108*k8[i] + a109*k9[i]); }  // 30
        fcn(n, *x + c10*h, y1, k10, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a111*k1[i] + a114*k4[i] + a115*k5[i] + a116*k6[i] + a117*k7[i] + a118*k8[i] + a119*k9[i] + a1110*k10[i]); }  // 31
        fcn(n, *x + c11*h, y1, k2, rpar, ipar);

        xph = *x + h;
        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(a121*k1[i] + a124*k4[i] + a125*k5[i] + a126*k6[i] + a127*k7[i] + a128*k8[i] + a129*k9[i] + a1210*k10[i] + a1211*k2[i]); }  // 32
        fcn(n, xph, y1, k3, rpar, ipar);

        *nfcn += 11;

        for (i = 0; i < n; i++)
        {
            k4[i] = b1*k1[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i] + b10*k10[i] + b11*k2[i] + b12*k3[i];
            k5[i] = y[i] + h*k4[i];
        }

        // error estimation
        err = 0.0;
        err2 = 0.0;

        if (itol == 0)
        {
            for (i = 0; i < n; i++)
            {
                sk = atoli + rtoli * fmax(fabs(y[i]), fabs(k5[i]));
                erri = k4[i] - bhh1 * k1[i] - bhh2 * k9[i] - bhh3 * k3[i];
                err2 = err2 + pow(erri / sk, 2);
                erri = er1 * k1[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] + er10 * k10[i] + er11 * k2[i] + er12 * k3[i];
                err = err + pow(erri / sk, 2.0);
            }
            // 41
        } else {
            // atol, rtol provided
            for (i = 0; i < n; i++)
            {
                sk = atol[i] + rtol[i] * fmax(fabs(y[i]), fabs(k5[i]));
                erri = k4[i] - bhh1 * k1[i] - bhh2 * k9[i] - bhh3 * k3[i];
                err2 = err2 + pow(erri / sk, 2.0);
                erri = er1 * k1[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] + er10 * k10[i] + er11 * k2[i] + er12 * k3[i];
                err = err + pow(erri / sk, 2.0);
            }
            // 42

        }

        deno = err + 0.01*err2;
        if (deno <= 0.0) { deno = 1.0; }

        err = fabs(h)*err*sqrt(1.0 / (n*deno));

        // --- computation of hnew
        fac11 = pow(err, expo1);
        // --- lund-stabilization
        fac = fac11 / pow(facold, beta);
        // --- we require fac1 <= hnew/h <= fac2
        fac = fmax(facc2, fmin(facc1,fac/safe));
        hnew = h/fac;

        if (err <= 1.0) {
            // --- step is accepted
            facold = fmax(err, 1.0e-4);
            *naccept += 1;
            fcn(n, xph, k5, k4, rpar, ipar);
            *nfcn += 1;

            // ------- stiffness detection
            if ((*naccept % nstiff == 0) || (iasti  > 0))
            {
                stnum = 0.0;
                stden = 0.0;
                for (i = 0; i < n; i++)
                {
                    stnum = stnum + pow(k4[i] - k3[i], 2.0);
                    stden = stden + pow(k5[i] - y1[i], 2.0);
                }
                // 64
                if (stden > 0.0) { hlamb = fabs(h) * sqrt(stnum / stden); }
                if (hlamb > 6.1)
                {
                    nonsti = 0;
                    iasti += 1;
                    if (iasti == 15) { *ierr = -4; return; }
                } else {
                    nonsti += 1;
                    if (nonsti == 6) { iasti = 0; }
                }
            }
            // ------- final preparation for dense output
            if (iout >= 2)
            {
                // save the first function evaluations
                for (j = 0; j < (*nrd); j++)
                {
                    i = icomp[j];
                    cont[j] = y[i];
                    ydiff = k5[i] - y[i];
                    bspl = h*k1[i] - ydiff;
                    cont[j + (*nrd)]   = ydiff;
                    cont[j + (*nrd)*2] = bspl;
                    cont[j + (*nrd)*3] = ydiff - h*k4[i] - bspl;
                    cont[j + (*nrd)*4] = d41*k1[i] + d46*k6[i] + d47*k7[i] + d48*k8[i] + d49*k9[i] + d410*k10[i] + d411*k2[i] + d412*k3[i];
                    cont[j + (*nrd)*5] = d51*k1[i] + d56*k6[i] + d57*k7[i] + d58*k8[i] + d59*k9[i] + d510*k10[i] + d511*k2[i] + d512*k3[i];
                    cont[j + (*nrd)*6] = d61*k1[i] + d66*k6[i] + d67*k7[i] + d68*k8[i] + d69*k9[i] + d610*k10[i] + d611*k2[i] + d612*k3[i];
                    cont[j + (*nrd)*7] = d71*k1[i] + d76*k6[i] + d77*k7[i] + d78*k8[i] + d79*k9[i] + d710*k10[i] + d711*k2[i] + d712*k3[i];
                }
                // the next three function evaluations
                for (i = 0; i < n; i++)
                {
                    y1[i] = y[i] + h*(a141*k1[i] + a147*k7[i] + a148*k8[i] + a149*k9[i] + a1410*k10[i] + a1411*k2[i] + a1412*k3[i] + a1413*k4[i]);
                }
                fcn(n, *x+c14*h, y1, k10, rpar, ipar);

                for (i = 0;i<n;i++)
                {
                    y1[i] = y[i] + h*(a151*k1[i] + a156*k6[i] + a157*k7[i] + a158*k8[i] + a1511*k2[i]+ a1512*k3[i] + a1513*k4[i] + a1514*k10[i]);
                }
                fcn(n, *x+c15*h, y1, k2, rpar, ipar);
                for (i = 0; i < n; i++)
                {
                    y1[i] = y[i] + h*(a161*k1[i] + a166*k6[i] + a167*k7[i] + a168*k8[i] + a169*k9[i] + a1613*k4[i] + a1614*k10[i] + a1615*k2[i]);
                }
                fcn(n, *x+c16*h, y1, k3, rpar, ipar);
                nfcn += 3;
                // final preparation
                for (j = 0; j < (*nrd); j++)
                {
                    i = icomp[j];
                    cont[j + (*nrd)*4] = h*(cont[j + (*nrd)*4] + d413*k4[i] + d414*k10[i] + d415*k2[i] + d416*k3[i]);
                    cont[j + (*nrd)*5] = h*(cont[j + (*nrd)*5] + d513*k4[i] + d514*k10[i] + d515*k2[i] + d516*k3[i]);
                    cont[j + (*nrd)*6] = h*(cont[j + (*nrd)*6] + d613*k4[i] + d614*k10[i] + d615*k2[i] + d616*k3[i]);
                    cont[j + (*nrd)*7] = h*(cont[j + (*nrd)*7] + d713*k4[i] + d714*k10[i] + d715*k2[i] + d716*k3[i]);
                }
                // 63
                hout = h;
            }
            for (i = 0; i < n; i++)
            {
                k1[i] = k4[i];
                y[i] = k5[i];
            }
            // 67
            xold = *x;
            *x = xph;
            if (iout >= 1)
            {
                solout(*naccept+1, xold, *x, y, n, cont, icomp, *nrd, rpar, ipar, &irtrn);
                if (irtrn < 0) { *ierr = 2; return; }
            }
            // Normal exit
            if (last)
            {
                h = hnew;
                *ierr = 1;
                return;
            }
            if (fabs(hnew) > hmax) { hnew = posneg*hmax; }
            if (reject) { hnew = posneg*fmin(fabs(hnew), fabs(h)); }
            reject = 0;
        } else {
            // step is rejected
            hnew = h / fmin(facc1, facc1 / safe);
            reject = 1;
            if (*naccept >= 1) { *nreject += 1; }
            last = 0;
        }
        h = hnew;
    }
}


static double
hinit853(int n, dopri_fcn* fcn, double *x, double* y, double posneg,
         double* f0, double* f1, double* y1, int iord, double hmax, double* atol,
         double* rtol, int itol, double* rpar, int* ipar)
{
    double dnf = 0.0, dny = 0.0, atoli = atol[0], rtoli = rtol[0];
    double h, sk, der2, der12, h1;
    int i;

    // ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
    // ----   H = 0.01 * NORM (Y0) / NORM (F0)
    // ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
    // ---- COMPARED TO THE SOLUTION
    if (itol == 0)
    {
        for (i = 0; i < n; i++)
        {
            sk = atoli + rtoli * fabs(y[i]);
            dnf = dnf + pow(f0[i] / sk, 2);
            dny = dny + pow(y[i] / sk, 2);
        }
        // 10
    } else {
        for (i = 0; i < n; i++)
        {
            sk = atol[i] + rtol[i] * fabs(y[i]);
            dnf = dnf + pow(f0[i] / sk, 2);
            dny = dny + pow(y[i] / sk, 2);
        }
        // 11
    }

    if ((dnf <= 1.0e-10) || (dny <= 1.0e-10)) {
        h = 1.0e-6;
    } else {
        h = sqrt(dny / dnf) * 0.01;
    }

    h = fmin(h, hmax);
    h = copysign(h, posneg);

    // perform an explicit euler step
    for (i = 0; i < n; i++)
    {
        y1[i] = y[i] + h * f0[i];
    }
    fcn(n, *x + h, y1, f1, rpar, ipar);

    // estimate the second derivative of the solution
    der2 = 0.0;
    if (itol == 0)
    {
        for (i = 0; i < n; i++)
        {
            sk = atoli + rtoli * fabs(y[i]);
            der2 = der2 + pow((f1[i] - f0[i]) / sk, 2.0);
        }
        // 15
    } else {
        for (i = 0; i < n; i++)
        {
            sk = atol[i] + rtol[i] * fabs(y[i]);
            der2 = der2 + pow((f1[i] - f0[i]) / sk, 2.0);
        }
        // 16
    }
    der2 = sqrt(der2) / h;

    // ---- STEP SIZE IS COMPUTED SUCH THAT
    // ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
    der12 = fmax(fabs(der2), sqrt(dnf));
    if (der12 <= 1.0e-15)
    {
        h1 = fmax(1.0e-6, fabs(h) * 1.0e-3);
    } else {
        h1 = pow(0.01 / der12, 1.0 / iord);
    }

    h = fmin(fmin(100.0 * fabs(h), h1), hmax);
    return copysign(h, posneg);
}


static double
contd8(const int ii, double *x, double* con, int* icomp, int nd, const double xold, const double h)
{
    int i = 0;
    double s, s1, conpar, res;
    // compute the index of the ii-th component
    for (int j = 0; j < nd; j++) { if (icomp[j] == ii) { i = j; } }

    if (i == 0) { return -1; }
    // calculate the dense output
    s = (*x - xold) / h;
    s1 = 1.0 - s;
    conpar = con[i + nd*4] + s*(con[i + nd*5] + s1*(con[i + nd*6] + s*con[i + nd*7]));
    res = con[i] + s*(con[i + nd] + s1*(con[i + nd*2] + s*(con[i + nd*3] + s1*conpar)));
    return res;
}

void
dopri5(const int n, dopri_fcn* fcn, double *x, double* y, double *xend, double* rtol,
       double* atol, const int itol, dopri_solout* solout, const int iout,
       double* work, int* iwork, double* rpar, int* ipar, int* ierr)
{
    int i, meth, nfcn, nmax, nrdens, nstep, nstiff, naccept, nreject, arret;
    double beta, h, hmax, fac1, fac2, safe, uround;
    nfcn = 0;
    nstep = 0;
    naccept = 0;
    nreject = 0;
    arret = 0;
    *ierr = 0;

    // -------- nmax, the maximal number of steps -----
    if (iwork[0] == 0)
    {
        nmax = 100000;
    } else {
        nmax = iwork[0];
        if (nmax <= 0)
        {
            arret = 11;
        }
    }

    // -------- meth, coefficients of the method
    if (iwork[1] == 0)
    {
        meth = 1;
    } else {
        meth = iwork[1];
        if ((meth <= 0) || (meth >= 4))
        {
            arret = 12;
        }
    }
    // -------- nstiff, parameter for stiffness detection
    nstiff = iwork[3];
    if (nstiff == 0) { nstiff = 1000; }
    if (nstiff < 0) { nstiff = nmax + 10; }

    // -------- nrdens, number of dense output components
    nrdens = iwork[4];
    if ((nrdens < 0) || (nrdens > n))
    {
        arret = 13;
    } else {
        if (nrdens == n) { for (i = 0; i < nrdens; i++) { iwork[i + 20] = i; } }
    }

    // -------- uround, smallest number satisfying 1.d0+uround>1.d0
    if (work[0] == 0.0)
    {
        uround = 2.220446049250313e-16;
    } else {
        uround = work[0];
        if ((uround <= 1e-35) || (uround >= 1.0))
        {
            arret = 14;
        }
    }

    // -------  safety factor -------------
    if (work[1] == 0.0)
    {
        safe = 0.9;
    } else {
        safe = work[1];
        if ((safe <= 1e-4) || (safe >= 1.0))
        {
            arret = 15;
        }
    }

    // -------  fac1,fac2, parameters for step size selection
    fac1 = (work[2] == 0.0 ? 0.2 : work[2]);
    fac2 = (work[3] == 0.0 ? 10.0 : work[3]);

    // --------- beta for step control stabilization -----------
    if (work[4] <= 0.0)
    {
        beta = (work[4] < 0.0 ? 0.0 : 0.04);
    } else {
        beta = work[4];
        if (beta > 0.2)
        {
            arret = 16;
        }
    }
    // -------- maximal step size
    if (work[5] == 0.0) { hmax = xend - x; } else { hmax = work[5]; }
    // -------- initial step size
    h = work[6];
    // -------- when a fail has occurred, we return with wrong input error code
    if (arret > 0) { *ierr = -arret; return; }

    // -------- call to core integrator ------------
    dopcor(n, fcn, x, y, xend, hmax, h, rtol, atol, itol, solout, iout, ierr,
           nmax, uround, meth, nstiff, safe, beta, fac1, fac2, work, iwork,
           &nrdens, rpar, ipar, &nfcn, &nstep, &naccept, &nreject);

    work[6] = h;
    iwork[16] = nfcn;
    iwork[17] = nstep;
    iwork[18] = naccept;
    iwork[19] = nreject;

    return;
}

static void
dopcor(const int n, dopri_fcn* fcn, double *x, double* y, double *xend,
       double hmax, double h, double* rtol, double* atol, const int itol,
       dopri_solout* solout, const int iout, int* ierr, const int nmax,
       const double uround, const int meth, int nstiff, const double safe,
       const double beta, const double fac1, const double fac2, double* work,
       int* iwork, int* nrd, double* rpar, int* ipar, int* nfcn, int* nstep,
       int* naccept, int* nreject)
{

    // core integrator for DOPRI5
    // parameters same as in DOPRI5 with workspace added

    // Pointers into work
    double* restrict y1   = &work[20];
    double* restrict k1   = &work[20 + n];
    double* restrict k2   = &work[20 + n*2];
    double* restrict k3   = &work[20 + n*3];
    double* restrict k4   = &work[20 + n*4];
    double* restrict k5   = &work[20 + n*5];
    double* restrict k6   = &work[20 + n*6];
    double* restrict ysti = &work[20 + n*7];
    double* restrict cont = &work[20 + n*8];
    int* restrict icomp = &iwork[20];
    (void)meth;
    const double C2  = 0.2;
    const double C3  = 0.3;
    const double C4  = 0.8;
    const double C5  = 8.0 / 9.0;
    const double A21 = 0.2;
    const double A31 = 3.0 / 40.0;
    const double A32 = 9.0 / 40.0;
    const double A41 = 44.0 / 45.0;
    const double A42 = -56.0 / 15.0;
    const double A43 = 32.0 / 9.0;
    const double A51 = 19372.0 / 6561.0;
    const double A52 = -25360.0 / 2187.0;
    const double A53 = 64448.0 / 6561.0;
    const double A54 = -212.0 / 729.0;
    const double A61 = 9017.0 / 3168.0;
    const double A62 = -355.0 / 33.0;
    const double A63 = 46732.0 / 5247.0;
    const double A64 = 49.0 / 176.0;
    const double A65 = -5103.0 / 18656.0;
    const double A71 = 35.0 / 384.0;
    const double A73 = 500.0 / 1113.0;
    const double A74 = 125.0 / 192.0;
    const double A75 = -2187.0 / 6784.0;
    const double A76 = 11.0 / 84.0;
    const double E1  = 71.0 / 57600.0;
    const double E3  = -71.0 / 16695.0;
    const double E4  = 71.0 / 1920.0;
    const double E5  = -17253.0 / 339200.0;
    const double E6  = 22.0 / 525.0;
    const double E7  = -1.0 / 40.0;
    //  ---- DENSE OUTPUT OF SHAMPINE (1986)
    const double D1  = -12715105075.0 / 11282082432.0;
    const double D3  = 87487479700.0 / 32700410799.0;
    const double D4  = -10690763975.0 / 1880347072.0;
    const double D5  = 701980252875.0 / 199316789632.0;
    const double D6  = -1453857185.0 / 822651844.0;
    const double D7  = 69997945.0 / 29380423.0;

    int i, j, last, iasti, iord, reject, irtrn, nonsti = 0;
    double err, facold, expo1, facc1, facc2, posneg, atoli, rtoli, xold, hlamb, hout;
    double fac11, fac, hnew, stnum, stden, yd0, ydiff, bspl;

    facold = 1.0e-4;
    expo1 = 0.2 - beta * 0.75;
    facc1 = 1.0 / fac1;
    facc2 = 1.0 / fac2;
    posneg = copysign(1.0, *xend - *x);

    atoli = atol[0];
    rtoli = rtol[0];
    last = 0;
    hlamb = 0.0;
    iasti = 0;
    fcn(n, *x, y, k1, rpar, ipar);
    hmax = fabs(hmax);
    iord = 5;

    if (h == 0.0)
    {
        h = hinit(n, fcn, x, y, posneg, k1, k2, k3, iord, hmax, atol, rtol, itol, rpar, ipar);
    }
    *nfcn += 2;
    reject = 0;
    xold = *x;

    if (iout != 0) {
        irtrn = 1;
        hout = h;
        solout(*naccept + 1, xold, *x, y, n, cont, icomp, *nrd, rpar, ipar, &irtrn);
        if (irtrn < 0) { *ierr = 2; return; }
    } else {
        irtrn = 0;
    }

    // basic integration step
    while (1) {
        if (*nstep > nmax) { *ierr = -2; return; }
        if (0.1 * fabs(h) <= fabs(*x)*uround) { *ierr = -3; return; }
        if ((*x + 1.01*h - *xend) * posneg > 0.0) {
            h = *xend - *x;
            last = 1;
        }
        *nstep += 1;

        // the first 6 stages
        if (irtrn >= 2) { fcn(n, *x, y, k1, rpar, ipar); }

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*A21*k1[i]; }  // 22
        fcn(n, *x + C2*h, y1, k2, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(A31*k1[i] + A32*k2[i]); }  // 23
        fcn(n, *x + C3*h, y1, k3, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(A41*k1[i] + A42*k2[i] + A43*k3[i]); }  // 24
        fcn(n, *x + C4*h, y1, k4, rpar, ipar);

        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(A51*k1[i] + A52*k2[i] + A53*k3[i] + A54*k4[i]); }  // 25
        fcn(n, *x + C5*h, y1, k5, rpar, ipar);

        for (i = 0; i < n; i++) { ysti[i] = y[i] + h*(A61*k1[i] + A62*k2[i] + A63*k3[i] + A64*k4[i] + A65*k5[i]); }  // 26

        double xph = *x + h;
        fcn(n, xph, ysti, k6, rpar, ipar);
        for (i = 0; i < n; i++) { y1[i] = y[i] + h*(A71*k1[i] + A73*k3[i] + A74*k4[i] + A75*k5[i] + A76*k6[i]); }  // 27
        fcn(n, xph, y1, k2, rpar, ipar);

        if (iout >= 2)
        {
            for (j = 0; j < *nrd; j++)
            {
                i = icomp[j];
                cont[4*(*nrd) + j] = h*(D1*k1[i] + D3*k3[i] + D4*k4[i] + D5*k5[i] + D6*k6[i] + D7*k2[i]);
            }
            // 40
        }

        for (i = 0; i < n; i++) { k4[i] = (E1*k1[i] + E3*k3[i] + E4*k4[i] + E5*k5[i] + E6*k6[i] + E7*k2[i])*h; }  // 28
        *nfcn += 6;

        // error estimation
        err = 0.0;
        if (itol == 0)
        {
            for (i = 0; i < n; i++)
            {
                double sk = atoli + rtoli * fmax(fabs(y[i]), fabs(y1[i]));
                err += pow(k4[i] / sk, 2);
            }
            // 41
        } else {
            for (i = 0; i < n; i++)
            {
                double sk = atol[i] + rtol[i] * fmax(fabs(y[i]), fabs(y1[i]));
                err += pow(k4[i] / sk, 2);
            }
            // 42
        }
        err = sqrt(err / n);

        // computation of hnew
        fac11 = pow(err, expo1);

        // lund-stabilization
        fac = fac11 / pow(facold, beta);

        // we require fac1 <= hnew/h <= fac2
        fac = fmax(facc2, fmin(facc1, fac / safe));
        hnew = h / fac;

        if (err <= 1.0) {
            // step is accepted
            facold = fmax(err, 1.0e-4);
            *naccept += 1;

            // stiffness detection
            if (((*naccept) % nstiff == 0) || (iasti > 0)) {
                stnum = 0.0;
                stden = 0.0;
                for (int i = 0; i < n; i++) {
                    stnum += pow(k2[i] - k6[i], 2);
                    stden += pow(y1[i] - ysti[i], 2);
                }
                if (stden > 0.0) hlamb = h * sqrt(stnum / stden);
                if (hlamb > 3.25)
                {
                    nonsti = 0;
                    iasti += 1;
                    if (iasti == 15) { *ierr = -4; return; }
                } else {
                    nonsti += 1;
                    if (nonsti == 6) { iasti = 0; }
                }
            }

            if (iout >= 2)
            {
                for (j = 0; j < *nrd; j++)
                {
                    i = icomp[j];
                    yd0 = y[i];
                    ydiff = y1[i] - yd0;
                    bspl = h * k1[i] - ydiff;
                    cont[j] = y[i];
                    cont[*nrd + j] = ydiff;
                    cont[2 * (*nrd) + j] = bspl;
                    cont[3 * (*nrd) + j] = -h * k2[i] + ydiff - bspl;
                }
                // 43
            }

            for (i = 0; i < n; i++)
            {
                k1[i] = k2[i];
                y[i] = y1[i];
            }
            // 44
            xold = *x;
            *x = xph;

            if (iout != 0)
            {
                hout = h;
                solout(*naccept + 1, xold, *x, y, n, cont, icomp, *nrd, rpar, ipar, &irtrn);
                if (irtrn < 0) { *ierr = 2; return; }
            }

            // normal exit
            if (last)
            {
                h = hnew;
                *ierr = 1;
                return;
            }

            if (fabs(hnew) > hmax) { hnew = posneg * hmax; }
            if (reject) { hnew = posneg * fmin(fabs(hnew), fabs(h)); }
            reject = 0;

        } else {
            // step is rejected
            hnew = h / fmin(facc1, fac11 / safe);
            reject = 1;
            if (*naccept >= 1) { *nreject += 1; }
            last = 0;
        }
        h = hnew;
    }
    return;
}


static double
hinit(int n, dopri_fcn* fcn, double *x, double* y, double posneg, double* f0,
      double* f1, double* y1, int iord, double hmax, double* atol, double* rtol,
      int itol, double* rpar, int* ipar)
{

    double dnf = 0.0, dny = 0.0;
    double atoli = atol[0], rtoli = rtol[0], sk, h, der2, der12, h1;

    // compute a first guess for explicit euler as
    // h = 0.01 * norm (y0) / norm (f0)
    // the increment for explicit euler is small compared to the solution
    if (itol == 0)
    {
        for (int i = 0; i < n; i++)
        {
            sk = atoli + rtoli * fabs(y[i]);
            dnf = dnf + pow(f0[i] / sk, 2);
            dny = dny + pow(y[i] / sk, 2);
        }
        // 10
    } else {
        for (int i = 0; i < n; i++)
        {
            sk = atol[i] + rtol[i] * fabs(y[i]);
            dnf = dnf + pow(f0[i] / sk, 2);
            dny = dny + pow(y[i] / sk, 2);
        }
        // 11
    }

    if ((dnf <= 1.0e-10) || (dny <= 1.0e-10))
    {
        h = 1.0e-6;
    } else {
        h = sqrt(dny / dnf) * 0.01;
    }

    h = fmin(h, hmax);
    h = copysign(h, posneg);

    // perform an explicit euler step
    for (int i = 0; i < n; i++)
    {
        y1[i] = y[i] + h * f0[i];
    }
    // 12
    fcn(n, *x + h, y1, f1, rpar, ipar);

    // estimate the second derivative of the solution
    der2 = 0.0;
    if (itol == 0)
    {
        for (int i = 0; i < n; i++)
        {
            sk = atoli + rtoli * fabs(y[i]);
            der2 = der2 + pow((f1[i] - f0[i]) / sk, 2);
        }
        // 15
    } else {
        for (int i = 0; i < n; i++)
        {
            sk = atol[i] + rtol[i] * fabs(y[i]);
            der2 = der2 + pow((f1[i] - f0[i]) / sk, 2);
        }
        // 16
    }
    der2 = sqrt(der2) / h;

    // step size is computed such that
    // h**iord * max(norm(f0), norm(der2)) = 0.01
    der12 = fmax(fabs(der2), sqrt(dnf));
    if (der12 <= 1.0e-15)
    {
        h1 = fmax(1.0e-6, fabs(h) * 1.0e-3);
    } else {
        h1 = pow(0.01 / der12, 1.0 / iord);
    }

    h = fmin(100.0 * fabs(h), fmin(h1, hmax));
    return copysign(h, posneg);
}
