#include <math.h>
#include <complex.h>
#include "_amos.h"


int acai(
    double complex z,
    double fnu,
    int kode,
    int mr,
    int n,
    double complex *y,
    double rl,
    double tol,
    double elim,
    double alim
) {
    double complex csgn, cspn, c1, c2, zn, cy[2];
    double arg, ascle, az, cpn, dfnu, fmr, sgn, spn, yy;
    int inu, iuf, nn, nw;
    double pi = 3.14159265358979324;
    int nz = 0;
    zn = -z;
    az = cabs(z);
    nn = n;
    dfnu = fnu + (n-1);
    if (az > 2.) {
        if (az*az*0.25 > dfnu+1.0) {
            goto L10;
        }
    }
    //
    // POWER SERIES FOR THE I FUNCTION
    //
    seri(zn, fnu, kode, nn, y, tol, elim, alim);
    goto L40;
L10:
    if (az >= rl) {
        //
        // ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
        //
        nw = asyi(zn, fnu, kode, nn, y, rl, tol, elim, alim);
        if (nw < 0) {
            nz = -1;
            if (nw == -2) { nz = -2; }
            return nz;
        }
    } else {
        //
        // MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
        //
        nw = mlri(zn, fnu, kode, nn, y, tol);
        if (nw < 0) {
            nz = -1;
            if (nw == -2) { nz = -2; }
            return nz;
        }
    }
L40:
//
// ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
//
    nw = bknu(zn, fnu, kode, 1, &cy[0], tol, elim, alim);
    if (nw != 0) {
        nz = -1;
        if (nw == -2) { nz = -2; }
        return nz;
    }
    fmr = mr;
    sgn = (fmr < 0.0 ? pi : -pi);
    csgn = CMPLX(0.0, sgn);
    if (kode != 1) {
        yy = -cimag(zn);
        cpn = cos(yy);
        spn = sin(yy);
        csgn *= CMPLX(cpn, spn);
    }
    //
    // CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    // WHEN FNU IS LARGE
    //
    inu = (int)fnu;
    arg = (fnu - inu)*sgn;
    cpn = cos(arg);
    spn = sin(arg);
    cspn = CMPLX(cpn, spn);
    if (inu % 2 == 1) { cspn = -cspn; }
    c1 = cy[0];
    c2 = y[0];
    if (kode != 1) {
        iuf = 0;
        ascle = 1e3 * d1mach[0] / tol;
        nw = s1s2(zn, &c1, &c2, ascle, alim, &iuf);
        nz += nw;
    }
    y[0] = cspn*c1 + csgn*c2;
    return nz;
}


int acon(
    double complex z,
    double fnu,
    int kode,
    int mr,
    int n,
    double complex *y,
    double rl,
    double fnul,
    double tol,
    double elim,
    double alim
) {
    double complex ck, cs, cscl, cscr, csgn, cspn, c1, c2, rz, sc1, sc2, st,\
                   s1, s2, zn;
    double arg, ascle, as2, bscle, c1i, c1m, c1r, fmr, sgn, yy;
    int i, inu, iuf, kflag, nn, nw, nz;
    double pi = 3.14159265358979324;
    double complex cy[2] = { 0.0 };
    double complex css[3] = { 0.0 };
    double complex csr[3] = { 0.0 };
    double bry[3] = { 0.0 };

    nz = 0;
    zn = -z;
    nn = n;
    nw = binu(zn, fnu, kode, nn, y, rl, fnul, tol, elim, alim);
    if (nw >= 0) {
        //
        // ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
        //
        nn = (n > 2 ? 2 : n);
        nw = bknu(zn, fnu, kode, nn, cy, tol, elim, alim);
        if (nw == 0) {
            s1 = cy[0];
            fmr = mr;
            sgn = ( fmr < 0 ? -pi : pi );
            csgn = sgn*I;
            if (kode != 1) {
                yy = -cimag(zn);
                csgn *= CMPLX(cos(yy), sin(yy));
            }
            inu = (int)fnu;
            arg = (fnu - inu)*sgn;
            cspn = CMPLX(cos(arg), sin(arg));
            if (inu % 2 == 1) { cspn = -cspn; }
            iuf = 0;
            c1 = s1;
            c2 = y[0];
            ascle = 1e3*d1mach[0]/tol;
            if (kode != 1) {
                nw = s1s2(zn, &c1, &c2, ascle, alim, &iuf);
                nz += nw;
                sc1 = c1;
            }
            y[0] = cspn*c1 + csgn*c2;
            if (n == 1) { return nz; }
            cspn = -cspn;
            s2 = cy[1];
            c1 = s2;
            c2 = y[1];
            if (kode != 1) {
                nw = s1s2(zn, &c1, &c2, ascle, alim, &iuf);
                nz += nw;
                sc2 = c1;
            }
            y[1] = cspn*c1 + csgn*c2;
            if (n == 2) { return nz; }
            cspn = -cspn;
            rz = 2.0 / zn;
            ck = (fnu + 1.0)*rz;
            //
            // SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
            //
            cscl = 1.0 / tol;
            cscr = tol;
            css[0] = cscl;
            css[1] = 1.0;
            css[2] = cscr;
            csr[0] = cscr;
            csr[1] = 1.0;
            csr[2] = cscl;
            bry[0] = ascle;
            bry[1] = 1.0 / ascle;
            bry[2] = d1mach[1];
            as2 = cabs(s2);
            kflag = 2;
            if (as2 <= bry[0] ) {
                kflag = 1;
            } else {
                if (as2 >= bry[1]) {
                    kflag = 3;
                }
            }
            bscle = bry[kflag-1];
            s1 *= css[kflag-1];
            s2 *= css[kflag-1];
            cs = csr[kflag-1];
            for (i = 3; i < (n+1); i++)
            {
                st = s2;
                s2 = ck*s2 + s1;
                s1 = st;
                c1 = s2*cs;
                st = c1;
                c2 = y[i-1];
                if (kode != 1) {
                    if (iuf >= 0) {
                        nw = s1s2(zn, &c1, &c2, ascle, alim, &iuf);
                        nz += nw;
                        sc1 = sc2;
                        sc2 = c1;
                        if (iuf == 3){
                            iuf = -4;
                            s1 = sc1 * css[kflag-1];
                            s2 = sc2 * css[kflag-1];
                            st = sc2;
                        }
                    }
                }
                y[i-1] = cspn*c1 + csgn*c2;
                ck += rz;
                cspn = -cspn;
                if (kflag < 3) {
                    c1r = fabs(creal(c1));
                    c1i = fabs(cimag(c1));
                    c1m = fmax(c1r, c1i);
                    if (c1m > bscle) {
                        kflag += 1;
                        bscle = bry[kflag-1];
                        s1 *= cs;
                        s2 = st;
                        s1 *= css[kflag-1];
                        s2 *= css[kflag-1];
                        cs = csr[kflag-1];
                    }
                }
            }
            return nz;
        }
    }
    nz = -1;
    if (nw == -2) { nz = -2; }
    return nz;
}   


double complex airy(
    double complex z,
    int id,
    int kode,
    int *nz,
    int *ierr
) {
    double complex ai, csq, cy[1], s1, s2, trm1, trm2, zta, z3;
    double aa, ad, ak, alim, atrm, az, az3, bk, ck, dig, dk, d1, d2,\
           elim, fid, fnu, rl, r1m5, sfac, tol, zi, zr, bb, alaz;
    int iflag, k, k1, k2, mr, nn;
    double tth = 2. / 3.;
    double c1 = 0.35502805388781723926;  /* 1/(Gamma(2/3) * 3**(2/3)) */
    double c2 = 0.25881940379280679841;  /* 1/(Gamma(1/3) * 3**(1/3)) */
    double coef = 0.18377629847393068317;  /* 1 / (sqrt(3) * PI) */

    *ierr = 0;
    *nz = 0;
    ai = 0.;
    if ((id < 0) || (id > 1)) { *ierr = 1; }
    if ((kode < 1) || (kode > 2)) { *ierr = 1; }
    if (*ierr != 0) return 0.;
    az = cabs(z);
    tol = d1mach[3];
    fid = id;

    if (az <= 1.0) {
        //
        // POWER SERIES FOR ABS(Z) <= 1.
        //
        s1 = 1.0;
        s2 = 1.0;
        if (az < tol) {
            aa = 1e3*d1mach[0];
            s1 = 0.;
            if (id != 1) {
                if (az > aa) { s1 = c2 * z; }
                ai = c1 - s1;
                return ai;
            }
            ai = -c2;
            aa = sqrt(aa);
            if (az > aa) { s1 = z * z * 0.5; }
            ai += s1 * c1;
            return ai;
        }
        aa = az*az;
        if (aa >= tol/az) {
            trm1 = 1.0;
            trm2 = 1.0;
            atrm = 1.0;
            z3 = z*z*z;
            az3 = az * aa;
            ak = 2.0 + fid;
            bk = 3.0 - fid - fid;
            ck = 4.0 - fid;
            dk = 3.0 + fid + fid;
            d1 = ak * dk;
            d2 = bk * ck;
            ad = (d1 > d2 ? d2 : d1);
            ak = 24.0 + 9.0*fid;
            bk = 30.0 - 9.0*fid;
            for (int k = 1; k < 26; k++)
            {
                trm1 *= z3/d1;
                s1 += trm1;
                trm2 *= z3/d2;
                s2 += trm2;
                atrm *= az3 / ad;
                d1 += ak;
                d2 += bk;
                ad = (d1 > d2 ? d2 : d1);
                if (atrm < tol*ad) { break; }
                ak += 18.0;
                bk += 18.0;
            }
        }
        if (id != 1) {
            ai = s1*c1 - z*s2*c2;
            if (kode == 1) { return ai; }
            zta = z*csqrt(z)*tth;
            ai *= cexp(zta);
            return ai;
        }
        ai = -s2*c2;
        if (az > tol) { ai += z*z*s1*c1/(1. + fid); }
        if (kode == 1) { return ai; }
        zta = z*csqrt(z)*tth;
        return ai*cexp(zta);
    }
    //
    // CASE FOR ABS(Z) > 1.0
    //
    fnu = (1.0 + fid) / 3.0;
    //
    // SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    // TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    // ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    // EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    // EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    // UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    // RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    // DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    //
    k1 = i1mach[14];
    k2 = i1mach[15];
    r1m5 = d1mach[4];
    k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
    elim = 2.303 * (k*r1m5 - 3.0);
    k1 = i1mach[13] - 1;
    aa = r1m5*k1;
    dig = (aa > 18.0 ? 18.0 : aa);
    aa *= 2.303;
    alim = elim + (-aa > -41.45 ? -aa : -41.45);
    rl = 1.2*dig + 3.0;
    alaz = log(az);
    // 
    // TEST FOR RANGE
    // 
    aa = 0.5 / tol;
    bb = i1mach[8] * 0.5;
    aa = (aa > bb ? bb : aa);
    aa = pow(aa, tth);
    if (az > aa) {
        *ierr = 4;
        *nz = 0;
        return 0.;
    }
    aa = sqrt(aa);
    if (az > aa) { *ierr = 3; }
    csq = csqrt(z);
    zta = z * csq * tth;
    //
    // RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
    // 
    iflag = 0;
    sfac = 1.0;
    zi = cimag(z);
    zr = creal(z);
    ak = cimag(zta);
    if (zr < 0.0) {
        bk = creal(zta);
        ck = -fabs(bk);
        zta = CMPLX(ck, ak);
    }
    if ((zi == 0.0) && (zr <= 0.0)) {
        zta = CMPLX(0.0, ak);
    }
    aa = creal(zta);
    if ((aa < 0.0) || (zr <= 0.0)) {
        if (kode != 2) {
            //
            // OVERFLOW TEST
            //
            if (aa <= -alim) {
                aa = -aa + 0.25*alaz;
                iflag = 1;
                sfac = tol;
                if (aa > elim) {
                    /* 270 */
                    *nz = 0;
                    *ierr = 2;
                    return ai;
                }
            }
        }
        //
        // CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
        //
        mr = 1;
        if (zi < 0.0) { mr = -1; }
        nn = acai(zta, fnu, kode, mr, 1, &cy[0], rl, tol, elim, alim);
        if (nn < 0) {
            if (nn == -1) {
                *nz = 1;
                return 0.;
            } else {
                *nz = 0;
                *ierr = 5;
                return 0.;
            }
        }
        *nz += nn;
    } else {
        if (kode != 2) {
            //
            // UNDERFLOW TEST
            // 
            if (aa >= alim) {
                aa = -aa - 0.25 * alaz;
                iflag = 2;
                sfac = 1.0 / tol;
                if (aa < -elim) {
                    *nz = 1;
                    return 0.;
                }
            }
        }
        *nz = bknu(zta, fnu, kode, 1, &cy[0], tol, elim, alim);
    }
    s1 = cy[0]*coef;
    
    if (iflag == 0) {
        if (id != 1) {
            return csq *s1;
        }
        return (-z*s1);
    }
    s1 *= sfac;
    if (id != 1) {
        s1 *= csq;
        return (s1/sfac);
    }
    s1 *= -z;
    return (s1/sfac);
}


int asyi(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    double rl,
    double tol,
    double elim,
    double alim
) {
    double complex ak1, ck, cs1, cs2, cz, dk, ez, p1, rz, s2;
    double aa, acz, aez, ak, arg, arm, atol, az, bb, bk, dfnu;
    double dnu2, fdn, rtr1, s, sgn, sqk, x, yy;
    int ib, il, inu, j, jl, k, koded, m, nn;
    double pi = 3.14159265358979324;
    double rpi = 0.159154943091895336; /* (1 / pi) */
    int nz = 0;
    az = cabs(z);
    x = creal(z);
    arm = 1e3*d1mach[0];
    rtr1 = sqrt(arm);
    il = (n > 2 ? 2 : n);
    dfnu = fnu + (n - il);
    // OVERFLOW TEST
    ak1 = csqrt(rpi / z);
    cz = z;
    if (kode == 2) { cz = z - x; }
    acz = creal(cz);
    if (fabs(acz) <= elim) {
        dnu2 = dfnu + dfnu;
        koded = 1;
        if (!((fabs(acz) > alim) && (n > 2))) {
            koded = 0;
            ak1 *= cexp(cz);
        }
        fdn = 0.;
        if (dnu2 > rtr1) { fdn = dnu2 * dnu2; }
        ez = z * 8.;
        // WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE
        // RELATIVE TO THE FIRST RECIPROCAL POWER SINCE THIS
        // IS THE LEADING TERM OF THE EXPANSION FOR THE
        // IMAGINARY PART.
        aez = 8. * az;
        s = tol / aez;
        jl = rl + rl + 2;
        yy = cimag(z);
        p1 = 0.;
        if (yy != 0.) {
            inu = (int)fnu;
            arg = (fnu - inu) * pi;
            inu += n - il;
            ak = -sin(arg);
            bk = cos(arg);
            if (yy < 0.) { bk = -bk; }
            p1 = ak + bk*I;
            if (inu % 2 == 1) { p1 = -p1; }
        }
        for (int k = 1; k < (il+1); k++)
        {
            sqk = fdn - 1.;
            atol = s*fabs(sqk);
            sgn = 1.;
            cs1 = 1.;
            cs2 = 1.;
            ck = 1.;
            ak = 0.;
            aa = 1.;
            bb = aez;
            dk = ez;
            j = 1;
            for (int j = 1; j < (jl+1); j++)
            {
                ck *= sqk / dk;
                cs2 += ck;
                sgn = -sgn;
                cs1 += ck*sgn;
                dk += ez;
                aa *= fabs(sqk) / bb;
                bb += aez;
                ak += 8.;
                sqk -= ak;
                if (aa <= atol) { break; }
            }
            if ((j == jl) && (aa > atol)) { return -2; }

            /* 20 */
            s2 = cs1;
            if (x + x < elim) { s2 += p1*cs2*cexp(-z-z); }
            fdn += 8. * dfnu + 4.;
            p1 = -p1;
            m = n - il + k;
            y[m - 1] = s2 * ak1;
        }
        if (n <= 2) { return nz; }
        nn = n;
        k = nn - 2;
        ak = k;
        rz = 2. / z;
        ib = 3;
        for (int i = ib; i < (nn+1); i++)
        {
            y[i-1] = (ak + fnu)*rz*y[k] + y[k+1];
            ak -= 1.;
            k -=1;
        }
        if (koded == 0) { return nz; }
        ck = cexp(cz);
        for (int i = 0; i < (nn + 1); i++) { y[i] *= ck; }
    }
    return -1;
}


int besh(
    double complex z,
    double fnu,
    int kode,
    int m,
    int n,
    double complex *cy,
    int *ierr
) {
    double complex zn, zt, csgn;
    double aa, alim, aln, arg, az, cpn, dig, elim, fmm, fn, fnul,
        rhpi, rl, r1m5, sgn, spn, tol, ufl, xn, xx, yn, yy,
        bb, ascle, rtol, atol;
    int i, inu, inuh, ir, k, k1, k2, mm, mr, nn, nuf, nw, nz;

    double hpi = 1.57079632679489662; /* 0.5 PI */

    nz = 0;
    xx = creal(z);
    yy = cimag(z);
    *ierr = 0;

    if ((xx == 0.0) && (yy == 0.0)) { *ierr = 1; }
    if (fnu < 0.0) { *ierr = 1; }
    if ((m < 1) || (m > 2)) { *ierr = 1; }
    if ((kode < 1) || (kode > 2)) { *ierr = 1; }
    if (n < 1) { *ierr = 1; }
    if (*ierr != 0) { return nz; }
    nn = n;
    //
    //  SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    //  TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    //  ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    //  EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    //  EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    //  UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    //  RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    //  DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    //  FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
    //
    tol = fmax(d1mach[3], 1e-18);
    k1 = i1mach[14];
    k2 = i1mach[15];
    r1m5 = d1mach[4];
    k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
    elim = 2.303 * (k*r1m5 - 3.0);
    k1 = i1mach[13] - 1;
    aa = r1m5*k1;
    dig = (aa > 18.0 ? 18.0 : aa);
    aa *= 2.303;
    alim = elim + (-aa > -41.45 ? -aa : -41.45);
    fnul = 10.0 + 6.0 * (dig - 3.0);
    rl = 1.2*dig + 3.0;
    fn = fnu + (nn - 1);
    mm = 3 - m - m;
    fmm = mm;
    zn = z * (-fmm * I);
    xn = creal(zn);
    yn = cimag(zn);
    az = cabs(z);

    //
    // TEST FOR RANGE
    //
    aa = 0.5 / tol;
    bb = d1mach[1] * 0.5;
    aa = fmin(aa, bb);
    if (az <= aa) {
        if (fn <= aa) {
            aa = sqrt(aa);
            if (az > aa) { *ierr = 3; }
            if (fn > aa) { *ierr = 3; }
            //
            // OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
            //
            ufl = d1mach[0] * 1.0e3;
            if (az >= ufl) {
                if (fnu <= fnul) {
                    if (fn > 1.0) {
                        if (fn <= 2.0) {
                            if (az > tol) { goto L10; }
                            arg = 0.5 * az;
                            aln = -fn * log(arg);
                            if (aln > elim) { goto L50; }
                        } else {
                            nuf = uoik(zn, fnu, kode, 2, nn, cy, tol, elim, alim);
                            if (nuf < 0) { goto L50; }
                            nz += nuf;
                            nn -= nuf;
                            //
                            // HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
                            // IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
                            //
                            if (nn == 0) goto L40;
                        }
                    }
L10:
                    if (!((xn < 0.0) || ((xn == 0.0) && (yn < 0.0) && (m == 2)))) {
                        //
                        // RIGHT HALF PLANE COMPUTATION, XN >= 0.  .AND.  (XN.NE.0.  .OR.
                        // YN >= 0.  .OR.  M=1)
                        //
                        nz = bknu(zn, fnu, kode, nn, cy, tol, elim, alim);
                        goto L20;
                    }
                    //
                    // LEFT HALF PLANE COMPUTATION
                    //
                    mr = -mm;
                    nw = acon(zn, fnu, kode, mr, nn, cy, rl, fnul, tol, elim, alim);
                    if (nw < 0) { goto L60; }
                    nz = nw;
                } else {
                    //
                    // UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
                    //
                    mr = 0;
                    if ((xn < 0.0) || ((xn == 0.0) && (yn < 0.0) && (m == 2))) {
                        mr = -mm;
                        if ((xn == 0.0) && (yn < 0.0)) { zn = -zn; }
                    }
                    nw = bunk(zn, fnu, kode, mr, nn, cy, tol, elim, alim);
                    if (nw < 0) { goto L60; }
                    nz += nw;
                }
                //
                // H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
                //
                // ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
                //
L20:
                sgn = copysign(hpi, -fmm);
                //
                // CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
                // WHEN FNU IS LARGE
                //
                inu = (int)fnu;
                inuh = inu / 2;
                ir = inu - 2 * inuh;
                arg = (fnu - (inu - ir)) * sgn;
                rhpi = 1.0 / sgn;
                cpn = rhpi * cos(arg);
                spn = rhpi * sin(arg);
                csgn = CMPLX(-spn, cpn);
                if (inuh % 2 == 1) { csgn = -csgn; }
                zt = -fmm*I;
                rtol = 1.0 / tol;
                ascle = ufl * rtol;
                for (i = 1; i < (nn+1); i++) {
                    zn = cy[i-1];
                    aa = creal(zn);
                    bb = cimag(zn);
                    atol = 1.0;
                    if (fmax(fabs(aa), fabs(bb)) <= ascle) {
                        zn *= rtol;
                        atol = tol;
                    }
                    zn *= csgn;
                    cy[i-1] = zn * atol;
                    csgn = zt;
                }
                return nz;
L40:
                if (xn >= 0.0) return nz;
            }
L50:
            *ierr = 2;
            return 0;
L60:
            if (nw == -1) goto L50;
            *ierr = 5;
            return 0;
        }
    }
    *ierr = 4;
    return 0;
}


int besi(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *cy,
    int *ierr
) {
    double complex csgn, zn;
    double aa, alim, arg, dig, elim, fnul, rl, r1m5, s1, s2, tol, xx, yy, az,\
           fn, bb, ascle, rtol, atol;
    int i, inu, k, k1, k2, nn, nz;
    double pi = 3.14159265358979324;

    *ierr = 0;
    nz = 0;
    if (fnu < 0.0) { *ierr = 1; }
    if ((kode < 1) || (kode > 2)) { *ierr = 1; }
    if (n < 1) { *ierr = 1; }
    if (*ierr != 0) { return nz; }
    xx = creal(z);
    yy = cimag(z);
    //
    //  SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    //  TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    //  ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    //  EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    //  EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    //  UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    //  RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    //  DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    //  FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
    //
    tol = fmax(d1mach[3], 1e-18);
    k1 = i1mach[14];
    k2 = i1mach[15];
    r1m5 = d1mach[4];
    k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
    elim = 2.303 * (k*r1m5 - 3.0);
    k1 = i1mach[13] - 1;
    aa = r1m5*k1;
    dig = (aa > 18.0 ? 18.0 : aa);
    aa *= 2.303;
    alim = elim + (-aa > -41.45 ? -aa : -41.45);
    rl = 1.2 * dig + 3.0;
    fnul = 10.0 + 6.0 * (dig - 3.0);
    az = cabs(z);
    //
    // TEST FOR RANGE
    // 
    aa = 0.5 / tol;
    bb = d1mach[1]*0.5;
    aa = fmin(aa, bb);
    if (az <= aa) {
        fn = fnu + (n - 1);
        if (fn <= aa) {
            aa = sqrt(aa);
            if (az > aa) { *ierr = 3; }
            if (fn > aa) { *ierr = 3; }
            zn = z;
            csgn = 1.0;
            if (xx < 0.0) {
                zn = -z;
                //
                // CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
                // WHEN FNU IS LARGE
                //
                inu = (int)fnu;
                arg = (fnu - inu)*pi;
                if (yy < 0.0) { arg = -arg; }
                s1 = cos(arg);
                s2 = sin(arg);
                csgn = CMPLX(s1, s2);
                if (inu % 2 == 1) { csgn = -csgn; }
            }
            nz = binu(zn, fnu, kode, n, cy, rl, fnul, tol, elim, alim);
            if (nz >= 0) {
                if (xx > 0.0) { return nz; }
                //
                // ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
                //
                nn = n - z;
                if (nn == 0) { return nz; }
                rtol = 1.0 / tol;
                ascle = d1mach[0]*rtol*1e3;
                for (i = 1; i < (nn+1); i++)
                {
                    zn = cy[i-1];
                    aa = creal(zn);
                    bb = cimag(zn);
                    atol = 1.0;
                    if (fmax(fabs(aa), fabs(bb)) <= ascle) {
                        zn *= rtol;
                        atol = tol;
                    }
                    zn *= csgn;
                    cy[i-1] = zn*atol;
                    csgn = -csgn;
                }
                return nz;
            }
            if (nz != 2) {
                *ierr = 2;
                return 0;
            }
            *ierr = 5;
            return 0;
        }
    }
    *ierr = 4;
    return 0;
}


int besj(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *cy,
    int *ierr
) {
    double complex ci, csgn, zn;
    double aa, alim, arg, dig, elim, fnul, rl, r1, r1m5, r2,
        tol, yy, az, fn, bb, ascle, rtol, atol;
    int i, inu, inuh, ir, k1, k2, nl, nz, k;
    double hpi = 1.570796326794896619;

    *ierr = 0;
    nz = 0;
    if (fnu < 0.0) *ierr = 1;
    if (kode < 1 || kode > 2) *ierr = 1;
    if (n < 1) *ierr = 1;
    if (*ierr != 0) return nz;

    tol = fmax(d1mach[3], 1e-18);
    k1 = i1mach[14];
    k2 = i1mach[15];
    r1m5 = d1mach[4];
    k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
    elim = 2.303 * (k*r1m5 - 3.0);
    k1 = i1mach[13] - 1;
    aa = r1m5*k1;
    dig = (aa > 18.0 ? 18.0 : aa);
    aa *= 2.303;
    alim = elim + (-aa > -41.45 ? -aa : -41.45);
    fnul = 10.0 + 6.0 * (dig - 3.0);
    rl = 1.2*dig + 3.0;
    ci = I;
    yy = cimag(z);
    az = cabs(z);

    //-----------------------------------------------------------------------
    //     TEST FOR RANGE
    //-----------------------------------------------------------------------
    aa = 0.5 / tol;
    bb = d1mach[1] * 0.5;
    aa = fmin(aa, bb);
    fn = fnu + (n - 1);
    if (az <= aa) {
        if (fn <= aa) {
            aa = sqrt(aa);
            if (az > aa) *ierr = 3;
            if (fn > aa) *ierr = 3;
            //-----------------------------------------------------------------------
            //     CALCULATE CSGN = EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
            //     WHEN FNU IS LARGE
            //-----------------------------------------------------------------------
            inu = (int)fnu;
            inuh = inu / 2;
            ir = inu - 2 * inuh;
            arg = (fnu - (inu - ir)) * hpi;
            r1 = cos(arg);
            r2 = sin(arg);
            csgn = CMPLX(r1, r2);
            if (inuh % 2 == 1) { csgn = -csgn; }
            //-----------------------------------------------------------------------
            //     ZN IS IN THE RIGHT HALF PLANE
            //-----------------------------------------------------------------------
            zn = -z * ci;
            if (yy < 0.0) {
                zn = -zn;
                csgn = conj(csgn);
                ci = conj(ci);
            }
            nz = binu(zn, fnu, kode, n, cy, rl, fnul, tol, elim, alim);
            if (nz >= 0) {
                nl = n - nz;
                if (nl == 0) { return nz; }
                rtol = 1.0 / tol;
                ascle = d1mach[0]*rtol*1e3;
                for (i = 1; i < (nl+1); i++) 
                {
                    zn = cy[i-1];
                    aa = creal(zn);
                    bb = cimag(zn);
                    atol = 1.0;
                    if (fmax(fabs(aa), fabs(bb)) <= ascle) {
                        zn *= rtol;
                        atol = tol;
                    }
                    zn *= csgn;
                    cy[i-1] = zn * atol;
                    csgn = csgn * ci;
                }
                return nz;
            }
            if (nz != -2) {
                *ierr = 2;
                return 0;
            }
            *ierr = 5;
            return 0;
        }
    }
    *ierr = 4;
    return 0;
}


int besk(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *cy,
    int *ierr
) {
    double xx = creal(z);
    double yy = cimag(z);
    double aa, alim, aln, arg, az, dig, elim, fn, fnul, rl, r1m5, tol, ufl, bb;
    int k, k1, k2, mr, nn, nuf, nw, nz;

    *ierr = 0;
    nz = 0;

    if ((yy == 0.0) && (xx == 0.0)) { *ierr = 1; }
    if (fnu < 0.0) { *ierr = 1; }
    if (kode < 1 || kode > 2) { *ierr = 1; }
    if (n < 1) { *ierr = 1; }
    if (*ierr != 0) { return nz; }

    nn = n;
    //
    //  SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    //  TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    //  ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    //  EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    //  EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    //  UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    //  RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    //  DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    //  FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
    //
    tol = fmax(d1mach[3], 1e-18);
    k1 = i1mach[14];
    k2 = i1mach[15];
    r1m5 = d1mach[4];
    k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
    elim = 2.303 * (k*r1m5 - 3.0);
    k1 = i1mach[13] - 1;
    aa = r1m5*k1;
    dig = (aa > 18.0 ? 18.0 : aa);
    aa *= 2.303;
    alim = elim + (-aa > -41.45 ? -aa : -41.45);
    fnul = 10.0 + 6.0 * (dig - 3.0);
    rl = 1.2 * dig + 3.0;
    az = cabs(z);
    fn = fnu + (nn - 1);
    //
    // TEST FOR RANGE
    //
    aa = 0.5 / tol;
    bb = d1mach[1] * 0.5;
    aa = fmin(aa, bb);
    if (az <= aa) {
        if (fn <= aa) {
            aa = sqrt(aa);
            if (az > aa) { *ierr = 3; }
            if (fn > aa) { *ierr = 3; }
            //
            // OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
            //
            ufl = d1mach[0] * 1.0E+3;
            if (az >= ufl) {
                if (fnu <= fnul) {
                    if (fn > 1.0) {
                        if (fn <= 2.0) {
                            if (az > tol) { goto L10; }
                            arg = 0.5 * az;
                            aln = -fn * log(arg);
                            if (aln > elim) { goto L30; }
                        } else {
                            nuf = uoik(z, fnu, kode, 2, nn, cy, tol, elim, alim);
                            if (nuf < 0) { goto L30; }
                            nz += nuf;
                            nn -= nuf;
                            //
                            // HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
                            // IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
                            //
                            if (nn == 0) { goto L20; }
                        }
                    }
L10:
                    if (xx >= 0.0) {
                        //
                        // RIGHT HALF PLANE COMPUTATION, REAL(Z) >= 0.
                        //
                        nw = bknu(z, fnu, kode, nn, cy, tol, elim, alim);
                        if (nw < 0) { goto L40; }
                        return nw;
                    }
                    //
                    // LEFT HALF PLANE COMPUTATION
                    // PI/2 < ARG(Z) <= PI AND -PI < ARG(Z) < -PI/2.
                    //
                    if (nz != 0) { goto L30; }
                    mr = 1;
                    if (yy < 0.0) { mr = -1; }
                    nw = acon(z, fnu, kode, mr, nn, cy, rl, fnul, tol, elim, alim);
                    if (nw < 0) { goto L40; }
                    return nw;
                }
                // Uniform asymptotic expansions for fnu > fnul
                mr = 0;
                if (xx < 0.0) {
                    mr = 1;
                    if (yy < 0.0) { mr = -1; }
                }
                nw = bunk(z, fnu, kode, mr, nn, cy, tol, elim, alim);
                if (nw < 0) { goto L40; }
                nz += nw;
                return nz;
L20:            
                if (xx > 0.0) { return nz ; }
            }
L30:
            *ierr = 2;
            return 0;
L40:
            if (nw == -1) { goto L30; }
            *ierr = 5;
            return 0;
        }
    }
    *ierr = 4;
    return 0;
}


int besy(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *cy,
    int *ierr
    ) {
    double complex ci, csgn, cspn, cwrk[n], ex, zu, zv, zz, zn;
    double arg, elim, ey, r1, r2, tay, xx, yy, ascle, rtol, atol, tol, aa, bb,
        ffnu, rhpi, r1m5;
    int i, ifnu, k, k1, k2, nz, nz1, nz2, i4;
    double complex cip[4] = { 1.0, I, -1.0, -I };
    double hpi = 1.57079632679489662; /* 0.5 PI */

    xx = creal(z);
    yy = cimag(z);
    *ierr = 0;
    nz = 0;
    if ((xx == 0.0) && (yy == 0.0)) { *ierr = 1; }
    if (fnu < 0.0) { *ierr = 1; }
    if ((kode < 1) || (kode > 2)) { *ierr = 1; }
    if (n < 1) { *ierr = 1; }
    if (*ierr != 0) { return nz; }
    ci = I;
    zz = z;
    if (yy < 0.0) { zz = conj(z); }
    zn = -ci * zz;
    nz1 = besi(zn, fnu, kode, n, cy, ierr);
    if (*ierr == 0 || *ierr == 3) {
        nz2 = besk(zn, fnu, kode, n, cwrk, ierr);
        if (*ierr == 0 || *ierr == 3) {
            nz = (nz1 < nz2 ? nz1 : nz2);
            ifnu = (int)fnu;
            ffnu = fnu - ifnu;
            arg = hpi * ffnu;
            csgn = CMPLX(cos(arg), sin(arg));
            i4 = (ifnu % 4) + 1;
            csgn = csgn * cip[i4-1];
            rhpi = 1.0 / hpi;
            cspn = conj(csgn) * rhpi;
            csgn = csgn * ci;
            if (kode != 2) {
                for (i = 1; i < (n+1); i++) {
                    cy[i-1] = csgn * cy[i-1] - cspn * cwrk[i-1];
                    csgn = ci * csgn;
                    cspn = -ci * cspn;
                }
                if (yy < 0.0)
                    for (i = 0; i < n; i++) { cy[i] = conj(cy[i]); }
                return nz;
            }
            r1 = cos(xx);
            r2 = sin(xx);
            ex = CMPLX(r1, r2);
            tol = fmax(d1mach[3], 1e-18);
            k1 = i1mach[14];
            k2 = i1mach[15];
            r1m5 = d1mach[4];
            k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
            elim = 2.303 * (k*r1m5 - 3.0);
            ey = 0.0;
            tay = fabs(yy + yy);
            if (tay < elim) { ey = exp(-tay); }
            cspn = ex * ey * cspn;
            nz = 0;
            rtol = 1.0 / tol;
            ascle = d1mach[0]*rtol*1e3;
            for (i = 1; i < (n+1); i++) {
                zv = cwrk[i-1];
                aa = creal(zv);
                bb = cimag(zv);
                atol = 1.0;
                if (fmax(fabs(aa), fabs(bb)) <= ascle) {
                    zv *= rtol;
                    atol = tol;
                }
                zv *= cspn;
                zv *= atol;
                zu = cy[i-1];
                aa = creal(zu);
                bb = cimag(zu);
                atol = 1.0;
                if (fmax(fabs(aa), fabs(bb)) <= ascle) {
                    zu *= rtol;
                    atol = tol;
                }
                zu *= csgn;
                zu *= atol;
                cy[i-1] = zu - zv;
                if (yy < 0.0) { cy[i-1] = conj(cy[i-1]); }
                if ((cy[i] == 0.0) && (ey == 0.0)) { nz = nz + 1; }
                csgn *= ci; 
                cspn *= -ci;
            }
        return nz;
        }
    }
    return 0;
}


int binu(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *cy,
    double rl,
    double fnul,
    double tol,
    double elim,
    double alim
) {
    double complex cw[2] = { 0. };
    double az, dfnu;
    int inw, nlast, nn, nui, nw, nz;

    nz = 0;
    az = cabs(z);
    nn = n;
    dfnu = fnu + n - 1;
    if ((az <= 2.) || (az*az*0.25 <= (dfnu + 1.0))) {
        /* GOTO 10 */
        nw = seri(z,fnu, kode, n, cy, tol, elim, alim);
        inw = abs(nw);
        nz += inw;
        nn -= inw;
        if (nn == 0) { return nz; }
        if (nw >= 0) { return nz; }
        dfnu = fnu + nn - 1;
    }
    /* GOTO 30 conditions*/
    //
    // ASYMPTOTIC EXPANSION FOR LARGE Z
    //
    if (az >= rl) {
        if ((dfnu <= 1.0) || ((dfnu > 1.0) && (az+az >= dfnu*dfnu))) {
            //
            // MILLER ALGORITHM NORMALIZED BY THE SERIES
            //
            nw = asyi(z, fnu, kode, n, cy, rl, tol, elim, alim);
            if (nw < 0) {
                nz = -1;
                if (nw == -2) {
                    nz = -2;
                }
                return nz;
            }
            return nz;
        }
    }
    /* 40 */
    if ((az < rl) || (az+az >= dfnu*dfnu)) {
        nw = mlri(z, fnu, kode, n, cy, tol);
        if (nw < 0) {
            nz = -1;
            if (nw == -2) {
                nz = -2;
            }
            return nz;
        }
        return nz;
    }
    /* 50 */
    //
    // OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
    //
    nw = uoik(z, fnu, kode, 1, nn, cy, tol, elim, alim);
    if (nw < 0) {
        nz = -1;
        if (nw == -2) { nz = -2; }
        return nz;
    }
    nz += nw;
    nn -= nw;
    if (nn == 0) { return nz; }
    dfnu = fnu + (nn -1);
    /* GOTO 110s handled here */
    if ((dfnu > fnul) || (az > fnul)) {
        nui = (int)(fnul-dfnu) + 1;
        nui = (nui > 0 ? nui : 0);
        nw = buni(z, fnu, kode, nn, cy, nui, &nlast, fnul, tol, elim, alim);
        if (nw < 0) {
            nz = -1;
            if (nw == -2) { nz = -2; }
            return nz;
        }
        nz += nw;
        if (nlast == 0) { return nz; }
        nn = nlast;
    }
    /* 60 */
    if (az <= rl) {
        /* 70 */
        nw = mlri(z, fnu, kode, n, cy, tol);
        nz = -1;
        if (nw == -2) { nz = -2; }
        return nz;
    }
    /* 80 */
    //
    // MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
    //
    //
    // OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
    //
    nw = uoik(z, fnu, kode, 2, 2, cw, tol, elim, alim);
    if (nw < 0) {
        nz = nn;
        /* 90 */
        for (int i=0; i < nn; i++) { cy[i] = 0.0; }
        return nz;
    }
    /* 100 */
    if (nw > 0) {
        return -1;
    }
    nw = wrsk(z, fnu, kode, nn, cy, cw, tol, elim, alim);
    if (nw < 0) {
        nz = -1;
        if (nw == -2) {
            nz = -2;
        }
        return nz;
    }
    return nz;
}


double complex biry(
    double complex z,
    int id,
    int kode,
    int *ierr
) {
    double complex bi, csq, s1, s2, trm1, trm2, zta, z3;
    double aa, ad, ak, alim, atrm, az, az3, bb, bk, ck, dig, dk, d1, d2,\
           elim, fid, fmr, fnu, fnul, rl, r1m5, sfac, tol, zi, zr;
    int k, k1, k2, nz;
    double complex cy[2] = { 0.0 };
    double tth = 2. / 3.;
    double c1 = 0.614926627446000735150922369;  /* 1/( 3**(1/6) Gamma(2/3)) */
    double c2 = 0.448288357353826357914823710;  /* 3**(1/6) / Gamma(1/3) */
    double coef = 0.577350269189625764509148780;  /* sqrt( 1 / 3) */
    double pi = 3.141592653589793238462643383;

    *ierr = 0;
    nz = 0;
    if ((id < 0) || (id > 1)) { *ierr= 1; }
    if ((kode < 1) || (kode > 2)) { *ierr= 1; }
    if ( *ierr != 0) { return 0.0;}
    az = cabs(z);
    tol = d1mach[3];
    fid = id;
    if (az <= 1.0) {
        //
        // POWER SERIES FOR ABS(Z) <= 1.
        //
        s1 = 1.0;
        s2 = 1.0;
        if (az < tol) {
            aa = c1 * (1.0 - fid) + fid * c2;
            return aa;
        }
        aa = az*az;
        if (aa >= tol/az) {
            trm1 = 1.0;
            trm2 = 1.0;
            atrm = 1.0;
            z3 = z*z*z;
            az3 = az * aa;
            ak = 2.0 + fid;
            bk = 3.0 - fid - fid;
            ck = 4.0 - fid;
            dk = 3.0 + fid + fid;
            d1 = ak * dk;
            d2 = bk * ck;
            ad = fmin(d1,d2);
            ak = 24.0 + 9.0*fid;
            bk = 30.0 - 9.0*fid;
            for (k = 1; k < 26; k++)
            {
                trm1 *= z3/d1;
                s1 += trm1;
                trm2 *= z3/d2;
                s2 += trm2;
                atrm *= az3 / ad;
                d1 += ak;
                d2 += bk;
                ad = fmin(d1, d2);
                if (atrm < tol*ad) { break; }
                ak += 18.0;
                bk += 18.0;
            }
        }
        if (id != 1) {
            bi = s1*c1 + z*s2*c2;
            if (kode == 1) { return bi; }
            zta = z*csqrt(z)*tth;
            aa = -fabs(creal(zta));
            bi *= exp(aa);
            return bi;
        }
        bi = s2*c2;
        if (az > tol) { bi += z*z*s1*c1/(1.0 + fid ); }
        if (kode == 1) { return bi; }
        zta = z*csqrt(z)*tth;
        aa = -fabs(creal(zta));
        bi += exp(aa);
        return bi;
    }
    //
    // CASE FOR ABS(Z) > 1.0
    //
    fnu = (1.0 + fid) / 3.0;
    //
    // SET PARAMETERS RELATED TO MACHINE CONSTANTS.
    // TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
    // ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
    // EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
    // EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
    // UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
    // RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
    // DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
    // FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
    //
    k1 = i1mach[14];
    k2 = i1mach[15];
    r1m5 = d1mach[4];
    k = (abs(k1) > abs(k2) ? abs(k2) : abs(k1) );
    elim = 2.303 * (k*r1m5 - 3.0);
    k1 = i1mach[13] - 1;
    aa = r1m5*k1;
    dig = (aa > 18.0 ? 18.0 : aa);
    aa *= 2.303;
    alim = elim + (-aa > -41.45 ? -aa : -41.45);
    rl = 1.2*dig + 3.0;
    fnul = 10.0 + 6.0*(dig - 3.0);
    // 
    // TEST FOR RANGE
    // 
    aa = 0.5 / tol;
    bb = d1mach[1] * 0.5;
    aa = fmin(aa, bb);
    aa = pow(aa, tth);
    if (az > aa) { *ierr = 4; nz = 0; return 0.0; }
    aa = sqrt(aa);
    if (az > aa) { *ierr = 3; }
    csq = csqrt(z);
    zta = z*csq*tth;
    //
    // RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
    //
    sfac = 1.0;
    zi = cimag(z);
    zr = creal(z);
    ak = cimag(zta);
    if (zr < 0.0) {
        bk = creal(zta);
        ck = -fabs(bk);
        zta = CMPLX(ck, ak);
    }
    if ((zi == 0.0) && (zr <= 0.0)) { zta = ak*I; }
    aa = creal(zta);
    if (kode != 2) {
        //
        // OVERFLOW TEST
        //
        bb = fabs(aa);
        if (bb >= alim) {
            bb += 0.5*log(az);
            sfac = tol;
            if (bb > elim) { nz = 0; *ierr = 2; return 0.0; }
        }
    }
    fmr = 0.0;
    if ((aa < 0.0) || (zr <= 0.0)) {
        fmr = pi;
        if (zi < 0.0) { fmr = -pi; }
        zta = -zta;
    }
    //
    // AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
    // KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBINU
    //
    nz = binu(zta, fnu, kode, 1, cy, rl, fnul, tol, elim, alim);
    if (nz < 0) {
        if (nz == -1) {
            nz = 0;
            *ierr = 2;
        } else {
        nz = 0;
        *ierr = 5;
        }
        return 0.0;
    }
    aa = fmr*fnu;
    z3 = sfac;
    s1 = cy[0] * CMPLX(cos(aa), sin(aa)) * z3;
    fnu = (2 - fid) / 3.0;
    nz = binu(zta, fnu, kode, 2, cy, rl, fnul, tol, elim, alim);
    cy[0] *= z3;
    cy[1] *= z3;
    //
    // BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
    //
    s2 = cy[0] * (fnu+fnu) / zta + cy[1];
    aa = fmr * (fnu - 1.0);
    s1 = (s1 + s2*CMPLX(cos(aa), sin(aa)))*coef;
    if (id != 1) {
        s1 *= csq;
        bi = s1 / sfac;
        return bi;
    }
    s1 *= z;
    bi = s1 / sfac;
    return bi;
}


int bknu(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    double tol,
    double elim,
    double alim
) {
    double complex cch, ck, coef, crsc, cs, cscl, csh, cz,\
                   f, fmu, p, pt, p1, p2, q, rz, smu, st, s1, s2, zd;
    double aa, ak, ascle, a1, a2, bb, bk, caz, dnu, dnu2, etest, fc, fhs,\
           fk, fks, g1, g2, p2i, p2m, p2r, rk, s, tm, t1, t2, xx, yy,\
           elm, xd, yd, alas, as;
    int i, iflag, inu, k, kflag, kk, koded, j, ic, inub;
    double complex cy[2];

    int kmax =30;
    double r1 = 2.;
    double pi = 3.14159265358979324;
    double rthpi = 1.25331413731550025;
    double spi = 1.90985931710274403;
    double hpi = 1.57079632679489662;
    double fpi = 1.89769999331517738;
    double tth = 2. / 3.;
    double cc[8] = {
        5.77215664901532861e-01, -4.20026350340952355e-02,
       -4.21977345555443367e-02,  7.21894324666309954e-03,
       -2.15241674114950973e-04, -2.01348547807882387e-05,
        1.13302723198169588e-06,  6.11609510448141582e-09
    };
    xx = creal(z);
    yy = cimag(z);
    caz = cabs(z);
    cscl = 1. / tol;
    crsc = tol;
    double complex css[3] = {cscl, 1., crsc};
    double complex csr[3] = {crsc, 1., cscl};
    double bry[3] = {1e3*d1mach[0]/tol, tol/(1e3*d1mach[0]), d1mach[1]};
    int nz = 0;
    iflag = 0;
    koded = kode;
    rz = 2. / z;
    inu = (int)(fnu + 0.5);
    dnu = fnu - inu;
    if (fabs(dnu) != 0.5) {
        dnu2 = 0.0;
        if (fabs(dnu) > tol) { dnu2 = dnu * dnu; }
        if (caz <= r1) {
            //
            //    SERIES FOR ABS(Z) <= R1
            //
            fc = 1.;
            smu = clog(rz);
            fmu = smu * dnu;
            csh = csinh(fmu);
            cch = ccosh(fmu);
            if (dnu != 0.0) {
                fc = dnu * pi;
                fc *= 1. / sin(fc);
                smu = csh / dnu;
            }
            a2 = 1. + dnu;
            //
            // GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
            //
            t2 = exp(-gamln(a2));
            t1 = 1. / (t2*fc);
            if (fabs(dnu) <= 0.1) {
                //
                // SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
                //
                ak = 1.;
                s = cc[0];
                for (int k = 2; k < 9; k++)
                {
                    ak *= dnu2;
                    tm = cc[k-1] * ak;
                    s += tm;
                    if (fabs(tm) < tol) { break; }
                }
                g1 = -s;
            } else {
                g1 = (t1-t2) / (dnu+dnu);
            }
            g2 = 0.5 * (t1+t2);
            f = fc*(g1*cch + smu*g2);
            pt = cexp(fmu);
            p = (0.5 / t2) * pt;
            q = (0.5 / t1) / pt;
            s1 = f;
            s2 = p;
            ak = 1.0;
            a1 = 1.0;
            ck = 1.0;
            bk = 1.0 - dnu2;
            if ((inu <= 0) && (n <= 1)) {
                //
                // GENERATE K(FNU,Z), 0.0D0  <=  FNU  <  0.5D0 AND N=1
                //
                if (caz >= tol) {
                    cz = z * z * 0.25;
                    t1 = 0.25 * caz * caz;
L30:
                    f = (f*ak + p + q) / bk;
                    p = p / (ak-dnu);
                    q = q / (ak+dnu);
                    rk = 1.0 / ak;
                    ck *= cz * rk;
                    s1 += ck * f;
                    a1 *= t1 * rk;
                    bk += ak + ak + 1.0;
                    ak += 1.0;
                    if (a1 > tol) { goto L30; }
                }
                y[0] = s1;
                if (koded == 1) { return nz; }
                y[0] = s1 * cexp(z);
                return nz;
            }
            //
            // GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
            //
            if (caz >= tol) {
                cz = z * z * 0.25;
                t1 = 0.25 * caz * caz;
L40:
                f = (f*ak + p + q) / bk;
                p *= 1.0 / (ak - dnu);
                q *= 1.0 / (ak + dnu);
                rk = 1. / ak;
                ck *= cz * rk;
                s1 += ck * f;
                s2 += ck * (p - f*ak);
                a1 *= t1 * rk;
                bk += ak + ak + 1.0;
                ak += 1.0;
                if (a1 > tol) { goto L40; }
            }
            kflag = 2;
            bk = creal(smu);
            a1 = fnu + 1.;
            ak = a1 * fabs(bk);
            if (ak > alim) { kflag = 3; }
            p2 = s2 * css[kflag-1];
            s2 = p2 * rz;
            s1 *= css[kflag-1];
            if (koded != 1) {
                f = cexp(z);
                s1 *= f;
                s2 *= f;
            }
            goto L100;
        }
    }
    //
    // IFLAG=0 MEANS NO UNDERFLOW OCCURRED
    // IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
    // KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
    // RECURSION
    //
    coef = rthpi / csqrt(z);
    kflag = 2;
    if (koded != 2) {
        if (xx > alim) { goto L200; }
        a1 = exp(-xx)*creal(css[kflag-1]);
        pt = a1*CMPLX(cos(yy), -sin(yy));
        coef *= pt;
    }

L50:
    if (fabs(dnu) == 0.5) {
        s1 = coef;
        s2 = coef;
        goto L100;
    }
//
//    MILLER ALGORITHM FOR ABS(Z) > R1
//
    ak = fabs(cos(pi*dnu));
    if (ak == 0.) {
        s1 = coef;
        s2 = coef;
        goto L100;
    }
    fhs = fabs(0.25 - dnu2);
    if (fhs == 0.) {
        s1 = coef;
        s2 = coef;
        goto L100;
    }
//
// COMPUTE R2=F(E). IF ABS(Z) >= R2, USE FORWARD RECURRENCE TO
// DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
// 12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-DIGITS(0.0_dp))=
// TOL WHERE B IS THE BASE OF THE ARITHMETIC.
//
    t1 = (i1mach[13] - 1)*d1mach[4]*(log(10)/log(2));
    t1 = fmin(fmax(t1, 12.0), 60.0);
    t2 = tth * t1 - 6.0;
    if (xx == 0.) {
        t1 = hpi;
    } else {
        t1 = fabs(atan(yy/xx));
    }
    if (t2 <= caz) {
        //
        // FORWARD RECURRENCE LOOP WHEN ABS(Z) >= R2
        //
        etest = ak / (pi*caz*tol);
        fk = 1.0;
        if (etest < 1.0) { goto L80; } 
        fks = 2.0;
        rk = caz + caz + 2.0;
        a1 = 0.0;
        a2 = 1.0;
        for (int i = 1; i < (kmax+1); i++)
        {
            ak = fhs / fks;
            bk = rk / (fk + 1.0);
            tm = a2;
            a2 = bk * a2 - ak * a1;
            a1 = tm;
            rk += 2.;
            fks += fk + fk + 2.0;
            fhs += fk + fk;
            fk += 1.0;
            tm = fabs(a2)*fk;
            if (etest < tm) {
                /* goto 160 */
                break;
            }
            if (i == kmax) {
                /* Didn't break so goes to 310 */
                return -2;
            }
        }

        /* 160 */
        fk += spi * t1 * sqrt(t2/caz);
        fhs = fabs(0.25 - dnu2);
    } else {
        //
        // COMPUTE BACKWARD INDEX K FOR ABS(Z) < R2
        //
        a2 = sqrt(caz);
        ak *= fpi / (tol*sqrt(a2));
        aa = 3.0 * t1 / (1.0 + caz);
        bb = 14.7 * t1 / (28.0 + caz);
        ak = (log(ak) + caz*cos(aa)/(1.0  + 0.008*caz)) / cos(bb);
        fk = 0.12125 * ak * ak / caz + 1.5;
    }
L80:
    //
    // BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
    //
    k = (int)fk;
    fk = (double)k;
    fks = fk * fk;
    p1 = 0.0;
    p2 = tol;
    cs = p2;
    for (i=1; i < (k+1); i++)
    {
        a1 = fks - fk;
        a2 = (fks+fk) / (a1+fhs);
        rk = 2.0 / (fk + 1.);
        t1 = (fk + xx) * rk;
        t2 = yy * rk;
        pt = p2;
        p2 = (p2 * CMPLX(t1, t2) - p1) * a2;
        p1 = pt;
        cs += p2;
        fks = a1 - fk + 1.0;
        fk -= 1.0;
    }

    //
    // COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER SCALING
    //
    tm = cabs(cs);
    pt = 1.0 / tm;
    s1 = pt * p2;
    cs = conj(cs) * pt;
    s1 *= coef * cs;
    if ((inu <= 0) && (n <= 1)) {
        zd = z;
        if (iflag == 1) { goto L190; }
        goto L130;
    }
    //
    // COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
    //
    tm = cabs(p2);
    pt = 1.0 / tm;
    p1 = pt * p1;
    p2 = conj(p2) * pt;
    pt = p1 * p2;
    s2 = s1 * (1. + (dnu+0.5 - pt)/z);
    //
    // FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
    // SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
    //
L100:
    ck = (dnu + 1.)*rz;
    if (n == 1) { inu -= 1; }
    if (inu <= 0) {
        if (n <= 1) { s1 = s2; }
        zd = z;
        if (iflag == 1) { goto L190; }
        goto L130;
    }
    inub = 1;
    if (iflag == 1) { goto L160; }
L110:
    p1 = csr[kflag-1];
    ascle = bry[kflag-1];
    for (int i = inub; i < inu+1; i++)
    {
        st = s2;
        s2 = ck*s2 + s1;
        s1 = st;
        ck += rz;
        if (kflag < 3) {
            p2 = s2*p1;
            p2m = fmax(fabs(creal(p2)), fabs(cimag(p2)));
            if (p2m > ascle) {
                kflag += 1;
                ascle = bry[kflag-1];
                s1 *= p1;
                s2 = p2;
                s1 *= css[kflag-1];
                s2 *= css[kflag-1];
                p1 = csr[kflag-1];
            }
        }
    }
    if (n == 1) { s1 = s2; }
    
L130:
    y[0] = s1 * csr[kflag-1];
    if (n == 1) { return nz; }
    y[1] = s2 * csr[kflag-1];
    if (n == 2) { return nz; }
    kk = 2;
L140:
    kk += 1;
    if (kk > n) { return nz; }
    p1 = csr[kflag-1];
    ascle = bry[kflag-1];
    for (int i = kk; i < (n+1); i++)
    {
        p2 = s2;
        s2 = ck*s2 + s1;
        s1 = p2;
        ck += rz;
        p2 = s2*p1;
        y[i-1] = p2;
        if (kflag < 3) {
            p2m = fmax(fabs(creal(p2)), fabs(cimag(p2)));
            if (p2m > ascle) {
                kflag += 1;
                ascle = bry[kflag-1];
                s1 *= p1;
                s2 = p2;
                s1 *= css[kflag-1];
                s2 *= css[kflag-1];
                p1 = csr[kflag-1];
            }
        }
    }
    return nz;
//
// IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
//
L160:
    elm = exp(-elim);
    ascle = bry[0];
    zd = z;
    xd = xx;
    yd = yy;
    ic = -1;
    j = 2;
    for (int i = 1; i < (inu+1); i++)
    {
        st = s2;
        s2 = ck*s2 + s1;
        s1 = st;
        ck += rz;
        as = cabs(s2);
        alas = log(as);
        p2r = alas - xd;
        if (p2r >= -elim) {
            p2 = -zd + clog(s2);
            p2r = creal(p2);
            p2i = cimag(p2);
            p2m = exp(p2r) / tol;
            p1 = p2m * CMPLX(cos(p2i), sin(p2i));
            if (!(uchk(p1, ascle, tol))) {
                j = 3 - j;
                cy[j-1] = p1;
                if (ic == i-1) { goto L180; }
                ic = i;
                continue;
            }
        }
        if (alas >= 0.5 * elim) {
            xd -= elim;
            s1 *= elm;
            s2 *= elm;
            zd = CMPLX(xd, yd);
        }
    }
    if (n == 1) { s1 = s2; }
    goto L190;
L180:
    kflag = 1;
    inub = i + 1;
    s2 = cy[j-1];
    j = 3 - j;
    s1 = cy[j-1];
    if (inub <= inu) { goto L110; }
    if (n == 1) { s1 = s2; }
    goto L130;
L190:
    y[0] = s1;
    if (n != 1) { y[1] = s2; }
    ascle = bry[0];
    nz = kscl(zd, fnu, n, &y[0], rz, &ascle, tol, elim);
    inu = n - nz;
    if (inu < 0) { return nz; }
    kk = nz + 1;
    s1 = y[kk-1];
    y[kk-1] = s1 * csr[0];
    if (inu == 1) { return nz; }
    kk = nz + 2;
    s2 = y[kk-1];
    y[kk-1] = s2 * csr[0];
    if (inu == 2) { return nz; }
    t2 = fnu + (kk-1);
    ck = t2 * rz;
    kflag = 1;
    goto L140;
L200:
    koded = 2;
    iflag = 1;
    kflag = 2;
    goto L50;
}


int buni(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    int nui,
    int *nlast,
    double fnul,
    double tol,
    double elim,
    double alim
) {
    double complex cscl, cscr, rz, st, s1, s2;
    double ax, ay, dfnu, fnui, gnu, xx, yy, ascle, str, sti, stm;
    int i, iflag, iform, k, nl, nw, nz;
    double complex cy[2] = { 0.0 };
    double bry[3] = { 0.0 };

    nz = 0;
    xx = creal(z);
    yy = cimag(z);
    ax = fabs(xx) + sqrt(3.);
    ay = fabs(yy);
    iform = 1;
    if (ay > ax) { iform = 2; }
    if (nui == 0) {
        if (iform != 2) {
            uni1(z, fnu, kode, n, y, &nw, nlast, fnul, tol, elim, alim);
        } else {
            uni2(z, fnu, kode, n, y, &nw, nlast, fnul, tol, elim, alim);
        }
        if (nw < 0) {
            nz = -1;
            if (nw == -2) { nz = -2; }
            return nz;
        }
        return nw;
    }

    fnui = nui;
    dfnu = fnu + (n - 1);
    gnu = dfnu + fnui;
    if (iform != 2) {
        //
        // ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
        // -PI/3 <= ARG(Z) <= PI/3
        //
        uni1(z, gnu, kode, 2, cy, &nw, nlast, fnul, tol, elim, alim);
    } else {
        uni2(z, gnu, kode, 2, cy, &nw, nlast, fnul, tol, elim, alim);
    }
    if (nw >= 0) {
        if (nw != 0) { *nlast = n; return nz; }
        ay = cabs(cy[0]);
        //
        // SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
        //
        bry[0] = 1e3*d1mach[0] / tol;
        bry[1] = tol / 1e3*d1mach[0];
        bry[2] = bry[1];
        iflag = 2;
        ascle = bry[1];
        ax = 1.0;
        cscl = ax;
        if (ay <= bry[0]) {
            iflag = 1;
            ascle = bry[0];
            ax = 1.0 / tol;
            cscl = ax;
        } else {
            if (ay >= bry[1]) {
                iflag = 3;
                ascle = bry[2];
                ax = tol;
                cscl = ax;

            }
        }
        ay = 1.0 / ax;
        cscr = ay;
        s1 = cy[1] * cscl;
        s2 = cy[0] * cscl;
        rz = 2.0 / z;
        for (i = 1; i < nui; i++)
        {
            st = s2;
            s2 = (dfnu +fnui)*rz*s2 + s1;
            s1 = st;
            fnui -= 1.0;
            if (iflag < 3) {
                st = s2 * cscr;
                str = fabs(creal(st));
                sti = fabs(cimag(st));
                stm = fmax(str, sti);
                if (stm > ascle) {
                    iflag += 1;
                    ascle = bry[iflag-1];
                    s1 *= cscr;
                    s2 = st;
                    ax *= tol;
                    ay = 1.0 / ax;
                    cscl = ax;
                    cscr = ay;
                    s1 *= cscl;
                    s2 *= cscl;
                }
            }
        }
        y[n-1] = s2*cscr;
        if (n == 1) { return nz; }
        nl = n-1;
        fnui = nl;
        k = nl;
        for (i = 0; i < (nl+1); i++)
        {
            st = s2;
            s2 = (fnu + fnui)*rz*s2 + s1;
            s1 = st;
            st = s2 * cscr;
            y[k-1] = st;
            fnui -= 1.0;
            k -= 1;
            if (iflag < 3) {
                st = s2 * cscr;
                str = fabs(creal(st));
                sti = fabs(cimag(st));
                stm = fmax(str, sti);
                if (stm > ascle) {
                    iflag += 1;
                    ascle = bry[iflag-1];
                    s1 *= cscr;
                    s2 = st;
                    ax *= tol;
                    ay = 1.0 / ax;
                    cscl = ax;
                    cscr = ay;
                    s1 *= cscl;
                    s2 *= cscl;
                }
            }
        }
        return nz;
    }
    nz = -1;
    if (nw == -2) { nz = -2; }
    return nz;
}


int bunk(
    double complex z,
    double fnu,
    int kode,
    int mr,
    int n,
    double complex *y,
    double tol,
    double elim,
    double alim
) {
    double ax, ay;
    
    int nz = 0;
    ax = fabs(creal(z)) * 1.7321;
    ay = fabs(cimag(z));
    
    if (ay <= ax) {
        //
        // Asymptotic expansion for K(FNU,Z) for large FNU applied in
        // -PI/3 <= ARG(Z) <= PI/3
        //
        nz = unk1(z, fnu, kode, mr, n, y, tol, elim, alim);
    } else {
        //
        // Asymptotic expansion for H(2, FNU, Z*EXP(M*HPI)) for large FNU
        // applied in PI/3 < ABS(ARG(Z)) <= PI/2 where M = +I or -I and HPI = PI/2
        //
        nz = unk2(z, fnu, kode, mr, n, y, tol, elim, alim);
    }
    return nz;
}


double gamln(double z) {
    int i1m, mz;
    double fln, fz, rln, s, tlg, trm, tst, t1, wdtol, zdmy, zinc, zm, zmin, zp, zsq;
    const double con = 1.83787706640934548;  /* LN(2*PI) */
    int nz = 0;
    if (z > 0.0) {
        if (z <= 101.0) {
            nz = (int)z;
            fz = z - nz;
            if (fz <= 0.0) {
                if (nz <= 100) {
                    return dgamln_gln[nz-1];
                }
            }
        }
        wdtol = fmax(d1mach[3], 1e-18);
        i1m = i1mach[13];
        rln = d1mach[4]*i1m;
        fln = fmax(fmin(rln, 20.), 3.0) - 3.0;
        zm = 1.8 + 0.3875*fln;
        mz = ((int)zm) + 1;
        zmin = mz;
        zdmy = z;
        zinc = 0.0;
        if (z < zmin){
            zinc = zmin - nz;
            zdmy = z + zinc;
        }
        zp = 1. / zdmy;
        t1 = dgamln_cf[0]*zp;
        s = t1;
        if (zp >= wdtol) {
            zsq = zp*zp;
            tst = t1*wdtol;
            for (int i = 2; i < 23; i++)
            {
                zp *= zsq;
                trm = dgamln_cf[i-1] * zp;
                if (fabs(trm) < tst) { break; }
                s += trm;
            }
        }

        if (zinc == 0.) {
            tlg = log(z);
            return z*(tlg-1.0) + 0.5*(con - tlg) + s;
        }
        zp = 1.0;
        nz = (int)zinc;
        for (int i = 0; i < nz; i++)
        {
            zp *= (z + i);
        }
        tlg = log(zdmy);
        return zdmy*(tlg-1.0) - log(zp) + 0.5*(con-tlg) + s;
    }
    // Zero or negative argument
    return NAN;
}


int mlri(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    double tol
) {
    double complex ck, cnorm, pt, p1, p2, rz, sum;
    double ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, rho,\
           rho2, scle, tfnf, tst, x;
    int i, iaz, ifnu, inu, itime, k, kk, km, m, nz;
    scle = d1mach[0] / tol;
    nz = 0;
    az = cabs(z);
    x = creal(z);
    iaz = (int)az;
    ifnu = (int)fnu;
    inu = ifnu + n - 1;
    at = iaz + 1;
    ck = at / z;
    rz = 2. / z;
    p1 = 0.;
    p2 = 1.;
    ack = (at + 1.0) / az;
    rho = ack + sqrt(ack*ack - 1.);
    rho2 = rho * rho;
    tst = (rho2 + rho2) / ((rho2 - 1.0)*(rho - 1.0));
    tst /= tol;
    //
    // COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
    //
    ak = at;
    i = 1;
    for (i = 1; i < 81; i++ )
    {
        pt = p2;
        p2 = p1 - ck * p2;
        p1 = pt;
        ck += rz;
        ap = cabs(p2);
        if (ap > tst*ak*ak) { break; }
        ak += 1.0;
        if (i == 80) {
            /* Exhausted loop without break */
            return -2;
        }
    }
    i += 1;
    k = 0;
    if (inu >= iaz) {
    //
    // COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
    //
        p1 = 0.0;
        p2 = 1.0;
        at = inu + 1;
        ck = at / z;
        ack = at / az;
        tst = sqrt(ack / tol);
        itime = 1;
        k = 1;
        for (k = 1; k < 81; k++ )
        {
            pt = p2;
            p2 = p1 - ck * p2;
            p1 = pt;
            ck += rz;
            ap = cabs(p2);
            if (ap >= tst) {
                if (itime == 2) { break; }
                ack = cabs(ck);
                flam = ack + sqrt(ack*ack - 1.0);
                fkap = ap / cabs(p1);
                rho = fmin(flam, fkap);
                tst *= sqrt(rho / (rho*rho - 1.0));
                itime = 2;
            }
            if (i == 80) {
                /* Exhausted loop without break */
                return -2;
            }
        }
    }
    //
    // BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
    //
    k += 1;
    kk = fmax(i+iaz, k+inu);
    fkk = kk;
    p1 = 0.0;
    //
    // SCALE P2 AND SUM BY SCLE
    //
    p2 = scle;
    fnf = fnu - ifnu;
    tfnf = fnf + fnf;
    bk = gamln(fkk + tfnf + 1.0) - gamln(fkk + 1.0) - gamln(tfnf + 1.0);
    bk = exp(bk);
    sum = 0.;
    km = kk - inu;
    for (i = 1; i < (km+1); i++)
    {
        pt = p2;
        p2 = p1 + (fkk + fnf)*rz*p2;
        p1 = pt;
        ak = 1. - tfnf / (fkk+tfnf);
        ack = bk*ak;
        sum += (ack + bk)*p1;
        bk = ack;
        fkk -= 1.;
    }
    y[n-1] = p2;
    if (n != 1) {
        for (i = 2; i < (n+1); i++)
        {
            pt = p2;
            p2 = p1 + (fkk + fnf)*rz*p2;
            p1 = pt;
            ak = 1. - tfnf / (fkk+tfnf);
            ack = bk*ak;
            sum += (ack + bk)*p1;
            bk = ack;
            fkk -= 1.;
            m = n - i + 1;
            y[m-1] = p2;
        }
    }
    if (ifnu > 0) {
        for (i = 1; i < (ifnu+1); i++)
        {
            pt = p2;
            p2 = p1 + (fkk + fnf)*rz*p2;
            p1 = pt;
            ak = 1. - tfnf / (fkk+tfnf);
            ack = bk*ak;
            sum += (ack + bk)*p1;
            bk = ack;
            fkk -= 1.;
        }
    }
    pt = z;
    if (kode == 2) { pt -= x; }
    p1 = -fnf * clog(rz) + pt;
    ap = gamln(1. + fnf);
    pt = p1 - ap;
    //
    // THE DIVISION EXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
    // IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
    //
    p2 += sum;
    ap = cabs(p2);
    p1 = 1. / ap;
    ck = cexp(pt) * p1;
    pt = conj(p2)*p1;
    cnorm = ck * pt;
    for (int i = 0; i < n; i++) { y[i] *= cnorm; }
    return nz;
}


int kscl(
    double complex zr,
    double fnu,
    int n,
    double complex *y,
    double complex rz,
    double *ascle,
    double tol,
    double elim
) {
    double complex cy[2] = { 0. };
    double as, acs, alas, fn, zri, xx;
    double complex s1, s2, cs, ck, zd;
    int nz = 0;
    int ic = 0;
    int nn = (3 <= n+1? 3 : n + 1);
    int kk = 0;
    double elm = exp(-elim);

    for (int i = 0; i < nn; i++)
    {
        s1 = y[i];
        cy[i] = s1;
        as = cabs(s1);
        acs = -creal(zr) + log(as);
        nz += 1;
        y[i] = 0.;
        if (acs < -elim) {
            continue;
        }
        cs = -zr + clog(s1);
        cs = (exp(creal(cs))/tol)*(cos(cimag(cs)) + sin(cimag(cs)*I));
        if (!uchk(cs, *ascle, tol)) {
            y[i] = cs;
            nz -= 1;
            ic = i;
        }
    }
    if (n == 1) {
        return nz;
    }
    if (ic <= 1) {
        y[0] = 0.;
        nz = 2;
    }
    if (n == 2) {
        return nz;
    }
    fn = fnu + 1.;
    ck = fn*rz;
    s1 = cy[0];
    s2 = cy[1];
    zri = cimag(zr);
    zd = zr;
    for (int i = 3; i < (n+1); i++)
    {
        kk = i;
        cs = s2;
        s2 *= ck;
        s2 += s1;
        s1 = cs;
        ck += rz;
        as = cabs(s2);
        alas = log(as);
        acs = alas - xx;
        nz += 1;
        y[i-1] = 0.;
        if (acs >= -elim) {
            cs = clog(s2);
            cs -= zd;
            cs = (exp(creal(cs))/tol)*(cos(cimag(cs)) + sin(cimag(cs)*I));
            if (!uchk(cs, *ascle, tol)) {
                y[i-1] = cs;
                nz -= 1;
                if (ic == kk-1) {
                    nz = kk - 2;
                    for (int i = 0; i < nz; i++)
                    {
                        y[i] = 0.;
                    }
                    return nz;
                }
                ic = kk;
                continue;
            }
            if (alas >= 0.5*elim){
                xx -= elim;
                s1 *= elm;
                s2 *= elm;
                zd = xx + zri*I;
            }
        }
    }
    nz = n;
    if (ic == n) {
        nz = n-1;
    } else {
        nz = kk - 2;
    }
    for (int i = 0; i < nz; i++)
    {
        y[i] = 0.;
    }
    return nz;
}


void rati(
    double complex z,
    double fnu,
    int n,
    double complex *cy,
    double tol
    ) {
    double complex cdfnu, pt, p1, p2, rz, t1;
    double ak, amagz, ap1, ap2, arg, az, dfnu, fdnu, flam, fnup, rap1, rho, test, test1;
    int i, id, idnu, inu, itime, k, kk, magz;
    
    az = cabs(z);
    inu = (int)fnu;
    idnu = inu + n - 1;
    fdnu = idnu;
    magz = az;
    amagz = magz + 1;
    fnup = fmax(amagz, fdnu);
    id = idnu - magz - 1;
    itime = 1;
    k = 1;
    rz = 2.0 / z;
    t1 = fnup * rz;
    p2 = -t1;
    p1 = 1.0;
    t1 += rz;
    if (id > 0) {
        id = 0;
    }
    ap2 = cabs(p2);
    ap1 = cabs(p1);
    
    //
    // THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
    // GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT P2
    // VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR PREMATURELY.
    //    
    arg = (ap2 + ap2) / (ap1 * tol);
    test1 = sqrt(arg);
    test = test1;
    rap1 = 1.0 / ap1;
    p1 *= rap1;
    p2 *= rap1;
    ap2 *= rap1;
    
    while (1) {
        k += 1;
        ap1 = ap2;
        pt = p2;
        p2 = p1 - t1*p2;
        p1 = pt;
        t1 += rz;
        ap2 = cabs(p2);
        if (ap1 > test) { break; }
        if (itime != 2) {
            ak = cabs(t1)*0.5;
            flam = ak + sqrt(ak*ak - 1.0);
            rho = fmin(ap2/ap1, flam);
            test = test1*sqrt(rho / (rho*rho - 1.0));
            itime = 2;
        }
    }
    kk = k + 1 - id;
    ak = kk;
    dfnu = fnu + n - 1;
    cdfnu = dfnu;
    t1 = ak;
    p1 = 1.0 / ap2;
    p2 = 0.0;
    
    for (i = 1; i < (kk+1); i++) {
        pt = p1;
        p1 = rz*(cdfnu+t1)*p1 + p2;
        p2 = pt;
        t1 -= 1.0;
    }
    
    if (p1 == 0.) {
        p1 = CMPLX(p1, p1);
    }
    
    cy[n-1] = p2 / p1;
    if (n == 1) { return; }
    k = n - 1;
    ak = k;
    t1 = ak;
    cdfnu = fnu*rz;
    
    for (i = 2; i < (n+1); i++) {
        pt = cdfnu + t1*rz*cy[k];
        if (pt == 0.0) {
            pt = CMPLX(tol, tol);
        }
        cy[k-1] = 1.0 / pt;
        t1 -= 1.0;
        k -= 1;
    }
    return;
}


int seri(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    double tol,
    double elim,
    double alim
) {
    double complex ak1, ck, coef, crsc, cz, hz, rz, s1, s2, w[2];
    double aa, acz, ak, arm, ascle, atol, az, dfnu, fnup, rak1,\
           rs, rtr1, s, ss, x;
    int ib, iflag, il, k, l, m, nn;

    int nz = 0;
    az = cabs(z);
    if (az == 0.0) {
        y[0] = 0.0;
        if (fnu == 0.) { y[0] = 1.0; }
        if (n == 1) { return nz; }
        for (int i = 1; i < n; i++) { y[i] = 0.0; }
        return nz;
    }
    x = creal(z);
    arm = 1e3*d1mach[0];
    rtr1 = sqrt(arm);
    crsc = 1.0;
    iflag = 0;
    if (az >= arm) {
        hz = z*0.5;
        cz = 0.;
        if (az > rtr1) { cz = hz*hz; }
        acz = cabs(cz);
        nn = n;
        ck = clog(hz);
L10:
        dfnu = fnu + (nn-1);
        fnup = dfnu + 1.0;
        //
        // UNDERFLOW TEST
        //
        ak1 = ck * dfnu;
        ak = gamln(fnup);
        ak1 -= ak;
        if (kode == 2) { ak1 -= x; }
        rak1 = creal(ak1);
        if (rak1 > -elim) { goto L30; }
L20:
        nz += 1;
        y[nn - 1] = 0.0;
        //
        // RETURN WITH NZ < 0 IF ABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
        // THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
        //
        if (acz > dfnu) { return -nz; }
        nn -= 1;
        if (nn == 0) { return nz; }
        goto L10;
L30:
        if (rak1 <= -alim) {
            iflag = 1;
            ss = 1.0 / tol;
            crsc = tol;
            ascle = arm * ss;
        }
        ak = cimag(ak1);
        aa = exp(rak1);
        if (iflag == 1) { aa *= ss; }
        coef = aa * (cos(ak) + sin(ak)*I);
        atol = tol * acz / fnup;
        il = (nn > 2 ? 2 : nn);
        for (int i = 1; i < (il +1); i++)
        {
            dfnu = fnu + (nn-i);
            fnup = dfnu + 1.0;
            s1 = 1.0;
            if (acz >= tol*fnup) {
                ak1 = 1.0;
                ak = fnup + 2.0;
                s = fnup;
                aa = 2.0;
L40:
                rs = 1.0 / s;
                ak1 *= cz * rs;
                s1 += ak1;
                s += ak;
                ak += 2.0;
                aa *= acz * rs;
                if (aa > atol) { goto L40; }
            }
            m = nn - i + 1;
            s2 = s1 * coef;
            w[i-1] = s2;
            if (iflag != 0) {
                if (uchk(s2, ascle, tol)) { goto L20; }
            }
            y[m-1] = s2 * crsc;
            if (i != il) { coef *= dfnu / hz; }
        }
        if (nn <= 2) { return nz; }
        k = nn - 2;
        ak = k;
        rz = 2.0 / z;
        if (iflag == 1) { goto L80; }
        ib = 3;
L60:
        for (int i = ib; i < (nn+1); i++)
        {
            y[k-1] = (ak+fnu)*rz*y[k] + y[k+1];
            ak -= 1.0;
            k -= 1;
        }
        return nz;
L80:
        s1 = w[0];
        s2 = w[1];
        l = 3;
        for (int l = 3; l < (nn+1); l++)
        {
            ck = s2;
            s2 = s1 + (ak+fnu)*rz*s2;
            s1= ck;
            ck = s2*crsc;
            y[k-1] = ck;
            ak -= 1.0;
            k -= 1;
            if (cabs(ck) > ascle) { goto L100; }
        }
        return nz;
L100:
        ib = l+1;
        if (ib > nn) { return nz; }
        goto L60;
    }
    nz = n;
    if (fnu == 0.0) { nz -= 1; }
    y[0] = 0.0;
    if (fnu == 0.) { y[0] = 1.0; }
    if (n == 1) { return nz; }
    for (int i = 1; i < n; i++) { y[i] = 0.0; }
    return nz;
}


int s1s2(
    double complex zr,
    double complex *s1,
    double complex *s2,
    double ascle,
    double alim,
    int *iuf
) {
    double complex c1, s1d;
    double aa, aln, as1, as2, xx;
    int nz = 0;
    as1 = cabs(*s1);
    as2 = cabs(*s2);
    aa = creal(*s1);
    aln = cimag(*s1);

    if ((aa != 0.) || (aln != 0.)) {
        if (as1 != 0.){
            xx = creal(zr);
            aln = -xx - xx + log(as1);
            s1d = *s1;
            *s1 = 0.;
            as1 = 0.;
            if (aln >= -alim) {
                c1 = clog(s1d) - zr - zr;
                *s1 = cexp(c1);
                as1 = cabs(*s1);
                *iuf += 1;
            }
        }
    }
    aa = fmax(as1, as2);
    if (aa > ascle) {
        return nz;
    }
    *s1 = 0.;
    *s2 = 0.;
    *iuf = 0;
    return 1;
}


int uchk(
    double complex y,
    double ascle,
    double tol
) {
    double yr = fabs(creal(y));
    double yi = fabs(cimag(y));
    double ss = fmax(yr, yi);
    double st = fmin(yr, yi);
    if (st > ascle) {
        return 0;
    } else {
        st /= tol;
        if (ss < st) {
            return 1;
        } else {
            return 0;
        }
    }
}


void unhj(
    double complex z,
    double fnu,
    int ipmtr,
    double tol,
    double complex *phi,
    double complex *arg,
    double complex *zeta1,
    double complex *zeta2,
    double complex *asum,
    double complex *bsum
) {
    double complex cfnu, przth, ptfn, rtzta, rzth, suma, sumb;
    double complex tfn, t2, w, w2, za, zb, zc, zeta, zth;
    double ang, atol, aw2, azth, btol, fn13, fn23, pp, rfn13;
    double rfnu, rfnu2, wi, wr, zci, zcr, zetai, zetar, zthi;
    double zthr, asumr, asumi, bsumr, bsumi, test, tstr, tsti, ac;
    double ex1 = 1./3.;
    double ex2 = 2./3.;
    double hpi = 1.57079632679489662;
    double pi = 3.14159265358979324;
    double thpi = 4.71238898038468986;
    int ias, ibs, j, ju, k, kmax, kp1, ks, l, lrp1, l1, l2, m;
    /* array vars */
    double complex cr[14] = { 0. };
    double complex dr[14] = { 0. };
    double complex up[14] = { 0. };
    double complex p[30] = { 0. };
    double ap[30] = { 0. };

    rfnu = 1. / fnu;
    tstr = creal(z);
    tsti = cimag(z);
    test = d1mach[0] * 1e3;
    ac = fnu*test;
    if ((fabs(tstr) <= ac) && (fabs(tsti) <= ac)) {
        ac = 2.*fabs(log(test)) + fnu;
        *zeta1 = ac;
        *zeta2 = fnu;
        *phi = 1.;
        *arg = 1.;
        return;
    }
    zb = z*rfnu;
    rfnu2 = rfnu*rfnu;

    fn13 = pow(fnu, ex1);
    fn23 = fn13 * fn13;
    rfn13 = 1./fn13;
    w2 = 1. - zb*zb;
    aw2 = cabs(w2);
    if (aw2 <= 0.25) {
        k = 1;
        p[0] = 1.;
        suma = zunhj_gama[0];
        ap[0] = 1.;
        if (aw2 >= tol) {
            for (int k = 2; k < 31; k++)
            {
                p[k-1] = p[k-2]*w2;
                suma += p[k-1]*zunhj_gama[k-1];
                ap[k-1] = ap[k-2]*aw2;
                if (ap[k-1] < tol) { break; }
            }
        }

        kmax = k;
        zeta = w2*suma;
        *arg = zeta*fn23;
        za = sqrt(suma);
        *zeta2 = sqrt(w2)*fnu;
        *zeta1 = *zeta2 * (1. + zeta*za*ex2);
        za = za + za;
        *phi = sqrt(za)*rfn13;
        if (ipmtr == 1) { return; }

        sumb = 0.;
        for (k = 1; k < kmax+1; k++) {
            sumb += p[k-1]*zunhj_beta[k-1];
        }
        *asum = 0.;
        *bsum = sumb;
        l1 = 0;
        l2 = 30;
        btol = tol * (fabs(creal(*bsum)) + fabs(cimag(*bsum)));
        atol = tol;
        pp = 1.;
        ias = 0;
        ibs = 0;
        if (rfnu2 < tol) {
            *asum += 1.;
            *bsum *= rfnu*rfn13;
            return;
        }
        for (int is = 2; is < 8; is++)
        {
            atol /= rfnu2;
            pp *= rfnu2;
            if (ias != 1) {
                suma = 0.;
                for (int k = 1; k < (kmax+1); k++)
                {
                    m = l1 + k;
                    suma += p[k-1]*zunhj_alfa[m-1];
                    if (ap[k-1] < atol) { return; }
                }
                *asum += suma*pp;
                if (pp < tol) { ias = 1; }
            }
            if (ibs != 1) {
                sumb = 0.;
                for (int k = 1; k < (kmax+1); k++)
                {
                    m = l2 + k;
                    sumb += p[k-1]*zunhj_beta[m-1];
                    if (ap[k-1] < atol) { return; }
                }
                *bsum += sumb*pp;
                if (pp < btol) { ibs = 1; }
            }
            if ((ias == 1) && (ibs == 1)) { return; }
            l1 += 30;
            l2 += 30;
        }
        *asum += 1.;
        *bsum *= rfnu*rfn13;
        return;
    } else {
        w = csqrt(w2);
        wr = creal(w);
        wi = cimag(w);
        if (wr < 0) { wr = 0.;}
        if (wi < 0) { wi = 0.;}
        w = wr + wi*I;
        za = (1. + w) / zb;
        zc = clog(za);
        zcr = creal(zc);
        zci = cimag(zc);
        if (zci < 0) { zci = 0.;}
        if (zci > hpi) { zci = hpi;}
        if (zcr < 0) { zcr = 0.;}
        zc = zcr + zci*I;
        zth = (zc-w)*1.5;
        cfnu = fnu;
        *zeta1 = zc*cfnu;
        *zeta2 = w*cfnu;
        azth = cabs(zth);
        zthr = creal(zth);
        zthi = cimag(zth);
        ang = thpi;
        if ((zthr < 0.) || (zthi >= 0.)) {
            ang = hpi;
            if (zthr != 0.) {
                ang = atan(zthi/zthr);
                if (zthr < 0.) { ang += pi; }
            }
        }
        pp = pow(azth, ex2);
        ang *= ex2;
        zetar = pp * cos(ang);
        zetai = pp * sin(ang);
        if (zetai < 0.) { zetai = 0.; }
        zeta = zetar + zetai*I;
        *arg = zeta*fn23;
        rtzta = zth / zeta;
        za = rtzta / w;
        *phi = csqrt(za + za) * rfn13;
        if (ipmtr == 1) { return; }
        tfn = rfnu / w;
        rzth = rfnu / zth;
        zc = rzth * zunhj_ar[1];
        t2 = 1. / w2;
        up[1] = (t2*zunhj_c[1] + zunhj_c[2])*tfn;
        *bsum = up[1] + zc;
        *asum = 0.;

        if (rfnu < tol) {
            *asum += 1.;
            *bsum *= -rfn13 / rtzta;
            return;
        }

        przth = rzth;
        ptfn = tfn;
        up[0] = 1.;
        pp = 1.;
        bsumr = creal(*bsum);
        bsumi = cimag(*bsum);
        btol = tol * (fabs(bsumr) + fabs(bsumi));
        ks = 0;
        kp1 = 2;
        l = 3;
        ias = 0;
        ibs = 0;

        for (int lr = 2; lr < 13; lr += 2)
        {
            lrp1 = lr + 1;
            for (int k = lr; k < (lrp1+1); k++)
            {
                ks += 1;
                kp1 += 1;
                l += 1;
                za = zunhj_c[l-1];
                for (int j = 2; j < (kp1+1); j++)
                {
                    l += 1;
                    za = za*t2 + zunhj_c[l-1];
                }
                ptfn *= tfn;
                up[kp1-1] = ptfn*za;
                cr[ks-1] = przth*zunhj_br[ks];
                przth *= rzth;
                dr[ks-1] = przth*zunhj_ar[ks+1];
            }
            pp *= rfnu2;
            if (ias != 1) {
                suma = up[lr];
                ju = lrp1;
                for (int jr = 1; j < lrp1; j++)
                {
                    ju -= 1;
                    suma += cr[jr-1] * up[ju-1];
                }
                *asum += suma;
                asumr = creal(*asum);
                asumi = cimag(*asum);
                test = fabs(asumr) + fabs(asumi);
                if ((pp < tol) && (test < tol)) { ias = 1; }
            }
            if (ibs != 1) {
                sumb = up[lr+1] + up[lr]*zc;
                ju = lrp1;
                for (int jr = 1; j < lrp1; j++)
                {
                    ju -= 1;
                    sumb += dr[jr-1] * up[ju-1];
                }
                *bsum += sumb;
                bsumr = creal(*bsum);
                bsumi = cimag(*bsum);
                test = fabs(bsumr) + fabs(bsumi);
                if ((pp < tol) && (test < tol)) { ibs = 1; }
            }
            if ((ias == 1) && (ibs == 1)) { break; }
        }
        *asum += 1.;
        *bsum *= -rfn13 / rtzta;
        return;
    }
}


void uni1(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    int *nz,
    int *nlast,
    double fnul,
    double tol,
    double elim,
    double alim
) {
    double complex cfn, crsc, cscl, c1, c2, phi, rz, sum, s1, s2, zeta1, zeta2;
    double aphi, ascle, c2i, c2m, c2r, fn, rs1, yy;
    int iflag, init, k, m, nd, nn, resetfor = 0;
    double complex cwrk[16] = { 0. };
    double complex cy[2] = { 0. };
    nz = 0;
    nd = n;
    nlast = 0;
    cscl = 1.;
    crsc = tol;
    double complex css[3] = {cscl, 1., crsc};
    double complex csr[3] = {crsc, 1., cscl};
    double bry[3] = {1e3*d1mach[0]/tol, 0., 0.};
    fn = fmax(fnu, 1.);
    init = 0;
    unik(z, fn, 1, 1, tol, &init, &phi, &zeta1, &zeta2, &sum, &cwrk[0]);
    if (kode != 1) {
        cfn = fn;
        s1 = -zeta1 + cfn * (cfn/(z+zeta2));
    } else {
        s1 = zeta2 - zeta1;
    }
    rs1 = creal(s1);
    if (fabs(rs1) > elim) {
        if (rs1 > 0) {
            *nz = -1;
            return;
        }
        *nz = n;
        for (int i = 0; i < n; i++) { y[i] = 0.; }
    }

    while (1) {
        if (resetfor == 1) { resetfor = 0; }
        nn = (nd > 2 ? 2 : nd);
        for (int i = 1; i < (nn+1); i++)
        {
            fn = fnu + (nd-i);
            init = 0;
            unik(z, fn, 1, 0, tol, &init, &phi, &zeta1, &zeta2, &sum, &cwrk[0]);
            if (kode != 1) {
                cfn = fn;
                yy = cimag(z);
                s1 = -zeta1 + cfn*(cfn/(z/zeta2)) + yy*I;
            } else {
                s1 = zeta2 - zeta1;
            }
            //
            // TEST FOR UNDERFLOW AND OVERFLOW
            //
            rs1 = creal(s1);
            if (fabs(rs1) > elim) {
                if (rs1 <= 0.) {
                    y[nd-1] = 0.;
                    nz += 1;
                    nd -= 1;
                    if (nd == 0) { return; }
                    int nuf = uoik(z, fnu, kode, 1, nd, &y[0], tol, elim, alim);
                    if (nuf >= 0) {
                        nd -= nuf;
                        nz += nuf;
                        if (nd == 0) { return; }
                        fn = fnu + (nd -1);
                        /* Resetting for loop (GOTO 30) */
                        if (fn >= fnul) {
                            resetfor = 1;
                            break;
                        }
                        *nlast = nd;
                        return;
                    }
                }
            }
            if (i == 1) { iflag = 2; }
            if (fabs(rs1) >= alim) {
                //
                // REFINE TEST AND SCALE
                //
                aphi = cabs(phi);
                rs1 += log(aphi);

                /* another go to 110 */
                if (fabs(rs1) > elim) {
                    if (rs1 <= 0.) {
                        y[nd-1] = 0.;
                        nz += 1;
                        nd -= 1;
                        if (nd == 0) { return; }
                        int nuf = uoik(z, fnu, kode, 1, nd, &y[0], tol, elim, alim);
                        if (nuf >= 0) {
                            nd -= nuf;
                            nz += nuf;
                            if (nd == 0) { return; }
                            fn = fnu + (nd -1);
                            /* Resetting for loop (GOTO 30) */
                            if (fn >= fnul) {
                                resetfor = 1;
                                break;
                            }
                            *nlast = nd;
                            return;
                        }
                    }
                }
                if (i == 1) { iflag = 1; }
                if ((rs1 >= 0.) && (i == 1)) {
                    iflag = 3;
                }
            }
            //
            // SCALE S1 IF ABS(S1) < ASCLE
            //
            s2 = phi * sum;
            c2r = creal(s1);
            c2i = cimag(s1);
            c2m = exp(c2r)*creal(css[iflag-1]);
            s1 = c2m * (cos(c2i) + sin(c2i)*I);
            s2 *= s1;
            if (iflag == 1) {
                if (!(uchk(s2, bry[0], tol))) {
                    /* another go to 110 */
                    if (rs1 <= 0.) {
                        y[nd-1] = 0.;
                        nz += 1;
                        nd -= 1;
                        if (nd == 0) { return; }
                        int nuf = uoik(z, fnu, kode, 1, nd, &y[0], tol, elim, alim);
                        if (nuf >= 0) {
                            nd -= nuf;
                            nz += nuf;
                            if (nd == 0) { return; }
                            fn = fnu + (nd -1);
                            /* Resetting for loop (GOTO 30) */
                            if (fn >= fnul) {
                                resetfor = 1;
                                break;
                            }
                            *nlast = nd;
                            return;
                        }
                    }
                }
            }
            m = nd - i + 1;
            cy[i-1] = s2;
            y[m-1] = s2*csr[iflag-1];
        }
        /* Get out of while loop */
        if (resetfor == 0) { break; }
    }    
    if (nd <= 2) { return; }

    rz = 2. / z;
    bry[1] = 1. / bry[0];
    bry[2] = d1mach[1];
    s1 = cy[0];
    s2 = cy[1];
    c1 = csr[iflag-1];
    ascle = bry[iflag-1];
    k = nd - 2;
    fn = k;
    for (int i = 3; i < nd+1; i++)
    {
        c2 = s2;
        s2 = s1 + (fnu+fn)*rz*s2;
        s1 = c2;
        c2 = s2*c1;
        y[k-1] = c2;
        k -= 1;
        fn -= 1.;
        if (iflag < 3) {
            c2r = fabs(creal(c2));
            c2i = fabs(cimag(c2));
            c2m = fmax(c2r, c2i);
            if (c2m > ascle) {
                iflag += 1;
                ascle = bry[iflag-1];
                s1 *= c1;
                s2 = c2;
                s1 *= css[iflag-1];
                s2 *= css[iflag-1];
                c1 = csr[iflag-1];
            }
        }
    }
    return;
}


void uni2(
    double complex z,
    double fnu,
    int kode,
    int n,
    double complex *y,
    int *nz,
    int *nlast,
    double fnul,
    double tol,
    double elim,
    double alim
) {
    double complex ai, arg, asum, bsum, cfn, cid, crsc, cscl, c1, c2, dai, phi, rz,\
                   s1, s2, zb, zeta1, zeta2, zn, zar;
    double aarg, ang, aphi, ascle, ay, c2i, c2m, c2r, fn, rs1, yy;
    int i, iflag, in, inu, j, k, nai, nd, ndai, nn, nuf, idum;
    double hpi = 1.57079632679489662; /* 0.5 pi */
    double aic = 1.265512123484645396; /* log(2 sqrt(pi)) */
    double complex cip[4] = { 1.0, I, -1.0, -I };
    double complex ci = I;
    //
    // COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAGNITUDE
    // ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
    // EXP(ALIM) = EXP(ELIM)*TOL
    //
    cscl = 1.0 / tol;
    crsc = tol;
    double complex csr[3] = { crsc, 1.0, cscl };
    double complex css[3] = { cscl, 1.0, crsc };
    double bry[3] = { 1.0+3*d1mach[0]/tol, 0.0, 0.0 };
    double complex cy[2] = { 0.0 };
    yy = cimag(z);
    *nz = 0;
    nd = 0;
    *nlast = 0;
    //
    // ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
    //
    zn = -z * ci;
    zb = z;
    cid = -ci;
    inu = (int)fnu;
    ang = hpi * (fnu - inu);
    c2 = CMPLX(cos(ang), sin(ang));
    zar = c2;
    in = inu + n - 1;
    in = in % 4;
    c2 = cip[in];
    if (yy <= 0.0) {
      zn = conj(-zn);
      zb = conj(zb);
      cid = -cid;
      c2 = conj(c2);
    }
    //
    // CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
    //
    fn = fmax(fnu, 1.0);
    unhj(zn, fn, 0, tol, &phi, &arg, &zeta1, &zeta2, &asum, &bsum);
    if (kode != 1) {
        cfn = fnu;
        s1 = -zeta1 + cfn*(cfn/(zb + zeta2));
    } else {
        s1 = -zeta1 + zeta2;
    }
    rs1 = creal(s1);
    if (fabs(rs1) > elim) {
        if (rs1 > 0.) {
            *nz = -1;
            return;
        }
        for (i = 0; i < n; i++) {
            y[i] = 0.0;
        }
        return;
    }
    
L10:
    nn = (nd > 2 ? 2 : nd);
    i = 1;
    for (i = 1; i < (nn+1); i++) {
        fn = fnu + (nd-i);
        unhj(zn, fn, 0, tol, &phi, &arg, &zeta1, &zeta2, &asum, &bsum);
        if (kode != 1) {
            cfn = fnu;
            ay = fabs(yy);
            s1 = -zeta1 + cfn*(cfn/(zb + zeta2)) + ay*I;
        } else {
            s1 = -zeta1 + zeta2;
        }
        //
        // TEST FOR UNDERFLOW AND OVERFLOW
        //
        rs1 = creal(s1);
        if (fabs(rs1) > elim) { goto L50; }
        if (i == 1) { iflag = 2; }
        if (fabs(rs1) >= alim) {
            //
            // REFINE TEST AND SCALE
            //
            aphi = cabs(phi);
            aarg = cabs(arg);
            rs1 += log(aphi) - 0.25*log(aarg) - aic;
            if (fabs(rs1) > elim) { goto L50; }
            if (i == 1) { iflag = 1; }
            if (rs1 >= 0.0){ if (i== 1) { iflag = 3; }}
        }
        //
        // SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
        // EXPONENT EXTREMES
        //
        ai = airy(arg, 0, 2, &nai, &idum);
        dai = airy(arg, 1, 2, &ndai, &idum);
        s2 = phi * (ai*asum + dai*bsum);
        c2r = creal(s1);
        c2i = cimag(s1);
        c2m = exp(c2r)*CMPLX(cos(c2i), sin(c2i));
        s2 *= s1;
        if (iflag == 1) { if (uchk(s1, bry[0], tol)) { goto L50; } }
        if (yy <= 0.0) { s2 = conj(s2); }
        j = nd - i + 1;
        s2 *= c2;
        cy[i-1] = s2;
        y[j-1] = s2*csr[iflag-1];
        c2 *= cid;
    }
    if (nd > 2) {
        rz = 2.0 / z;
        bry[1] = 1.0 / bry[0];
        bry[2] = d1mach[1];
        s1 = cy[0];
        s2 = cy[1];
        c1 = csr[iflag-1];
        ascle = bry[iflag-1];
        k = nd - 2;
        fn = k;
        for (i = 3; i < (nd+1); i++) {
            c2 = s2;
            s2 = s1 + (fnu+fn)*rz*s2;
            s1 = c2;
            c2 = s2*c1;
            y[k-1] = c2;
            k -= 1;
            fn -= 1.0;
            if (iflag < 3) {
                c2r = fabs(creal(c2));
                c2i = fabs(cimag(c2));
                c2m = fmax(c2r, c2i);
                if (c2m > ascle) {
                    iflag += 1;
                    ascle = bry[iflag-1];
                    s1 *= c1;
                    s2 = c2;
                    s1 *= css[iflag-1];
                    s2 *= css[iflag-1];
                    c1 = csr[iflag-1];
                }
            }
        }
    }
    return;

L50:
    if (rs1 <= 0.0) {
        //
        // SET UNDERFLOW AND UPDATE PARAMETERS
        //
        y[nd-1] = 0.0;
        nz += 1;
        nd -= 1;
        if (nd == 0) { return; }
        nuf = uoik(z, fnu, kode, 1, nd, y, tol, elim, alim);
        if (nuf >= 0) {
            nd -= nuf;
            nz += nuf;
            if (nd == 0) { return; }
            fn = fnu + nd - 1;
            if (fn >= fnul) {
                // The following was commented out in the original F77 code
                // C      FN = CIDI
                // C      J = NUF + 1
                // C      K = MOD(J,4) + 1
                // C      S1R = CIPR(K)
                // C      S1I = CIPI(K)
                // C      IF (FN.LT.0.0D0) S1I = -S1I
                // C      STR = C2R*S1R - C2I*S1I
                // C      C2I = C2R*S1I + C2I*S1R                 
                // C      C2R = STR
                in = (inu + nd - 1) % 4;
                c2 = zar*cip[in];
                if (yy <= 0.0) { c2 = conj(c2); }
                goto L10;
            }
            *nlast = nd;
            return;
        }
    }
    *nz = -1;
    return;
}


void unik(
    double complex zr,
    double fnu,
    int ikflg,
    int ipmtr,
    double tol,
    int *init,
    double complex *phi,
    double complex *zeta1,
    double complex *zeta2,
    double complex *total,
    double complex *cwrk
) {
    double complex cfn, crfn, s, sr, t, t2, zn;
    double ac, rfn, test, tstr, tsti;
    int k, l;
    double con[2] = { 3.98942280401432678, 1.25331413731550025 };

    if (init == 0) {
        rfn = 1. / fnu;
        crfn = rfn;

        tstr = creal(zr);
        tsti = cimag(zr);
        test = d1mach[0] * 1e3;
        ac = fnu * test;
        if ((fabs(tstr) <= ac) && (fabs(tsti) <= ac)) {
            ac = 2.*fabs(log(test)) + fnu;
            *zeta1 = ac;
            *zeta2 = fnu;
            *phi = 1.;
        }
        t = zr * crfn;
        s = 1. + t*t;
        sr = sqrt(s);
        cfn = fnu;
        zn = (1. + sr) / t;
        *zeta1 = cfn * log(zn);
        *zeta2 = cfn * sr;
        t = 1. / sr;
        sr = t*crfn;
        cwrk[15] = sqrt(sr);
        *phi = sqrt(sr)*con[ikflg-1];
        if (ipmtr != 0) { return; }
        t2 = 1. / s;
        cwrk[0] = 1.;
        crfn = 1.;
        ac = 1.;
        l = 1;
        k = 2;
        for (int k = 2; k < 16; k++)
        {
            s = 0.;
            for (int j = 1; j < (k+1); j++)
            {
                l += 1;
                s = s*t2 + zunik_c[l-1];
            }
            crfn *= sr;
            cwrk[k-1] = crfn*s;
            ac *= rfn;
            tstr = creal(cwrk[k-1]);
            tsti = cimag(cwrk[k-1]);
            test = fabs(tstr) + fabs(tsti);
            if ((ac < tol) && (test < tol)) {
                break;
            }
        }
        *init = k;
    }
    
    if (ikflg != 2) {
        *total = 0.;
        for (int i = 0; i < (k+1); i++) { *total += cwrk[i]; }
        *phi = cwrk[15] * con[0];
        return;
    }

    s = 0.;
    t = 1.;
    for (int i = 1; i < (k+1); i++) {
        s += t*cwrk[i];
        t = -t;
    }
    *total = s;
    *phi = cwrk[15] * con[1];
    return;
}


int unk1(double complex z,
    double fnu,
    int kode,
    int mr,
    int n,
    double complex *y,
    double tol,
    double elim,
    double alim) {
    double complex cfn, ck, crsc, cs, cscl, csgn, cspn, c1, c2, rz, s1, s2,zr,\
                   phid, zeta1d, zeta2d, sumd;
    double ang, aphi, asc, ascle, c2i, c2m, c2r, fmr, fn, fnf, rs1, sgn, x;
    int i, ib, iflag, ifn, il, inu, iuf, k, kdflg, kflag, kk, m, nw, nz, j,\
        jc, ipard, initd, ic;

    cscl = 1.0 / tol;
    crsc = tol;
    double complex css[3] = {cscl, 1.0, crsc };
    double complex csr[3] = {crsc, 1.0, cscl };
    double complex cwrk[3][16] = {{ 0.0 }};
    double complex phi[2] = { 0.0 };
    double complex sum[2] = { 0.0 };
    double complex zeta1[2] = { 0.0 };
    double complex zeta2[2] = { 0.0 };
    double complex cy[2] = { 0.0 };
    double bry[3] = { 1e3*d1mach[0] / tol, tol / 1e3*d1mach[0], d1mach[1]};
    int init[2] = { 0 };
    double pi = 3.14159265358979324;

    kdflg = 1;
    nz = 0;
    x = creal(z);
    zr = z;
    if (x < 0.0) { zr = -z; }
    j = 2;
    for (i = 1; i < (n+1); i++)
    {
        j = 3 - j; /* j flip flops between 1, 2 */
        jc = j - 1; /* dummy index for 0-indexing */
        fn = fnu + (i - 1);
        unik(zr, fn, 2, 0, tol, &init[j-1], &phi[jc], &zeta1[jc], &zeta2[jc], &sum[jc], cwrk[jc]);
        if (kode != 1) {
            cfn = fn;
            s1 = zeta1[jc] - cfn*(cfn / (zr + zeta2[jc]));
        } else {
            s1 = zeta1[jc] - zeta2[jc];
        }
        //
        // TEST FOR UNDERFLOW AND OVERFLOW
        //
        rs1 = creal(s1);
        if (fabs(rs1) <= elim) {
            if (kdflg == 1) { kflag = 2; }
            if (fabs(rs1) >= alim) {
                //
                // REFINE TEST AND SCALE
                //
                aphi = cabs(phi[jc]);
                rs1 += log(aphi);
                if (fabs(rs1) > elim) { goto L10;}
                if (kdflg == 1) { kflag = 1; }
                if (rs1 >= 0.0) { if (kdflg == 1) { kflag = 3; } }
            }

            //
            // SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
            // EXPONENT EXTREMES
            //
            s2 = phi[jc]*sum[jc];
            c2r = creal(s1);
            c2i = cimag(s1);
            c2m = exp(c2r)*creal(css[kflag-1]);
            s1 = c2m * CMPLX(cos(c2i), sin(c2i));
            s2 *= s1;
            if (kflag == 1) { if (uchk(s2, bry[0], tol)) { goto L10; }}
            cy[kdflg-1] = s2;
            y[i-1] = s2*csr[kflag-1];
            if (kdflg == 2) { goto L30; }
            kdflg = 2;
            continue;
        }
L10:
        if (rs1 > 0.0 ) { return -1; }
        if (x < 0.0) { return -1; }
        kdflg = 1;
        y[i-1] = 0.0;
        nz += 1;
        if (i != 1) {
            if (y[i-2] != 0.0) {
                y[i-2] = 0.0;
                nz += 1;
            }
        }
    }
    i = n;

L30:
    rz = 2.0 / zr;
    ck = fn * rz;
    ib = i + 1;
    if (n >= ib) {
        //
        // TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
        // ON UNDERFLOW
        //
        fn = fnu + (n-1);
        ipard = 1;
        if (mr != 0) { ipard = 0; }
        initd = 0;
        unik(zr, fn, 2, ipard, tol, &initd, &phid, &zeta1d, &zeta2d, &sumd, cwrk[2]);
        if (kode != 1) {
            cfn = fn;
            s1 = zeta1d - cfn*(cfn / (zr + zeta2d));
        } else {
            s1 = zeta1d - zeta2d;
        }
        rs1 = creal(s1);
        if (fabs(rs1) <= elim) {
            if (fabs(rs1) < alim) { goto L50; }
            //
            // REFINE ESTIMATE AND TEST
            //
            aphi = cabs(phid);
            rs1 += log(aphi);
            if (fabs(rs1) < elim) { goto L50; }
        }
        if (rs1 > 0.0) { return -1; }
        //
        // FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
        //
        if (x < 0.0) { return -1; }
        nz = n;
        for (i = 0; i < (n+1); i++) { y[i] = 0.0; }
        return nz;
L50:
        //
        // RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
        //
        s1 = cy[0];
        s2 = cy[1];
        c1 = csr[kflag-1];
        ascle = bry[kflag-1];
        for (i = ib; i < (n+1); i++)
        {
            c2 = s2;
            s2 = ck*s2 + s1;
            s1 = c2;
            ck += rz;
            c2 = s2*c1;
            y[i-1] = c2;
            if (kflag < 3) {
                c2m = fmax(fabs(creal(c2)), fabs(cimag(c2)));
                if (c2m > ascle) {
                    kflag += 1;
                    ascle = bry[kflag-1];
                    s1 *= c1;
                    s2 = c2;
                    s1 *= css[kflag-1];
                    s2 *= css[kflag-1];
                    c1 = csr[kflag-1];
                }
            }
        }
    }
    if (mr == 0) { return nz; }
    //
    // ANALYTIC CONTINUATION FOR RE(Z) < 0.0
    //
    nz = 0;
    fmr = mr;
    sgn = (fmr < 0.0 ? pi : -pi );
    //
    // CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
    //
    csgn = sgn*I;
    inu = (int)fnu;
    fnf = fnu - inu;
    ifn = inu + n - 1;
    ang = fnf * sgn;
    cspn = CMPLX(cos(ang), sin(ang));
    if (ifn % 2 == 1) { cspn = -cspn; }
    asc = bry[0];
    kk = n;
    iuf = 0;
    kdflg = 1;
    ib -= 1;
    ic = ib - 1;
    k = 1;
    for (k = 1; k < (n+1); k++)
    {
        fn = fnu + (kk-1);
        //
        // LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
        // FUNCTION ABOVE
        //
        m = 3;
        if (n > 2) { goto L80; }
L70:
        initd = init[j-1];
        phid = phi[j -1];
        zeta1d = zeta1[j-1];
        zeta2d = zeta2[j-1];
        sumd = sum[j-1];
        m = j;
        j = 3 - j;
        goto L90;
L80:
        if (!((kk == n) && (ib < n))) {
            if ((kk == ib) || (kk == ic)){ goto L70; }
            initd = 0;
        }
L90:
        unik(zr, fn, 1, 0, tol, &initd, &phid, &zeta1d, &zeta2d, &sumd, cwrk[m-1]);
        if (kode != 1) {
            cfn = fn;
            s1 = -zeta1d + cfn * (cfn/(zr + zeta2d));
        } else {
            s1 = -zeta1d + zeta2d;
        }
        //
        // TEST FOR UNDERFLOW AND OVERFLOW
        //
        rs1 = creal(s1);
        if (fabs(rs1) > elim) { goto L110; }
        if (kdflg == 1) { iflag = 2; }
        if (fabs(rs1) >= alim) {
            //
            // REFINE TEST AND SCALE
            //
            aphi = cabs(phid);
            rs1 += log(aphi);
            if (fabs(rs1) > elim) { goto L110; }
            if (kdflg == 1) { iflag = 1; }
            if (rs1 >= 0.0) { if (kdflg == 1) { iflag = 3; } }
        }

        s2 = csgn * phid * sumd;
        c2r = creal(s1);
        c2i = cimag(s1);
        c2m = exp(c2r) * creal(css[iflag-1]);
        s1 = c2m * CMPLX(cos(c2i), sin(c2i));
        s2 = s2 * s1;
        if (iflag == 1) { if (uchk(s2, bry[0], tol)) { s2 = 0.0; } }
L100:
        cy[kdflg -1] = s2;
        c2 = s2;
        s2 *= csr[iflag-1];
        //
        // ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
        //
        s1 = y[kk-1];
        if (kode != 1) {
            nw = s1s2(zr, &s1, &s2, asc, alim, &iuf);
            nz += nw;
        }
        y[kk-1] = s1*cspn + s2;
        kk -= 1;
        cspn = -cspn;
        if (c2 == 0.0) {
            kdflg = 1;
            continue;
        }
        if (kdflg == 2) { goto L130; }
        kdflg = 2;
        continue;
L110:
        if (rs1 > 0.0) { return -1; }
        s2 = 0.0;
        goto L100;
    }
    k = n;

L130:
    il = n-k;
    if (il == 0) { return nz; };
    s1 = cy[0];
    s2 = cy[1];
    cs = csr[iflag-1];
    ascle = bry[iflag-1];
    fn = inu + il;
    for (i = 1; i < (il+1); i++)
    {
        c2 = s2;
        s2 = s1 + (fn + fnf) * rz * s2;
        s1 = c2;
        fn -= 1.0;
        c2 = s2 * cs;
        ck = c2;
        c1 = y[kk-1];
        if (kode != 1) {
            nw = s1s2(zr, &c1, &c2, asc, alim, &iuf);
            nz = nz + nw;
        }
        y[kk-1] = c1 * cspn + c2;
        kk -= 1;
        cspn = -cspn;
        if (iflag < 3) {
            c2m = fmax(fabs(creal(c2)), fabs(cimag(c2)));
            if (c2m > ascle) {
                iflag += 1;
                ascle = bry[iflag-1];
                s1 *= cs;
                s2 = ck;
                s1 *= css[iflag-1];
                s2 *= css[iflag-1];
                cs = csr[iflag-1];
            }
        }
    }
    return nz;
};


int unk2(double complex z,
    double fnu,
    int kode,
    int mr,
    int n,
    double complex *y,
    double tol,
    double elim,
    double alim
) {
    double complex ai, cfn, ck, cs, csgn, cspn, c1, c2, dai, rz, s1, s2,\
                 zb, zn, zr, phid, argd, zeta1d, zeta2d, asumd, bsumd;
    double aarg, ang, aphi, asc, ascle, car, cpn, c2i, c2m, c2r, crsc, cscl,\
           fmr, fn, fnf, rs1, sar, sgn, spn, x, yy;
    int i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, kk, nai, ndai,\
        nw, nz, idum, j, ipard, ic;

    double complex cr1 = CMPLX(1.0, 1.73205080756887729);
    double complex cr2 = CMPLX(-0.5, -8.66025403784438647e-1);
    double hpi = 1.57079632679489662;  /* 0.5 pi */
    double pi = 3.14159265358979324;
    double aic = 1.26551212348464539;  /* log(2 sqrt(pi)) */
    double complex cip[4] = {1.0, I, -1.0, -I};
    cscl = 1.0 / tol;
    crsc = tol;
    double complex css[3] = {cscl, 1.0, crsc };
    double complex csr[3] = {crsc, 1.0, cscl };
    double complex phi[2] = { 0.0 };
    double complex arg[2] = { 0.0 };
    double complex zeta1[2] = { 0.0 };
    double complex zeta2[2] = { 0.0 };
    double complex asum[2] = { 0.0 };
    double complex bsum[2] = { 0.0 };
    double complex cy[2] = { 0.0 };
    double bry[3] = { 1e3*d1mach[0] / tol, tol / 1e3*d1mach[0], d1mach[1]};

    kdflg = 1;
    nz = 0;
    //
    // EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
    // THE UNDERFLOW LIMIT
    //
    x = creal(z);
    zr = z;
    if (x < 0.0) { zr = -z; }
    yy = cimag(zr);
    zn = -zr*I;
    zb = zr;
    inu = (int)fnu;
    fnf = fnu - inu;
    ang = -hpi * fnf;
    car = cos(ang);
    sar = sin(ang);
    cpn = -hpi * car;
    spn = -hpi * sar;
    c2 = CMPLX(-spn, cpn);
    kk = (inu % 4) + 1;
    cs = cr1 * c2 * cip[kk - 1];
    if (yy <= 0.0) {
        zn = conj(-zn);
        zb = conj(zb);
    }
    //
    // K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
    // QUADRANT.  FOURTH QUADRANT VALUES (YY <= 0.0_dp) ARE COMPUTED BY
    // CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
    //
    j = 2;
    for (i = 1; i < (n+1); i++)
    {
        j = 3 - j;
        fn = fnu + (i-1);
        unhj(zn, fn, 0, tol, &phi[j-1], &arg[j-1], &zeta1[j-1], &zeta2[j-1], &asum[j-1],&bsum[j-1]);
        if (kode != 1) {
            cfn = fn;
            s1 = zeta1[j-1] - cfn*(cfn/(zb + zeta2[j-1]));
        } else {
            s1 = zeta1[j-1] - zeta2[j-1];
        }
        //
        // TEST FOR UNDERFLOW AND OVERFLOW
        //
        rs1 = creal(s1);
        if (fabs(rs1) <= elim) {
            if (kdflg == 1) { kflag = 2; }
            if (fabs(rs1) >= alim) {
                //
                // REFINE TEST AND SCALE
                aphi = cabs(phi[j-1]);
                aarg = cabs(arg[j-1]);
                rs1 += log(aphi) - 0.25 * log(aarg) - aic;
                if (fabs(rs1) > elim) { goto L10; }
                if (kdflg == 1) { kflag = 1; }
                if (rs1 >= 0.0) { if (kdflg == 1) { kflag = 3; } }
            }
            //
            // SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
            // EXPONENT EXTREMES
            //
            c2 = arg[j-1] * cr2;
            ai = airy(c2, 0, 2, &nai, &idum);
            dai = airy(c2, 1, 2, &ndai, &idum);
            s2 = cs * phi[j-1] * (ai*asum[j-1] + cr2*dai*bsum[j-1]);
            c2r = creal(s1);
            c2i = cimag(s1);
            c2m = exp(c2r) * creal(css[kflag-1]);
            s1 = c2m * CMPLX(cos(c2i), sin(c2i));
            s2 *= s1;
            if (kflag == 1) { if (uchk(s2, bry[0], tol)) { goto L10; } };
            if (yy <= 0.0) { s2 = conj(s2); }
            cy[kdflg-1] = s2;
            y[i-1] = s2 * csr[kflag-1];
            cs *= -I;
            if (kdflg == 2) { goto L30; }
            kdflg = 2;
            continue;
        }
L10:
        if (rs1 > 0.0) { return -1; }
        //
        // FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
        //
        if (x < 0.0) { return -1; }
        kdflg = 1;
        y[i-1] = 0.0;
        cs *= -I;
        nz += 1;
        if (i != 1) {
            if (y[i-2] != 0.0) {
                y[i-2] = 0.0;
                nz += 1;
            }
        }
    }
    i = n;
L30:
    rz = 2.0 / zr;
    ck = fn * rz;
    ib = i + 1;
    if (n >= ib) {
        fn = fnu + (n - 1);
        ipard = 1;
        if (mr != 0) { ipard = 0; }
        unhj(zn, fn, ipard, tol, &phid, &argd, &zeta1d, &zeta2d, &asumd, &bsumd);
        if (kode != 1) {
            cfn = fn;
            s1 = zeta1d - cfn * (cfn / (zb + zeta2d));
        } else {
            s1 = zeta1d - zeta2d;
        }
        rs1 = creal(s1);
        if (fabs(rs1) <= elim) {
            if (fabs(rs1) < alim) { goto L50; }
            //
            // REFINE ESTIMATE AND TEST
            //
            aphi = cabs(phid);
            aarg = cabs(argd);
            rs1 += log(aphi) - 0.25 * log(aarg) - aic;
            if (fabs(rs1) < elim) { goto L50; }
        }
        if (rs1 > 0.0) { return -1; }
        //
        //  FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
        //
        if (x < 0.0) { return -1; }
        nz = n;
        for (i = 0; i < n; i++) { y[i] = 0.0; }
        return nz;
L50:
        //
        // SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
        //
        s1 = cy[0];
        s2 = cy[1];
        c1 = csr[kflag-1];
        ascle = bry[kflag-1];
        for (i = ib; i < (n+1); i++) 
        {
            c2 = s2;
            s2 = ck * s2 + s1;
            s1 = c2;
            ck += rz;
            c2 = s2 * c1;
            y[i-1] = c2;
            if (kflag < 3) {
                c2m = fmax(fabs(creal(c2)), fabs(cimag(c2)));
                if (c2m > ascle) {
                    kflag += 1;
                    ascle = bry[kflag-1];
                    s1 *= c1;
                    s2 = c2;
                    s1 *= css[kflag-1];
                    s2 *= css[kflag-1];
                    c1 = csr[kflag-1];
                }
            }
        }
    }
    if (mr == 0) { return nz; }
    //
    // ANALYTIC CONTINUATION FOR RE(Z) < 0.0_dp
    //
    nz = 0;
    fmr = mr;
    sgn = ( fmr < 0.0 ? -pi : pi);
    //
    // CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
    //
    csgn = sgn*I;
    if (yy <= 0.0) { csgn = -csgn; }
    ifn = inu + n - 1;
    ang = fnf*sgn;
    cspn = CMPLX(cos(ang), sin(ang));
    if (ifn % 2 == 1) { cspn = -cspn; }
    //
    // CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION.  I(FNU,Z) IS
    // COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
    // QUADRANT.  FOURTH QUADRANT VALUES (YY <= 0.0_dp) ARE COMPUTED BY
    // CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
    //
    cs = CMPLX(car, -sar) * csgn;
    in = (ifn % 4) + 1;
    c2 = cip[in-1];
    cs *= conj(c2);
    asc = bry[0];
    kk = n;
    kdflg = 1;
    ib -= 1;
    ic = ib - 1;
    iuf = 0;
    for (k = 1; k <= (n+1); k++) {
        fn = fnu + (kk-1);
        if (n > 2) { goto L80; }
L70:
        phid = phi[j-1];
        argd = arg[j-1];
        zeta1d = zeta1[j-1];
        zeta2d = zeta2[j-1];
        asumd = asum[j-1];
        bsumd = bsum[j-1];
        j = 3 - j;
        goto L90;
L80:
        if (!((kk == n) && (ib < n))) {
            if ((kk == ib) || (kk == ic)) { goto L70; }
            unhj(zn, fn, 0, tol, &phid, &argd, &zeta1d, &zeta2d, &asumd, &bsumd);
        }
L90:
        if (kode != 1) {
            cfn = fn;
            s1 = -zeta1d + cfn * (cfn/(zb + zeta2d));
        } else {
            s1 = -zeta1d + zeta2d;
        }
        //
        // TEST FOR UNDERFLOW AND OVERFLOW
        //
        rs1 = creal(s1);
        if (fabs(rs1) > elim) { goto L110; }
        if (kdflg == 1) { iflag = 2; }
        if (fabs(rs1) >= alim) {
            aphi = cabs(phid);
            aarg = cabs(argd);
            rs1 += log(aphi) - 0.25f * log(aarg) - aic;
            if (fabs(rs1) > elim) { goto L110; }
            if (kdflg == 1) { iflag = 1; }
            if (rs1 >= 0.0) { if (kdflg == 1) {iflag = 3;} }
        }

        ai = airy(argd, 0, 2, &nai, &idum);
        dai = airy(argd, 1, 2, &ndai, &idum);
        s2 = cs * phid * (ai*asumd + dai*bsumd);
        c2r = creal(s1);
        c2i = cimag(s1);
        c2m = exp(c2r) * creal(css[iflag-1]);
        s1 = c2m * CMPLX(cos(c2i), sin(c2i));
        s2 *= s1;
        if (iflag == 1) { if (uchk(s2, bry[0], tol)) { s2 = 0.0; } }

L100:
        if (yy <= 0.0) { s2 = conj(s2); }
        cy[kdflg-1] = s2;
        c2 = s2;
        s2 *= csr[iflag-1];
        //
        // ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
        //
        s1 = y[kk-1];
        if (kode != 1) {
            nw = s1s2(zr, &s1, &s2, asc, alim, &iuf);
            nz += nw;
        }
        y[kk-1] = s1 * cspn + s2;
        kk -= 1;
        cspn = -cspn;
        cs *= -I;
        if (c2 == 0.0) {
            kdflg = 1;
            continue;
        }
        if (kdflg == 2) { goto L130; }
        kdflg = 2;
        continue;

L110:
        if (rs1 > 0.0) { return -1; }
        s2 = 0.0;
        goto L100;

    }
    k = n;

L130:
    il = n - k;
    if (il == 0) { return nz; }
    //
    // RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
    // K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
    // INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
    //
    s1 = cy[0];
    s2 = cy[1];
    cs = csr[iflag-1];
    ascle = bry[iflag-1];
    fn = inu + il;
    for (i = 1; i < (il+1); i++) {
        c2 = s2;
        s2 = s1 + (fn + fnf) * rz * s2;
        s1 = c2;
        fn -= 1.0;
        c2 = s2 * cs;
        ck = c2;
        c1 = y[kk-1];
        if (kode != 1) {
            nw = s1s2(zr, &c1, &c2, asc, alim, &iuf);
            nz = nz + nw;
        }
        y[kk-1] = c1 * cspn + c2;
        kk -= 1;
        cspn = -cspn;
        if (iflag < 3) {
            c2m = fmax(fabs(creal(ck)), fabs(cimag(ck)));
            if (c2m > ascle) {
                iflag += 1;
                ascle = bry[iflag-1];
                s1 *= cs;
                s2 = ck;
                s1 *= css[iflag-1];
                s2 *= css[iflag-1];
                cs = csr[iflag-1];
            }
        }
    }
    return nz;
};


int uoik(
    double complex z,
    double fnu,
    int kode,
    int ikflg,
    int n,
    double complex *y,
    double tol,
    double elim,
    double alim
) {
    double complex arg, asum, bsum, cz, phi, sum, zb, zeta1;
    double complex zeta2, zn, zr;
    double aarg, aphi, ascle, ax, ay, fnn, gnn, gnu, rcz, x, yy;
    int iform, init, nn;
    double aic = 1.265512123484645396;
    double complex cwrk[16] = { 0. };

    int nuf = 0;
    nn = n;
    x = creal(z);
    zr = z;
    if (x < 0.) { zr = -z; }
    zb = zr;
    yy = cimag(zr);
    ax = fabs(x) * sqrt(3.);
    ay = fabs(yy);
    iform = 1;
    if (ay > ax) { iform = 2; }
    gnu = fmax(fnu, 1.);
    if (ikflg != 1) {
        fnn = nn;
        gnn = fnu + fnn -1;
        gnu = fmax(gnn, fnn);
    }

    if (iform != 2) {
        init = 0;
        unik(zr, gnu, ikflg, 1, tol, &init, &phi, &zeta1, &zeta2, &sum, &cwrk[0]);
        cz = -zeta1 + zeta2;
    } else {
        zn = -zr * I;
        if (yy <= 0.) {
            zn = conj(zn);
        }
        unhj(zn, gnu, 1, tol, &phi, &arg, &zeta1, &zeta2, &asum, &bsum);
        cz = zeta2 - zeta1;
        aarg = cabs(arg);
    }
    if (kode == 2) { cz -= zb; }
    if (ikflg == 2) { cz = -cz; }
    aphi = cabs(phi);
    rcz = creal(cz);

    /*  OVERFLOW TEST  */
    if (rcz > elim) { return -1; }
    if (rcz >= alim) {
        rcz += log(aphi);
        if (iform == 2) { rcz -= 0.25*log(aarg) + aic; }
        if (rcz > elim) { return -1; }
    } else {
        /*  UNDERFLOW TEST  */
        if (rcz >= -elim) {
            if (rcz > -alim) {
                /* pass */
            } else {
                rcz += log(aphi);
                if (iform == 2) { rcz -= 0.25*log(aarg) + aic; }
                if (rcz > -elim) {
                    /* goto 30 */
                    ascle = 1e3*d1mach[0] / tol;
                    cz += clog(phi);
                    if (iform != 1) { cz -= 0.25*log(arg) + aic;}
                    ax = exp(rcz) / tol;
                    ay = cimag(cz);
                    cz = ax*(cos(ay)+sin(ay)*I);
                    if (uchk(cz, ascle, tol)) {
                        for (int i = 0; i < nn; i++){ y[i] = 0.; }
                        return nn;
                    }
                } else {
                    for (int i = 0; i < nn; i++){ y[i] = 0.; }
                    return nn;
                }
            }
        } else {
            for (int i = 0; i < nn; i++){ y[i] = 0.; }
            return nn;
        }
    }
    if ((ikflg == 2) || (n == 1)) { return nuf; }
    /* 140 */
    while (1) {
        gnu = fnu + (nn -1);
        if (iform != 2) {
            init = 0;
            unik(zr, gnu, ikflg, 1, tol, &init, &phi, &zeta1, &zeta2, &sum, &cwrk[0]);
            cz = zeta2 - zeta1;
        } else {
            unhj(zn, gnu, 1, tol, &phi, &arg, &zeta1, &zeta2, &asum, &bsum);
            cz = zeta2 - zeta1;
            aarg = cabs(phi);
        }
        if (kode == 2) { cz -= zb; }

        /* 170 */
        aphi = cabs(phi);
        rcz = creal(cz);

        if (rcz >= -elim) {
            if (rcz > -alim) { return nuf; }
            rcz += log(aphi);
            if (iform == 2) { rcz -= 0.25*log(aarg) + aic; }
            if (rcz > -elim) {
                ascle = 1e3 * d1mach[0] / tol;
                cz = clog(phi);
                if (iform != 1) { cz -= 0.25*clog(arg) + aic; }
                ax = exp(rcz)/tol;
                ay = cimag(cz);
                cz = ax*(cos(ay)+sin(ay*I));
                if (!(uchk(cz, ascle, tol))) { return nuf; }
            }
        }

        y[nn-1] = 0.;
        nn -= 1;
        nuf += 1;
        if (nn == 0) { return nuf; }
    }
    return -1;
}


int wrsk(
    double complex zr,
    double fnu,
    int kode,
    int n,
    double complex *y,
    double complex *cw,
    double tol,
    double elim,
    double alim
) {
   double complex cinu, cscl, ct, c1, c2, rct, st;
   double act, acw, ascle, yy;
   int i, nw, nz;
 
    //
    // I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
    // Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
    // WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
    //
    nz = 0;
    nw = bknu(zr, fnu, kode, 2, cw, tol, elim, alim);
    if (nw != 0) {
        nz = -1;
        if (nw == -2) {
            nz = -2;
        }
        return nz;
    }
    
    rati(zr, fnu, 2, y, tol);
    //
    // RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    // R(FNU+J-1,Z)=Y(J),  J=1,...,N
    //
    cinu = 1.0;
    if (kode != 1) {
        yy = cimag(zr);
        cinu = CMPLX(cos(yy), sin(yy));
    }
    //
    // ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH THE
    // UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE SCALED TO
    // PREVENT OVER OR UNDERFLOW.  CUOIK HAS DETERMINED THAT THE RESULT
    // IS ON SCALE.
    //
    acw = cabs(cw[1]);
    ascle = 1.0 + 3*d1mach[0]/tol;
    cscl = 1.0;
    
    if (acw <= ascle) {
        cscl = 1.0 / tol;
    } else {
        ascle = 1.0 / ascle;
        if (acw >= ascle) {
            cscl = tol;
        }
    }
    
    c1 = cw[0]*cscl;
    c2 = cw[1]*cscl;
    st = y[0];
    //
    // CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0_dp/ABS(CT) PREVENTS
    // UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
    //
    ct = zr * (c2 + st*c1);
    act = cabs(ct);
    rct = 1.0 / act;
    ct = conj(ct)*rct;
    cinu *= ct*rct;
    y[0] = cinu*cscl;
    if (n == 1) { return nz; }
    for (i = 2; i < (n+1); i++) {
        cinu *= st;
        st = y[i-1];
        y[i-1] = cinu*cscl;
    }
    return nz;
}
