# cython: boundscheck=False
# cython: initializedcheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: cpow=True

import numpy as np
cimport numpy as cnp
cnp.import_array()
from libc.math cimport sin, tan, log, exp, sqrt, floor
from numpy.math cimport INFINITY, PI

cdef double[3] spmpar = [np.finfo(np.float64).eps,
                         np.finfo(np.float64).tiny,
                         np.finfo(np.float64).max]


# %%----------------------------------------- algdiv
cdef inline double algdiv(double a, double b) noexcept nogil:
    cdef double c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2
    cdef double[6] carr = [0.833333333333333e-01, -0.277777777760991e-02,
                           0.793650666825390e-03, -0.595202931351870e-03,
                           0.837308034031215e-03, -0.165322962780713e-02]

    if a > b:
        h = b / a
        c = 1./(1. + h)
        x = h/(1. + h)
        d = a + (b - 0.5)
    else:
        h = a / b
        c = h/(1. + h)
        x = 1./(1. + h)
        d = b + (a - 0.5)

    x2 = x*x
    s3 = 1. + (x + x2)
    s5 = 1. + (x + x2*s3)
    s7 = 1. + (x + x2*s5)
    s9 = 1. + (x + x2*s7)
    s11 = 1. + (x + x2*s9)

    t = (1. / b)**2
    w = (((((carr[5]*s11
            )*t + carr[4]*s9
           )*t + carr[3]*s7
          )*t + carr[2]*s5
         )*t + carr[1]*s3
        )*t + carr[0]
    w *= c / b

    u = d * alnrel(a / b)
    v = a * (log(b) - 1.)
    return (w - v) - u if (u>v) else (w - u) - v

# %%----------------------------------------- alngam
cdef inline double alngam(double x) noexcept nogil:
    cdef double prod, xx, result, offset
    cdef int i, n
    cdef double[9] scoefn = [0.62003838007127258804e2, 0.36036772530024836321e2,
                             0.20782472531792126786e2, 0.6338067999387272343e1,
                             0.215994312846059073e1, 0.3980671310203570498e0,
                             0.1093115956710439502e0, 0.92381945590275995e-2,
                             0.29737866448101651e-2]
    cdef double[4] scoefd = [0.62003838007126989331e2, 0.9822521104713994894e1,
                             -0.8906016659497461257e1, 0.1000000000000000000e1]
    cdef double[5] coef = [0.83333333333333023564e-1, -0.27777777768818808e-2,
                           0.79365006754279e-3, -0.594997310889e-3, 0.8065880899e-3]

    if x <= 6.:
        prod = 1.
        xx = x

        if x > 3.:
            while xx > 3.:
                xx -= 1.
                prod *= xx

        if x < 2.:
            while xx < 2.:
                prod /= xx
                xx += 1.

        result = devlpl(scoefn, 9, xx - 2.) / devlpl(scoefd, 4, xx - 2.)
        return log(result * prod)

    offset = 0.5*log(2.*PI)
    if x <= 12.:
        n = <int>(12. - x)
        if n > 0:
            prod = 1.
            for i in range(n):
                prod *= x + i
            offset -= log(prod)
            xx = x + n
        else:
            xx = x
    else:
        xx = x

    result = devlpl(coef, 5, (1./xx)**2) / xx
    result += offset + (xx - 0.5)*log(xx) - xx
    return result

# %%----------------------------------------- alnrel
cdef inline double alnrel(double a) noexcept nogil:
    cdef double[3] p = [-0.129418923021993e+01, 0.405303492862024e+00,
                        -0.178874546012214e-01]
    cdef double[3] q = [-0.162752256355323e+01, 0.747811014037616e+00,
                        -0.845104217945565e-01]
    cdef double t, t2, w

    if abs(a) > 0.375:
        return log(1. + a)
    else:
        t = a / (a + 2.)
        t2 = t*t
        w = ((p[2]*t2 + p[1])*t2 + p[0])*t2 + 1.
        w /= ((q[2]*t2 + q[1])*t2 + q[0])*t2 + 1.
        return 2*t*w

# %%----------------------------------------- apser
cdef inline double apser(double a, double b, double x, double eps) noexcept nogil:
    cdef double aj, bx, c, j, s, t, tol
    cdef double g = 0.577215664901532860606512090082

    bx = b*x
    t = x - bx
    if b*eps > 0.02:
        c = log(bx) + g + t
    else:
        c = log(x) + psi(b) + g + t

    tol = 5.*eps*abs(c)
    j = 1.
    s = 0.

    while True:
        j += 1
        t *= (x - bx/j)
        aj = t / j
        s += aj
        if abs(aj) <= tol:
            break

    return -a * (c + s)

# %%-------------------------------------- basym
cdef inline double basym(double a, double b, double lmbda, double eps) noexcept nogil:
    cdef double a0[21]
    cdef double b0[21]
    cdef double c[21]
    cdef double d[21]
    cdef double bsum, dsum, f, h, h2, hn, j0, j1, r, r0, r1, s, ssum
    cdef double t, t0, t1, u, w, w0, z, z0, z2, zn, znm1
    cdef double e0 = 2. / sqrt(PI)
    cdef double e1 = 2.**(-3./2.)
    cdef int i, imj, j, m, mmj, n, num

    num = 20

    if a < b:
        h = a / b
        r0 = 1./(1.+h)
        r1 = (b - a) / b
        w0 = 1. / sqrt(a * (1. + h))
    else:
        h = b / a
        r0 = 1./(1.+h)
        r1 = (b - a) / a
        w0 = 1. / sqrt(b * (1. + h))

    f = a*rlog1(-lmbda/a) + b*rlog1(lmbda/b)
    t = exp(-f)
    if t == 0.:
        return 0.
    z0 = sqrt(f)
    z = 0.5*(z0/e1)
    z2 = f + f

    a0[0] = (2./3.)*r1
    c[0] = -0.5*a0[0]
    d[0] = -c[0]
    j0 = (0.5/e0)*erfc1(1, z0)
    j1 = e1
    ssum = j0 + d[0]*w0*j1

    s, h2, hn, w, znm1, zn = 1., h*h, 1., w0, z, z2

    for n in range(2, num+1, 2):
        hn *= h2
        a0[n-1] = 2.*r0*(1.+ h*hn)/(n + 2.)
        s += hn
        a0[n] = 2.*r1*s/(n + 3.)

        for i in range(n, n+2):
            r = -0.5*(i + 1.)
            b0[0] = r*a0[0]

            for m in range(2, i+1):
                bsum = 0.
                for j in range(1, m):
                    mmj = m - j
                    bsum += (j*r - mmj)*a0[j-1]*b0[mmj-1]
                b0[m-1] = r*a0[m-1] + bsum/m

            c[i-1] = b0[i-1] / (i + 1.)
            dsum = 0.

            for j in range(1, i):
                imj = i - j
                dsum += d[imj-1]*c[j-1]
            d[i-1] = - (dsum+c[i-1])

        j0 = e1*znm1 + (n-1.)*j0
        j1 = e1*zn + n*j1
        znm1 *= z2
        zn *= z2
        w *= w0
        t0 = d[n-1]*w*j0
        w *= w0
        t1 = d[n]*w*j1
        ssum += t0 + t1
        if (abs(t0) + abs(t1)) <= eps*ssum:
            break

    u = exp(-bcorr(a, b))

    return e0*t*u*ssum

# %%-------------------------------------- bcorr
cdef inline double bcorr(double a0, double b0) noexcept nogil:
    cdef double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2
    cdef double[6] carr = [0.833333333333333e-01, -0.277777777760991e-02,
                           0.793650666825390e-03, -0.595202931351870e-03,
                           0.837308034031215e-03, -0.165322962780713e-02]

    a, b = min(a0, b0), max(a0, b0)
    h = a / b
    c = h/(1. + h)
    x = 1./(1. + h)
    x2 = x*x
    s3 = 1. + (x + x2)
    s5 = 1. + (x + x2*s3)
    s7 = 1. + (x + x2*s5)
    s9 = 1. + (x + x2*s7)
    s11 = 1. + (x + x2*s9)

    t = (1. / b)**2
    w = (((((carr[5]*s11
            )*t + carr[4]*s9
           )*t + carr[3]*s7
          )*t + carr[2]*s5
         )*t + carr[1]*s3
        )*t + carr[0]
    w *= c / b

    t = (1. / a)**2
    return ((((((carr[5])*t + carr[4]
               )*t + carr[3]
              )*t + carr[2]
             )*t + carr[1]
            )*t + carr[0]
           )/a + w

# %% betaln ----------------------------------------- betaln
cdef inline double betaln(double a0, double b0) noexcept nogil:
    cdef double a, b, c, h, u, v, w, z
    cdef double e = .918938533204673
    cdef int i, n

    a, b = min(a0, b0), max(a0, b0)
    if a >= 8.:
        w = bcorr(a, b)
        h = a / b
        c = h/(1. + h)
        u = -(a - 0.5)*log(c)
        v = b*alnrel(h)
        if u > v:
            return (((-0.5*log(b)+e)+w)-v) - u
        else:
            return (((-0.5*log(b)+e)+w)-u) - v

    if a < 1:
        if b > 8:
            return gamln(a) + algdiv(a,b)
        else:
            return gamln(a) + (gamln(b) - gamln(a+b))


    if a <= 2:
        if b <= 2:
            return gamln(a) + gamln(b) - gsumln(a, b)

        if b >= 8:
            return gamln(a) + algdiv(a, b)

        w = 0.


    if a > 2:
        if b <= 1000:
            n = <int>(a - 1.)
            w = 1.
            for i in range(n):
                a -= 1.
                h = a / b
                w *= h/(1.+h)
            w = log(w)
            if b >= 8.:
                return w + gamln(a) + algdiv(a, b)
        else:
            n = <int>(a - 1.)
            w = 1.
            for i in range(n):
                a -= 1.
                w *= a/(1. + (a/b))
            return (log(w) - n*log(b)) + (gamln(a) + algdiv(a, b))

    n = <int>(b - 1.)
    z = 1.
    for i in range(n):
        b -= 1.
        z *= b / (a + b)
    return w + log(z) + (gamln(a) + gamln(b) - gsumln(a, b))

# %%----------------------------------------- bfrac
cdef inline double bfrac(double a, double b, double x, double y,
                         double lmbda, double eps) noexcept nogil:
    cdef double alpha, beta, e, r0, t, w, result
    cdef:
        double c = 1. + lmbda
        double c0 = b / a
        double c1 = 1. + (1. / a)
        double yp1 = y + 1.
        double n = 0.
        double p = 1.
        double s = a + 1.
        double an = 0.
        double bn = 1.
        double anp1 = 1.
        double bnp1 = c / c1
        double r = c1 / c

    result = brcomp(a, b, x, y)

    if result == 0.:
        return 0

    while True:
        n += 1.
        t = n / a
        w = n * (b - n)*x
        e = a / s
        alpha = (p*(p + c0)*e*e) * (w*x)
        e = (1. + t) / (c1 + t + t)
        beta = n + (w / s) + e*(c + n*yp1)
        p = 1. + t
        s += 2.
        t = alpha*an + beta*anp1
        an = anp1
        anp1 = t
        t = alpha*bn + beta*bnp1
        bn = bnp1
        bnp1 = t
        r0 = r
        r = anp1 / bnp1
        if (abs(r - r0) <= eps*r):
            break

        an /= bnp1
        bn /= bnp1
        anp1 = r
        bnp1 = 1.

    return result*r

# %%----------------------------------------- bgrat
cdef inline (double, int) bgrat(double a, double b, double x , double y, double w,
                                double eps) noexcept nogil:
    cdef double bp2n, cn, coef, dj, j, l, n2, q, r, s, ssum, t, t2, u, v
    cdef double c[30]
    cdef double d[30]
    cdef double bm1 = (b - 0.5) - 0.5
    cdef double nu = a + bm1*0.5
    cdef double lnx = log(x) if (y > 0.375) else alnrel(-y)
    cdef double z = -nu*lnx
    cdef int i, n

    if (b*z) == 0.:
        return (w, 1)

    r = b * (1. + gam1(b)) * exp(b*log(z))
    r *= exp(a*lnx) * exp(0.5*bm1*lnx)
    u = algdiv(b, a) + b*log(nu)
    u = r*exp(-u)
    if u == 0.:
        return (w, 1)
    _, q = grat1(b, z, r, eps)
    v = 0.25 * (1 / nu)**2
    t2 = 0.25*lnx*lnx
    l = w / u
    j = q / r
    ssum = j
    t = 1.
    cn = 1.
    n2 = 0.

    for n in range(30):
        bp2n = b + n2
        j = (bp2n*(bp2n + 1.)*j + (z + bp2n + 1.)*t)*v
        n2 += 2.
        t *= t2
        cn /= (n2*(n2 + 1.))
        c[n] = cn
        s = 0.
        if n > 0:
            coef = b - (n+1)
            for i in range(n):
                s += coef*c[i]*d[n-i-1]
                coef += b

        d[n] = bm1*cn + s/(n+1)
        dj = d[n]*j
        ssum += dj
        if ssum <= 0.:
            return (w, 1)
        if abs(dj) <= eps*(ssum+l):
            break

    return (w + u*ssum, 0)

# %%----------------------------------------- bpser
cdef inline double bpser(double a, double b, double x, double eps) noexcept nogil:
    cdef double a0, apb, b0, c, n, ssum, t, tol, u, w, z, result
    cdef int i, m

    if x == 0.:
        return 0.

    a0 = min(a, b)

    if a0 < 1.:
        b0 = max(a, b)
        if b0 <= 1.:
            result = x**a
            if result == 0.:
                return 0.

            apb = a + b
            if apb > 1.:
                u = a + b - 1.
                z = (1. + gam1(u)) / apb
            else:
                z = 1. + gam1(apb)

            c = (1. + gam1(a)) * (1. + gam1(b))/z
            result *= c * (b / apb)

        if b0 < 8:
            u = gamln1(a0)
            m = <int>(b0 - 1.)
            if m > 0:
                c = 1.
                for i in range(m):
                    b0 -= 1.
                    c *= (b0 / (a0 + b0))
                u += log(c)

            z = a*log(x) - u
            b0 -= 1.
            apb = a0 + b0
            if apb > 1.:
                u = a0 + b0 - 1.
                t = (1. + gam1(u)) / apb
            else:
                t = 1. + gam1(apb)

            result = exp(z) * (a0 / a) * (1. + gam1(b0)) / t

        if b0 >= 8:
            u = gamln1(a0) + algdiv(a0, b0)
            z = a*log(x) - u
            result = (a0 / a) * exp(z)
    else:
        z = a*log(x) - betaln(a, b)
        result = exp(z) / a

    if (result == 0.) or (a <= 0.1*eps):
                return result

    ssum = 0.
    n = 0.
    c = 1.
    tol = eps / a
    while True:
        n += 1.
        c *= (0.5 + (0.5 - b/n))*x
        w = c / (a + n)
        ssum += w
        if abs(w) <= tol:
            break
    return result * (1. + a*ssum)

# %%----------------------------------------- bratio
cdef inline (double, double, int) bratio(double a, double b,
                                         double x, double y) noexcept nogil:
    cdef double a0, b0, lmbda, x0, y0
    cdef int ierr1, ind, n
    cdef double w = 0.
    cdef double w1 = 0.
    cdef double eps = max(spmpar[0], 1e-15)

    if (a < 0.) or (b < 0.):
        return (w, w1, 1)
    elif (a == 0.) and (b == 0.):
        return (w, w1, 2)
    elif (x < 0.) or (x > 1.):
        return (w, w1, 3)
    elif (y < 0.) or (y > 1.):
        return (w, w1, 4)
    elif abs(((x+y)-0.5) - 0.5) > 3.*eps:
        return (w, w1, 5)
    elif (x == 0.):
        return (w, w1, 6) if (a == 0.) else (0., 1., 0)
    elif (y == 0.):
        return (w, w1, 7) if (b == 0.) else (1., 0., 0)
    elif (a == 0.):
        return (1., 0., 0)
    elif (b == 0.):
        return (0., 1., 0)
    elif max(a, b) < 1e-3*eps:
        return (b/(a+b), a/(a+b), 0)

    ind = 0
    a0, b0 = a, b
    x0, y0 = x, y

    if min(a0, b0) <= 1.:
        if x > 0.5:
            ind = 1
            a0, b0 = b, a
            x0, y0 = y, x

        if b0 < min(eps, eps*a0):
            w = fpser(a0, b0, x0, eps)
            w1 = 0.5 + (0.5 - w)
            return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        if (a0 < min(eps, eps*b0)) and (b0*x0 <= 1.):
            w1 = apser(a0, b0, x0, eps)
            w = 0.5 + (0.5 - w1)
            return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        if max(a0, b0) <= 1.:
            if ((x0**a0) <= 0.9) or (a0 >= min(0.2, b0)):
                w = bpser(a0, b0, x0, eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            if (x0 >= 0.3):
                w1 = bpser(b0, a0, y0, eps)
                w = 0.5 + (0.5 - w1)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            # 140
            n = 20
            w1 = bup(b0, a0, y0, x0, n, eps)
            b0 += n
            w1, ierr1 = bgrat(b0, a0, y0, x0, w1, 15.*eps)
            w = 0.5 + (0.5 - w1)
            return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        else:
            # 20
            if b0 <= 1.:
                w = bpser(a0, b0, x0, eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            if (x0 >= 0.3):
                w1 = bpser(b0, a0, y0, eps)
                w = 0.5 + (0.5 - w1)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            if x0 >= 0.1:
                # 30
                if (b0 > 15.):
                    w1, ierr1 = bgrat(b0, a0, y0, x0, w1, 15.*eps)
                    w = 0.5 + (0.5 - w1)
                    return (w, w1, 0) if (ind == 0) else (w1, w, 0)
                else:
                    # 140
                    n = 20
                    w1 = bup(b0, a0, y0, x0, n, eps)
                    b0 += n
                    w1, ierr1 = bgrat(b0, a0, y0, x0, w1, 15.*eps)
                    w = 0.5 + (0.5 - w1)
                    return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            if (x0*b0)**a0 <= 0.7:
                w = bpser(a0, b0, x0, eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            if (b0 > 15.):
                w1, ierr1 = bgrat(b0, a0, y0, x0, w1, 15.*eps)
                w = 0.5 + (0.5 - w1)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

            # 140
            n = 20
            w1 = bup(b0, a0, y0, x0, n, eps)
            b0 += n
            w1, ierr1 = bgrat(b0, a0, y0, x0, w1, 15.*eps)
            w = 0.5 + (0.5 - w1)
            return (w, w1, 0) if (ind == 0) else (w1, w, 0)
    else:
        lmbda = ((a + b)*y - b) if (a > b) else (a - (a + b)*x)
        if lmbda < 0.:
            ind = 1
            a0, b0 = b, a
            x0, y0 = y, x
            lmbda = abs(lmbda)

        if (b0 < 40.) and (b0*x0 <= 0.7):
            w = bpser(a0, b0, x0, eps)
            w1 = 0.5 + (0.5 - w)
            return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        if (b0 < 40.):
            n = <int>b0
            b0 -= n
            if (b0 == 0.):
                n -= 1
                b0 = 1.

            w = bup(b0, a0, y0, x0, n, eps)
            if x0 <= 0.7:
                w += bpser(a0,b0,x0,eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)
            else:
                if a0 <= 15.:
                    n = 20
                    w += bup(a0, b0, x0, y0, n, eps)
                    a0 += n
                w, ierr1 = bgrat(a0, b0, x0, y0, w, 15.*eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        if a0 > b0:
            if (b0 <= 100) or (lmbda > 0.03*b0):
                w = bfrac(a0, b0, x0, y0, lmbda, 15.0*eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)
            else:
                # 200
                w = basym(a0, b0, lmbda, 100.*eps)
                w1 = 0.5 + (0.5 - w)
                return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        if (a0 <= 100.) or (lmbda > 0.03*a0):
            w = bfrac(a0, b0, x0, y0, lmbda, 15.0*eps)
            w1 = 0.5 + (0.5 - w)
            return (w, w1, 0) if (ind == 0) else (w1, w, 0)

        #200
        w = basym(a0, b0, lmbda, 100.*eps)
        w1 = 0.5 + (0.5 - w)
        return (w, w1, 0) if (ind == 0) else (w1, w, 0)

# %%----------------------------------------- brcmp1
cdef inline double brcmp1(int mu, double a, double b,
                          double x, double y) noexcept nogil:
    cdef double a0,apb,b0,c,e,h,lmbda,lnx,lny,t,u,v,x0,y0,z
    cdef int i, n
    cdef double const = 1./sqrt(2.*PI)

    a0 = min(a, b)
    if a0 >= 8.:
        if a > b:
            h = b / a
            x0 = 1. / (1. + h)
            y0 = h / (1. + h)
            lmbda = (a + b)*y - b
        else:
            h = a / b
            x0 = h / (1. + h)
            y0 = 1. / (1. + h)
            lmbda = a - (a + b)*x

        e = -lmbda / a
        if abs(e) > 0.6:
            u = e - log(x / x0)
        else:
            u = rlog1(e)

        e = lmbda / b
        if abs(e) > 0.6:
            v = e - log(y / y0)
        else:
            v = rlog1(e)

        z = esum(mu, -(a*u + b*v))
        return const*sqrt(b*x0)*z*exp(-bcorr(a, b))

    if x > 0.375 and y > 0.375:
        lnx = log(x)
        lny = log(y)
    elif x > 0.375 and y <= 0.375:
        lnx = alnrel(-y)
        lny = log(y)
    elif x <= 0.375:
        lnx = log(x)
        lny = alnrel(-x)

    z = a*lnx + b*lny

    if a0 >= 1.:
        z -= betaln(a, b)
        return esum(mu, z)

    b0 = max(a, b)
    if b0 >= 8.:
        u = gamln1(a0) + algdiv(a0, b0)
        return a0*esum(mu, z - u)

    if b0 > 1.:
        u = gamln1(a0)
        n = <int>(b0 - 1.)
        if n >= 1:
            c = 1.
            for i in range(n):
                b0 -= 1.
                c *= b0 / (a0 + b0)
            u += log(c)

        z -= u
        b0 -= 1.
        apb = a0 + b0

        if apb > 1.:
            u = a0 + b0 - 1.
            t = (1. + gam1(u)) / apb
        else:
            t = 1. + gam1(apb)
        return a0*esum(mu, z)*(1. + gam1(b0)) / t

    if esum(mu, z) == 0.:
        return 0.

    apb = a + b
    t = exp(z)
    if apb > 1.:
        u = a + b - 1.
        z = (1. + gam1(u)) / apb
    else:
        z = 1. + gam1(apb)

    c = (1. + gam1(a)) * (1. + gam1(b)) / z
    return t*(a0*c) / (1. + a0 / b0)

# %%----------------------------------------- brcomp
cdef inline double brcomp(double a, double b, double x, double y) noexcept nogil:
    cdef double a0, apb, b0, c, e, h, lmbda, lnx, lny, t, u, v, x0, y0, z
    cdef double const = 1. / sqrt(2 * PI)
    cdef int i, n

    if x == 0. or y == 0.:
        return 0.

    a0 = min(a, b)

    if a0 >= 8.:
        if a > b:
            h = b / a
            x0 = 1. / (1. + h)
            y0 = h / (1. + h)
            lmbda = (a + b)*y - b
        else:
            h = a / b
            x0 = h / (1. + h)
            y0 = 1. / (1. + h)
            lmbda = a - (a + b)*x

        e = -lmbda / a
        if abs(e) > 0.6:
            u = e - log(x / x0)
        else:
            u = rlog1(e)

        e = lmbda / b
        if abs(e) > 0.6:
            v = e - log(y / y0)
        else:
            v = rlog1(e)

        z = exp(-(a*u + b*v))
        return const*sqrt(b*x0)*z*exp(-bcorr(a, b))

    if x <= 0.375:
        lnx = log(x)
        lny = alnrel(-x)
    else:
        lnx = log(x) if (y > 0.375) else alnrel(-y)
        lny = log(y)

    z = a*lnx + b*lny
    if a0 >= 1.:
        z -= betaln(a, b)
        return exp(z)

    b0 = max(a, b)
    if b0 >= 8.:
        u = gamln1(a0) + algdiv(a0, b0)
        return a0*exp(z - u)

    if b0 > 1.:
        u = gamln1(a0)
        n = <int>(b0 - 1.)
        if n >= 1:
            c = 1.
            for i in range(n):
                b0 -= 1.
                c *= b0 / (a0 + b0)
            u += log(c)

        z -= u
        b0 -= 1.
        apb = a0 + b0

        if apb > 1.:
            u = a0 + b0 - 1.
            t = (1. + gam1(u)) / apb
        else:
            t = 1. + gam1(apb)
        return a0*exp(z)*(1. + gam1(b0)) / t

    if exp(z) == 0.:
        return 0.

    apb = a + b
    t = exp(z)
    if apb > 1.:
        u = a + b - 1.
        z = (1. + gam1(u)) / apb
    else:
        z = 1. + gam1(apb)

    c = (1. + gam1(a)) * (1. + gam1(b)) / z
    return t * (a0*c) / (1. + a0 / b0)

# %%----------------------------------------- bup
cdef inline double bup(double a, double b, double x, double y,
                       int n, double eps) noexcept nogil:
    cdef:
        double apb = a + b
        double ap1 = a + 1.
        double d = 1.
        double r, t, w, result
        int i, nm1
        int k = 0
        int mu = 0

    if n != 1 and a >= 1.:
        if (apb >= 1.1*ap1):
            mu = 708
            t = mu
            d = exp(-708)

    result = brcmp1(mu, a, b, x, y) / a
    if (n == 1) or (result == 0.):
        return result
    nm1 = n - 1
    w = d

    k = 0
    if b <= 1.:
        # 50
        for i in range(n-1):
            d *= ((apb + i) / (ap1 + i))*x
            w += d
            if d <= eps*w:
                break
        return result*w

    if y >= 1e-4:
        # 20
        r = (b - 1.)*x/y - a
        if r < 1.:
            # 50
            for i in range(n-1):
                d *= ((apb + i) / (ap1 + i))*x
                w += d
                if d <= eps*w:
                    break
            return result*w

        k = nm1
        t = nm1
        if r < t:
            k = <int>r
    else:
        k = nm1

    # 30
    for i in range(k):
        d *= ((apb + i) / (ap1 + i))*x
        w += d

    if k == nm1:
        return result*w

    # 50
    for i in range(k, n-1):
        d *= ((apb + i) / (ap1 + i))*x
        w += d
        if d <= eps*w:
            break

    return result*w

# %% ---------------------------------------- cdfbet_whichX
cdef inline (double, double, int, double) cdfbet_which1(
    double x, double y,double a, double b
    ) noexcept nogil:
    cdef double p, q

    if x < 0.:
        return (0., 0., -1, 0.)
    if x > 1.:
        return (0., 0., -1, 1.)
    if y < 0:
        return (0., 0., -2, 0.)
    if y > 1.:
        return (0., 0., -2, 1.)
    if a < 0.:
        return (0., 0., -3, 0.)
    if b < 0.:
        return (0., 0., -4, 0.)

    if ((abs(x+y)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 0., 4, (0. if (x+y) < 0 else 1.))

    p, q = cumbet(x, y, a, b)
    return (p, q, 0, 0.)

cdef inline (double, double, int, double) cdfbet_which2(
    double p, double q, double a, double b) noexcept nogil:

    cdef double ccum, cum, xx, yy
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq
    # Cython doesn't allow for default values in structs
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0)
    DZ.xlo = 0.
    DZ.xhi = 1.
    DZ.atol = atol
    DZ.rtol = tol
    DZ.x = 0.
    DZ.b = 0.

    if p < 0.:
        return (0., 0., -1, 0.)
    if p > 1.:
        return (0., 0., -1, 1.)
    if q < 0.:
        return (0., 0., -2, 0.)
    if q > 1.:
        return (0., 0., -2, 1.)
    if a < 0.:
        return (0., 0., -3, 0.)
    if b < 0.:
        return (0., 0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 0., 3, (0. if (p+q) < 0 else 1.))

    qporq = p <= q
    if qporq:
        dzror(&DZ)
        yy = 1. - DZ.x
        while DZ.status == 1:
            cum, ccum = cumbet(DZ.x, yy, a, b)
            DZ.fx = cum - p
            dzror(&DZ)
            yy = 1. - DZ.x
        xx = DZ.x
    else:
        dzror(&DZ)
        xx = 1. - DZ.x
        while DZ.status == 1:
            cum, ccum = cumbet(xx, DZ.x, a, b)
            DZ.fx = ccum - q
            dzror(&DZ)
            xx = 1. - DZ.x
        yy = DZ.x

    if DZ.status == -1:
        return (xx, yy, 1 if DZ.qleft else 2, 0. if DZ.qleft else 1.)
    else:
        return (xx, yy, 0, 0.)

cdef inline (double, int, double) cdfbet_which3(
    double p, double q, double x, double y, double b) noexcept nogil:
    cdef double tol = 1e-8
    cdef double atol = 1e-10
    cdef bint qporq = p <= q

    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1:
        return (0., -2, 1.)
    if x < 0.:
        return (0., -3, 0.)
    if x > 1.:
        return (0., -3, 1.)
    if y < 0:
        return (0., -4, 0.)
    if y > 1.:
        return (0., -4, 1.)
    if b < 0.:
        return (0., -5, 0.)

    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))
    if ((abs(x+y)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 4, (0. if (x+y) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumbet(x, y, DS.x, b)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)


cdef inline (double, int, double) cdfbet_which4(
    double p, double q, double x, double y, double a) noexcept nogil:
    cdef double tol = 1e-8
    cdef double atol = 1e-10
    cdef bint qporq = p <= q

    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1:
        return (0., -2, 1.)
    if x < 0.:
        return (0., -3, 0.)
    if x > 1.:
        return (0., -3, 1.)
    if y < 0:
        return (0., -4, 0.)
    if y > 1.:
        return (0., -4, 1.)
    if a < 0.:
        return (0., -5, 0.)

    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))
    if ((abs(x+y)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 4, (0. if (x+y) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumbet(x, y, a, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)


    if DS.status == -1:
        return (DS.x, (1 if DS.qleft else 2), (0. if DS.qleft else INFINITY))
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdfbin_whichX
cdef inline (double, double, int, double) cdfbin_which1(
    double s, double xn, double pr, double ompr) noexcept nogil:
    cdef double p, q
    if xn < 0.:
        return (0., 0., -2, 0.)
    if (s < 0.):
        return (0., 0., -1, 0.)
    if (s > xn):
        return (0., 0., -1, xn)
    if pr < 0.:
        return (0., 0., -3, 0.)
    if pr > 1.:
        return (0., 0., -3, 1.)
    if ompr < 0.:
        return (0., 0., -4, 0.)
    if ompr > 1.:
        return (0., 0., -4, 1.)
    if ((abs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 0., 4, (0. if (pr+ompr) < 0 else 1.))

    p, q = cumbin(s, xn, pr, ompr)
    return (p, q, 0, 0.)


cdef inline (double, int, double) cdfbin_which2(
    double p, double q, double xn, double pr, double ompr) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = xn
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = xn/2.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (xn < 0.):
        return (0., -3, 0.)
    if pr < 0.:
        return (0., -4, 0.)
    if pr > 1.:
        return (0., -4, 1.)
    if ompr < 0.:
        return (0., -5, 0.)
    if ompr > 1.:
        return (0., -5, 1.)

    if ((abs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 4, (0. if (pr+ompr) < 0 else 1.))
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumbin(DS.x, xn, pr, ompr)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else xn)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdfbin_which3(
    double p, double q, double s, double pr, double ompr) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1.e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (s < 0.):
        return (0., -3, 0.)
    if pr < 0.:
        return (0., -4, 0.)
    if pr > 1.:
        return (0., -4, 1.)
    if ompr < 0.:
        return (0., -5, 0.)
    if ompr > 1.:
        return (0., -5, 1.)

    if ((abs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 4, (0. if (pr+ompr) < 0 else 1.))
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumbin(s, DS.x, pr, ompr)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, double, int, double) cdfbin_which4(
    double p, double q, double s, double xn) noexcept nogil:
    cdef double ccum, cum, pr, ompr
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    # Cython doesn't allow for default values in structs
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0)
    DZ.xlo = 0.
    DZ.xhi = 1.
    DZ.atol = atol
    DZ.rtol = tol
    DZ.x = 0.
    DZ.b = 0.

    if p < 0.:
        return (0., 0., -1, 0.)
    if p > 1.:
        return (0., 0., -1, 1.)
    if q < 0.:
        return (0., 0., -2, 0.)
    if q > 1.:
        return (0., 0., -2, 1.)
    if (s < 0.):
        return (0., 0., -3, 0.)
    if (s > xn):
        return (0., 0., -3, xn)
    if (xn < 0.):
        return (0., 0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 0., 4, (0. if (p+q) < 0 else 1.))

    if qporq:
        dzror(&DZ)
        ompr = 1. - DZ.x
        while DZ.status == 1:
            cum, ccum = cumbin(s, xn, DZ.x, ompr)
            DZ.fx = cum - p
            dzror(&DZ)
            ompr = 1. - DZ.x
        pr = DZ.x
    else:
        dzror(&DZ)
        pr = 1. - DZ.x
        while DZ.status == 1:
            cum, ccum = cumbin(s, xn, pr, DZ.x)
            DZ.fx = ccum - q
            dzror(&DZ)
            pr = 1. - DZ.x
        ompr = DZ.x

    if DZ.status == -1:
        return (pr, ompr, 1 if DZ.qleft else 2, 0. if DZ.qleft else 1.)
    else:
        return (pr, ompr, 0, 0.)

# %% ---------------------------------------- cdfchi_whichX
cdef inline (double, double, int, double) cdfchi_which1(double x, double df) noexcept nogil:
    cdef double p, q
    if x < 0.:
        return (0., 0., -1, 0.)
    if df < 0.:
        return (0., 0., -2, 0.)
    p, q = cumchi(x, df)
    return (p, q, 0, 0)

cdef inline (double, int, double) cdfchi_which2(double p, double q, double df) noexcept nogil:
    cdef bint qporq = p <= q
    cdef double porq = p if qporq else q
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if df < 0.:
        return (0., -3, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    DS.small = 0.
    DS.big = 1.e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumchi(DS.x, df)
        DS.fx = cum - p if qporq else ccum - q
        if DS.fx + porq <= 1.5:
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdfchi_which3(double p, double q, double x) noexcept nogil:
    cdef bint qporq = p <= q
    cdef double porq = p if qporq else q
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if x < 0.:
        return (0., -3, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    DS.small = 0.
    DS.big = 1.e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumchi(x, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        if DS.fx + porq <= 1.5:
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdfchn_whichX
cdef inline (double, double, int, double) cdfchn_which1(
    double x, double df, double pnonc) noexcept nogil:
    cdef double p, q
    if x < 0.:
        return (0., 0., -1, 0.)
    if df < 0.:
        return (0., 0., -2, 0.)
    if pnonc < 0.:
        return (0., 0., -3, 0.)
    p, q = cumchn(x, df, pnonc)
    return (p, q, 0, 0.)

cdef inline (double, int, double) cdfchn_which2(
    double p, double df, double pnonc) noexcept nogil:
    cdef double cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    df = min(df, spmpar[2])
    pnonc = min(pnonc, 1.e9)

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if df < 0.:
        return (0., -2, 0.)
    if pnonc < 0.:
        return (0., -3, 0.)

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, _ = cumchn(DS.x, df, pnonc)
        DS.fx = cum - p
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdfchn_which3(
    double p, double x, double pnonc) noexcept nogil:
    cdef double cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    x = min(x, spmpar[2])
    pnonc = min(pnonc, 1.e9)

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if x < 0.:
        return (0., -2, 0.)
    if pnonc < 0.:
        return (0., -3, 0.)

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, _ = cumchn(x, DS.x, pnonc)
        DS.fx = cum - p
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdfchn_which4(
    double p, double x, double df) noexcept nogil:
    cdef double cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1.e9
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    x = min(x, spmpar[2])
    df = min(df, spmpar[2])

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if x < 0.:
        return (0., -2, 0.)
    if df < 0.:
        return (0., -3, 0.)

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, _ = cumchn(x, df, DS.x)
        DS.fx = cum - p
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdff_whichX
cdef inline (double, double, int, double) cdff_which1(
    double f, double dfn, double dfd) noexcept nogil:
    cdef double p, q
    if f < 0.:
        return (0., 0., -1, 0.)
    if dfn < 0.:
        return (0., 0., -2, 0.)
    if dfd < 0.:
        return (0., 0., -3, 0.)
    p, q = cumf(f, dfn, dfd)
    return (p, q, 0, 0.)

cdef inline (double, int, double) cdff_which2(
    double p, double q, double dfn, double dfd) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if dfn < 0.:
        return (0., -3, 0.)
    if dfd < 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumf(DS.x, dfd, dfd)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdff_which3(
    double p, double q, double f, double dfd) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if f < 0.:
        return (0., -3, 0.)
    if dfd < 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumf(f, DS.x, dfd)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdff_which4(
    double p, double q, double f, double dfn) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if f < 0.:
        return (0., -3, 0.)
    if dfn < 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumf(f, dfn, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdffnc_whichX
cdef (double, double, int, double) cdffnc_which1(
    double f, double dfn, double dfd, double phonc) noexcept nogil:
    cdef double p, q
    cdef int ierr
    if f < 0.:
        return (0., 0., -1, 0.)
    if dfn < 0.:
        return (0., 0., -2, 0.)
    if dfd < 0.:
        return (0., 0., -3, 0.)
    if phonc < 0.:
        return (0., 0., -4, 0.)
    p, q, ierr = cumfnc(f, dfn, dfd, phonc)
    if ierr != 0:
        return (p, q, 10, 0.)
    else:
        return (p, q, 0, 0.)

cdef (double, int, double) cdffnc_which2(
    double p, double q, double dfn, double dfd, double phonc) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef int ierr
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if dfn < 0.:
        return (0., -3, 0.)
    if dfd < 0.:
        return (0., -4, 0.)
    if phonc < 0.:
        return (0., -5, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum, ierr = cumfnc(DS.x, dfn, dfd, phonc)
        DS.fx = cum - p if qporq else ccum - q
        if ierr != 0:
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef (double, int, double) cdffnc_which3(
    double p, double q, double f, double dfd, double phonc) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef int ierr
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if f < 0.:
        return (0., -3, 0.)
    if dfd < 0.:
        return (0., -4, 0.)
    if phonc < 0.:
        return (0., -5, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum, ierr = cumfnc(f, DS.x, dfd, phonc)
        DS.fx = cum - p if qporq else ccum - q
        if ierr != 0:
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef (double, int, double) cdffnc_which4(
    double p, double q, double f, double dfn, double phonc) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef int ierr
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if f < 0.:
        return (0., -3, 0.)
    if dfn < 0.:
        return (0., -4, 0.)
    if phonc < 0.:
        return (0., -5, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum, ierr = cumfnc(f, dfn, DS.x, phonc)
        DS.fx = cum - p if qporq else ccum - q
        if ierr != 0:
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef (double, int, double) cdffnc_which5(
    double p, double q, double f, double dfn, double dfd) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef int ierr
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1.e9
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if f < 0.:
        return (0., -3, 0.)
    if dfn < 0.:
        return (0., -4, 0.)
    if dfd < 0.:
        return (0., -5, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum, ierr = cumfnc(f, dfn, dfd, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        if ierr != 0:
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdfgam_whichX
cdef (double, double, int, double) cdfgam_which1(
    double x, double shape, double scale) noexcept nogil:
    cdef double p, q

    if x < 0.:
        return (0., 0., -1, 0.)
    if shape < 0.:
        return (0., 0., -2, 0.)
    if scale < 0.:
        return (0., 0., -3, 0.)

    p, q = cumgam(x*scale, shape)
    if p >= 1.5:
        return (p, q, 10, 0.)
    else:
        return (p, q, 0, 0.)

cdef (double, int, double) cdfgam_which2(
    double p, double q, double shape, double scale) noexcept nogil:
    cdef double xx
    cdef int ierr
    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if shape <= 0.:
        return (0., -3, 0.)
    if scale <= 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    xx, ierr = gaminv(shape, p, q, -1)
    if ierr < 0:
        return (0., 10, 0.)
    else:
        return (xx/scale, 0, 0.)

cdef (double, int, double) cdfgam_which3(
    double p, double q, double x, double scale) noexcept nogil:
    cdef double cum, ccum
    cdef double xscale = x*scale
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q

    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if x <= 0.:
        return (0., -3, 0.)
    if scale <= 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumgam(xscale, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        if (qporq and (cum > 1.5) or ((not qporq) and (ccum > 1.5))):
            return (DS.x, 10, 0.)
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)


cdef (double, int, double) cdfgam_which4(
    double p, double q, double x, double shape) noexcept nogil:
    cdef double xx
    cdef int ierr

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if x <= 0.:
        return (0., -3, 0.)
    if shape <= 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    xx, ierr = gaminv(shape, p, q, -1)
    if ierr < 0:
        return (0., 10, 0.)
    else:
        return (xx/x, 0, 0.)

# %% ---------------------------------------- cdfnbn_whichX
cdef inline (double, double, int, double) cdfnbn_which1(
    double s, double xn, double pr, double ompr) noexcept nogil:
    cdef double p, q
    if xn < 0.:
        return (0., 0., -2, 0.)
    if (s < 0.):
        return (0., 0., -1, 0.)
    if (s > xn):
        return (0., 0., -1, xn)
    if pr < 0.:
        return (0., 0., -3, 0.)
    if pr > 1.:
        return (0., 0., -3, 1.)
    if ompr < 0.:
        return (0., 0., -4, 0.)
    if ompr > 1.:
        return (0., 0., -4, 1.)
    if ((abs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 0., 4, (0. if (pr+ompr) < 0 else 1.))

    p, q = cumnbn(s, xn, pr, ompr)
    return (p, q, 0, 0.)


cdef inline (double, int, double) cdfnbn_which2(
    double p, double q, double xn, double pr, double ompr) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (xn < 0.):
        return (0., -3, 0.)
    if pr < 0.:
        return (0., -4, 0.)
    if pr > 1.:
        return (0., -4, 1.)
    if ompr < 0.:
        return (0., -5, 0.)
    if ompr > 1.:
        return (0., -5, 1.)

    if ((abs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 4, (0. if (pr+ompr) < 0 else 1.))
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumnbn(DS.x, xn, pr, ompr)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else xn)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdfnbn_which3(
    double p, double q, double s, double pr, double ompr) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (s < 0.):
        return (0., -3, 0.)
    if pr < 0.:
        return (0., -4, 0.)
    if pr > 1.:
        return (0., -4, 1.)
    if ompr < 0.:
        return (0., -5, 0.)
    if ompr > 1.:
        return (0., -5, 1.)

    if ((abs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 4, (0. if (pr+ompr) < 0 else 1.))
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumbin(s, DS.x, pr, ompr)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, double, int, double) cdfnbn_which4(
    double p, double q, double s, double xn) noexcept nogil:
    cdef double ccum, cum, pr, ompr
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    # Cython doesn't allow for default values in structs
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0)
    DZ.xlo = 0.
    DZ.xhi = 1.
    DZ.atol = atol
    DZ.rtol = tol
    DZ.x = 0.
    DZ.b = 0.

    if p < 0.:
        return (0., 0., -1, 0.)
    if p > 1.:
        return (0., 0., -1, 1.)
    if q < 0.:
        return (0., 0., -2, 0.)
    if q > 1.:
        return (0., 0., -2, 1.)
    if (s < 0.):
        return (0., 0., -3, 0.)
    if (xn < 0.):
        return (0., 0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 0., 4, (0. if (p+q) < 0 else 1.))

    if qporq:
        dzror(&DZ)
        ompr = 1. - DZ.x
        while DZ.status == 1:
            cum, ccum = cumbin(s, xn, DZ.x, ompr)
            DZ.fx = cum - p
            dzror(&DZ)
            ompr = 1. - DZ.x
        pr = DZ.x
    else:
        dzror(&DZ)
        pr = 1. - DZ.x
        while DZ.status == 1:
            cum, ccum = cumbin(s, xn, pr, DZ.x)
            DZ.fx = ccum - q
            dzror(&DZ)
            pr = 1. - DZ.x
        ompr = DZ.x

    if DZ.status == -1:
        return (pr, ompr, 1 if DZ.qleft else 2, 0. if DZ.qleft else 1.)
    else:
        return (pr, ompr, 0, 0.)

# %% ---------------------------------------- cdfnor_whichX
cdef inline (double, double, int, double) cdfnor_which1(
    double x, double mean, double sd) noexcept nogil:
    cdef double z, p, q
    if sd <= 0.:
        return (0., 0., -3, 0.)
    z = (x-mean)/sd
    p, q = cumnor(z)
    return (p, q, 0, 0.)

cdef inline (double, int, double) cdfnor_which2(
    double p, double q, double mean, double sd) noexcept nogil:
    cdef double z
    if sd <= 0.:
        return (0., -4, 0.)
    z = dinvnr(p, q)
    return (sd*z + mean, 0, 0.)

cdef inline (double, int, double) cdfnor_which3(
    double p, double q, double x, double sd) noexcept nogil:
    cdef double z
    if sd <= 0.:
        return (0., -4, 0.)
    z = dinvnr(p, q)
    return (x - sd*z, 0, 0.)

cdef inline (double, int, double) cdfnor_which4(
    double p, double q, double x, double mean) noexcept nogil:
    cdef double z
    z = dinvnr(p, q)
    return ((x-mean)/z, 0, 0.)

# %% ---------------------------------------- cdfpoi_whichX
cdef inline (double, double, int, double) cdfpoi_which1(
    double s, double xlam) noexcept nogil:
    cdef double p, q
    if s < 0.:
        return (0., 0., -1, 0.)
    if xlam < 0.:
        return (0., 0., -2, 0.)
    p, q = cumpoi(s, xlam)
    return (p, q, 0, 0.)

cdef inline (double, int, double) cdfpoi_which2(
    double p, double q, double xlam) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (xlam < 0.):
        return (0., -3, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    if (xlam < 0.01) and (p < 0.975):
        return (0., 0, 0.)

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumpoi(DS.x, xlam)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdfpoi_which3(
    double p, double q, double s) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e300
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (s < 0.):
        return (0., -3, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumpoi(s, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdft_whichX
cdef inline (double, double, int, double) cdft_which1(
    double t, double df) noexcept nogil:
    cdef double p, q
    if df <= 0.:
        return (0., 0., -2, 0.)
    p, q = cumt(t, df)
    return (p, q, 0, 0.)

cdef inline (double, int, double) cdft_which2(
    double p, double q, double df) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = -1.e100
    DS.big = 1.e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = dt1(p, q, df)

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if (df < 0.):
        return (0., -3, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumt(DS.x, df)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, -1.e100 if DS.qleft else 1.e100)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdft_which3(
    double p, double q, double t) noexcept nogil:
    cdef double ccum, cum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)
    DS.small = 0.
    DS.big = 1e10
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumt(t, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cdftnc_whichX
cdef inline (double, double, int, double) cdftnc_which1(
    double t, double df, double pnonc) noexcept nogil:
    cdef double p, q
    if df < 0.:
        return (0., 0., -2, 0.)

    df = min(df, 1.e10)
    t = max(min(t, spmpar[2]), -spmpar[2])

    if t != t:
        return (0., 0., -1, 0.)
    if not (-1.e-6 <= pnonc <= 1.e6):
        return (0., 0., -3, 1.e6 if pnonc > -1e6 else -1.e6)


    p, q = cumtnc(t, df, pnonc)
    return (p, q, 0, 0.)

cdef inline (double, int, double) cdftnc_which2(
    double p, double q, double df, double pnonc) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = -1.e100
    DS.big = 1.e100
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    df = min(df, 1.e10)
    if not (-1.e-6 <= pnonc <= 1.e6):
        return (0., -4, 1.e6 if pnonc > -1e6 else -1.e6)
    if df < 0.:
        return (0., -3, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumf(DS.x, df, pnonc)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdftnc_which3(
    double p, double q, double t, double pnonc) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = 0
    DS.big = 1.e6
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    t = max(min(t, spmpar[2]), -spmpar[2])
    if t != t:
        return (0., -3, 0.)
    if not (-1.e-6 <= pnonc <= 1.e6):
        return (0., -4, 1.e6 if pnonc > -1e6 else -1.e6)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumf(t, DS.x, pnonc)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

cdef inline (double, int, double) cdftnc_which4(
    double p, double q, double t, double df) noexcept nogil:
    cdef double cum, ccum
    cdef double tol = 1e-8
    cdef double atol = 1e-50
    cdef bint qporq = p <= q
    cdef DinvrState DS = DinvrState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    cdef DzrorState DZ = DzrorState(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0)

    DS.small = -1.e6
    DS.big = 1.e6
    DS.absstp = 0.5
    DS.relstp = 0.5
    DS.stpmul = 5.
    DS.abstol = atol
    DS.reltol = tol
    DS.x = 5.

    if p < 0.:
        return (0., -1, 0.)
    if p > 1.:
        return (0., -1, 1.)
    if q < 0.:
        return (0., -2, 0.)
    if q > 1.:
        return (0., -2, 1.)
    if t < 0.:
        return (0., -3, 0.)
    if df < 0.:
        return (0., -4, 0.)
    if ((abs(p+q)-0.5)-0.5) > 3*spmpar[0]:
        return (0., 3, (0. if (p+q) < 0 else 1.))

    t = max(min(t, spmpar[2]), -spmpar[2])
    df = min(df, 1.e10)
    if t != t:
        return (0., -3, 0.)
    if df < 0.:
        return (0., -4, 0.)

    dinvr(&DS, &DZ)
    while DS.status == 1:
        cum, ccum = cumf(t, df, DS.x)
        DS.fx = cum - p if qporq else ccum - q
        dinvr(&DS, &DZ)

    if DS.status == -1:
        return (DS.x, 1 if DS.qleft else 2, 0. if DS.qleft else INFINITY)
    else:
        return(DS.x, 0, 0.)

# %% ---------------------------------------- cumbet
cdef inline (double, double) cumbet(double x, double y,
                                    double a, double b) noexcept nogil:
    cdef double u, v
    cdef int ierr

    if x <= 0.:
        return (0., 1.)

    if y <= 0.:
        return (1., 0.)
    u, v, _ = bratio(a, b, x, y)
    return (u, v)

# %% ---------------------------------------- cumbin
cdef inline (double, double) cumbin(double s, double xn,
                                    double pr, double ompr) noexcept nogil:
    cdef double cum, ccum
    if s < xn:
        ccum, cum = cumbet(pr, ompr, s + 1., xn - s)
    else:
        cum, ccum = 1., 0.
    return ccum, cum

# %% ---------------------------------------- cumchi
cdef inline (double, double) cumchi(double x, double df) noexcept nogil:
    return cumgam(0.5*x, 0.5*df)

# %% ---------------------------------------- cumchn
cdef inline (double, double) cumchn(double x, double df, double pnonc) noexcept nogil:
    cdef double adj, centaj, centwt, chid2, dfd2, lcntaj, lcntwt
    cdef double lfact, pcent, pterm, ssum, sumadj, term, wt, xnonc
    cdef double eps = 1.e-15
    cdef double abstol = 1.e-300
    cdef int i, icent

    if x <= 0.:
        return (0., 1.)
    if pnonc <= 1e-10:
        return cumchi(x, df)
    xnonc = pnonc/2.
    icent = <int>xnonc
    if icent == 0:
        icent = 1
    chid2 = x / 2.
    lfact = alngam(icent + 1)
    lcntwt = -xnonc + icent*log(xnonc) - lfact
    centwt = exp(lcntwt)
    pcent, _ = cumchi(x, df + 2.*icent)
    dfd2 = (df + 2.*icent)/2.
    lfact = alngam(1. + dfd2)
    lcntaj = dfd2*log(chid2) - chid2 - lfact
    centaj = exp(lcntaj)
    ssum = centwt*pcent
    sumadj = 0.
    adj = centaj
    wt = centwt
    i = icent

    while True:
        dfd2 = (df + 2.*i)/2.
        adj *= dfd2 / chid2
        sumadj += adj
        pterm = pcent + sumadj
        wt *= (i / xnonc)
        term = wt*pterm
        ssum += term
        i -= 1
        if ((ssum < abstol) or (term < eps*ssum)) or (i == 0):
            break

    sumadj, adj, wt, i = centaj, centaj, centwt, icent

    while True:
        wt *= (xnonc / (i+1.))
        pterm = pcent - sumadj
        term = wt*pterm
        ssum += term
        i += 1
        dfd2 = (df + 2.*i)/2.
        adj *= chid2/dfd2
        sumadj += adj
        if ((ssum < abstol) or (term < eps*ssum)):
            break

    return (ssum, 0.5 + (0.5 - ssum))

# %% ---------------------------------------- cumf
cdef inline (double, double) cumf(double f, double dfn, double dfd) noexcept nogil:
    cdef double dsum,prod,xx,yy
    cdef double cum, ccum
    cdef int ierr
    if f <= 0.:
        return (0., 1.)
    prod = dfn*f
    dsum = dfd + prod
    xx = dfd / dsum
    if xx > 0.5:
        yy = prod / dsum
        xx = 1. - yy
    else:
        yy = 1. - xx

    cum, ccum, ierr = bratio(dfd*0.5, dfn*0.5, xx, yy)
    return cum, ccum

# %% ---------------------------------------- cumfnc
cdef inline (double, double, int) cumfnc(double f, double dfn,
                                    double dfd, double pnonc) noexcept nogil:
    cdef double dsum, prod, xx, yy, adn, aup, b
    cdef double betdn, betup, centwt, dnterm, ssum
    cdef double upterm, xmult, xnonc
    cdef double cum, ccum
    cdef double eps = 1e-4
    cdef double abstol = 1e-300
    cdef int status = 0
    cdef int i, icent

    if f <= 0.:
        return 0., 1., status

    if pnonc < 1e-10:
        cum, ccum = cumf(f, dfn, dfd)
        return cum, ccum, status

    xnonc = pnonc / 2.
    icent = <int>xnonc
    if abs(xnonc - icent) >= 1.:
        return (0., 0., 1)
    if icent == 0:
        icent = 1

    centwt = exp(-xnonc + icent*log(xnonc) - alngam(icent + 1))
    prod = dfn * f
    dsum = dfd + prod
    yy = dfd / dsum
    if yy > 0.5:
        xx = prod / dsum
        yy = 1. - xx
    else:
        xx = 1. - yy

    betdn, _, ierr = bratio(dfn*0.5 + icent, dfd*0.5, xx, yy)
    adn = dfn/2. + icent
    aup = adn
    b = dfd / 2.
    betup = betdn
    ssum = centwt*betdn

    xmult = centwt
    i = icent

    if adn < 2.:
        dnterm = exp(alngam(adn+b) - alngam(adn+1.) - alngam(b) +
                     adn*log(xx) + b*log(yy))
    else:
        dnterm = exp(-betaln(adn, b) - log(adn) + adn*log(xx) +
                     b*log(yy))
    while True:
        xmult *= (i/xnonc)
        i -= 1
        adn -= 1
        dnterm *= (adn + 1) / ((adn + b)*xx)
        betdn += dnterm
        ssum += xmult*betdn
        if ((ssum < abstol) and (xmult*betdn < eps*ssum)) or i <= 0:
            break
    i = icent + 1
    xmult = centwt
    if (aup - 1 + b) == 0:
          upterm = exp(-alngam(aup) - alngam(b) + (aup-1)*log(xx) + b*log(yy))

    else:
        if (aup < 2):
            upterm = exp(alngam(aup-1+b) - alngam(aup)
                         - alngam(b) + (aup-1)*log(xx) + b*log(yy))
        else:
            # Same expression, but avoid problems for large aup
            upterm = exp(-betaln(aup-1, b) - log(aup-1) +
                         (aup-1)*log(xx) + b*log(yy))

    while True:
        xmult *= xnonc / i
        i += 1
        aup += 1
        upterm *= (aup + b - 2.)*xx/(aup - 1.)
        betup -= upterm
        ssum += xmult
        if ((ssum < abstol) and (xmult*betup < eps*ssum)):
            break

    return (ssum, 0.5 + (0.5 - ssum), status)

# %% ---------------------------------------- cumgam
cdef inline (double, double) cumgam(double x, double a) noexcept nogil:
    return gratio(a, x, 0) if (x > 0.) else (0., 1.)

# %% ---------------------------------------- cumnbn
cdef inline (double, double) cumnbn(double s, double xn,
                                    double pr, double ompr) noexcept nogil:
    return cumbet(pr, ompr, xn, s+1.)

# %% ---------------------------------------- cumnor
cdef inline (double, double) cumnor(double x) noexcept nogil:
    cdef double[5] a = [2.2352520354606839287e00, 1.6102823106855587881e02,
                        1.0676894854603709582e03, 1.8154981253343561249e04,
                        6.5682337918207449113e-2]
    cdef double[4] b = [4.7202581904688241870e01, 9.7609855173777669322e02,
                        1.0260932208618978205e04, 4.5507789335026729956e04]
    cdef double[9] c = [3.9894151208813466764e-1, 8.8831497943883759412e00,
                        9.3506656132177855979e01, 5.9727027639480026226e02,
                        2.4945375852903726711e03, 6.8481904505362823326e03,
                        1.1602651437647350124e04, 9.8427148383839780218e03,
                        1.0765576773720192317e-8]
    cdef double[8] d = [2.2266688044328115691e01, 2.3538790178262499861e02,
                        1.5193775994075548050e03, 6.4855582982667607550e03,
                        1.8615571640885098091e04, 3.4900952721145977266e04,
                        3.8912003286093271411e04, 1.9685429676859990727e04]
    cdef double[6] p = [2.1589853405795699e-1, 1.274011611602473639e-1,
                        2.2235277870649807e-2, 1.421619193227893466e-3,
                        2.9112874951168792e-5, 2.307344176494017303e-2]
    cdef double[5] q = [1.28426009614491121e00, 4.68238212480865118e-1,
                        6.59881378689285515e-2, 3.78239633202758244e-3,
                        7.29751555083966205e-5]
    cdef double eps = spmpar[0] * 0.5
    cdef double tiny = spmpar[1]
    cdef double y = abs(x)
    cdef double threshold = 0.66291
    cdef double dl, result, xden, xnum, xsq, ccum
    cdef int i

    if y <= threshold:
        xsq = x*x if (y > eps) else 0.
        xnum = a[4] * xsq
        xden = xsq
        for i in range(3):
            xnum += a[i]
            xnum *= xsq
            xden += b[i]
            xden *= xsq
        result = x*(xnum + a[3]) / (xden + b[3])
        ccum = 0.5 - result
        result += 0.5
    elif y < sqrt(32):
        xnum = c[8]*y
        xden = y
        for i in range(7):
            xnum += c[i]
            xnum *= y
            xden += d[i]
            xden *= y
        result = (xnum + c[7]) / (xden + d[7])
        xsq = <int>(y*1.6) / 1.6
        dl = (y - xsq) * (y + xsq)
        result *= exp(-xsq*xsq*0.5)*exp(-0.5*dl)
        ccum = 1. - result
        if x > 0:
            ccum, result = result, ccum
    else:
        result = 0.
        xsq = (1 / x)**2
        xnum = p[5] * xsq
        xden = xsq
        for i in range(4):
            xnum += p[i]
            xnum *= xsq
            xden += q[i]
            xden *= xsq
        result = xsq*(xnum + p[4]) / (xden + q[4])
        result = (sqrt(1./(2.*PI)) - result) / y
        xsq = <int>(x*1.6) / 1.6
        dl = (x - xsq) * (x + xsq)
        result *= exp(-xsq*xsq*0.5)*exp(-0.5*dl)
        ccum = 1. - result
        if x > 0:
            ccum, result = result, ccum

    if result < tiny:
        result = 0.
    if ccum < tiny:
        ccum = 0.

    return result, ccum

# %%----------------------------------------- cumpoi
cdef inline (double, double) cumpoi(double s, double xlam) noexcept nogil:
    return cumchi(2*xlam, 2.*(s + 1.))

# %%----------------------------------------- cumt
cdef inline (double, double) cumt(double t, double df) noexcept nogil:
    cdef double a, oma, tt, dfptt, xx, yy, cum, ccum

    tt = t*t
    dfptt = df + tt
    xx = df / dfptt
    yy = tt / dfptt
    a, oma = cumbet(xx, yy, 0.5*df, 0.5)
    if t > 0.:
        ccum = 0.5 * a
        cum = oma + ccum
    else:
        cum = 0.5 * a
        ccum = oma + cum
    return cum, ccum

# %%----------------------------------------- cumtnc
cdef inline (double, double) cumtnc(double t, double df, double pnonc) noexcept nogil:
    cdef double alghdf, b, bb, bbcent, bcent, cent, d, dcent, dpnonc
    cdef double dum1, dum2, e, ecent, lmbda, lnomx, lnx, omx
    cdef double pnonc2, s, scent, ss, sscent, t2, term, tt, twoi, x
    cdef double xi, xlnd, xlne
    cdef double conv = 1.e-7
    cdef double tiny = 1.e-10
    cdef int ierr
    cdef bint qrevs

    if abs(pnonc) <= tiny:
        return cumt(t, df)
    qrevs = t < 0.
    if qrevs:
        tt = -t
        dpnonc = -pnonc

    pnonc2 = dpnonc**2
    t2 = tt**2
    if abs(tt) <= tiny:
        return cumnor(-pnonc)

    lmbda = 0.5*pnonc2
    x = df / (df + t2)
    omx = 1. - x
    lnx = log(x)
    lnomx = log(omx)
    alghdf = gamln(0.5*df)
    cent = max(floor(lmbda), 1.)
    xlnd = cent*log(lmbda) - gamln(cent + 1.) - lmbda
    dcent = exp(xlnd)
    xlne = (cent + 0.5)*log(lmbda) - gamln(cent + 1.5) - lmbda
    ecent = exp(xlne)
    if dpnonc < 0.:
        ecent = -ecent
    bcent, dum1, ierr = bratio(0.5*df, cent + 0.5, x, omx)
    bbcent, dum2, ierr = bratio(0.5*df, cent + 1., x, omx)

    if (bbcent + bcent) < tiny:
        return (0., 1.) if qrevs else (1., 0.)

    if (dum1 + dum2) < tiny:
        return cumnor(-pnonc)

    ccum = dcent*bcent + ecent*bbcent
    scent = exp(gamln(0.5*df + cent + 0.5) - gamln(cent + 1.5) - alghdf
                + 0.5*df*lnx + (cent + 0.5)*lnomx)
    sscent = exp(gamln(0.5*df + cent + 1.) - gamln(cent + 2.) - alghdf
                 + 0.5*df*lnx + (cent + 1.)*lnomx)
    xi = cent + 1.
    twoi = 2.*xi
    d, e, b, bb, s, ss = dcent, ecent, bcent, bbcent, scent, sscent
    while True:
        b += s
        bb += ss
        d = (lmbda/xi)*d
        e = (lmbda/(xi + 0.5))*e
        term = d*b + e*bb
        ccum += term
        s *= omx*(df + twoi - 1.)/(twoi + 1.)
        ss *= omx*(df + twoi)/(twoi + 2.)
        xi += 1.
        twoi *= xi
        if abs(term) <= conv*ccum:
            break

    xi = cent
    twoi = 2.*xi
    d, e, b, bb = dcent, ecent, bcent, bbcent
    s = scent*(1. + twoi)/((df + twoi - 1.)*omx)
    ss = sscent*(2. + twoi)/((df + twoi)*omx)

    while True:
        b -= s
        bb -= ss
        d *= (xi / lmbda)
        e *= (xi + 0.5)/lmbda
        term = d*b + e*bb
        ccum += term
        xi -= 1.
        if xi < 0.5:
            break
        twoi *= xi
        s *= (1. + twoi) / ((df + twoi - 1.)*omx)
        ss *= (2. + twoi) / ((df + twoi)*omx)
        if abs(term) <= conv*ccum:
            break

    if qrevs:
        cum = max(min(0.5*ccum, 1.), 0.)
        ccum = max(min(1.-cum, 1.), 0.)
    else:
        ccum = max(min(0.5*ccum, 1.), 0.)
        cum = max(min(1.-ccum, 1.), 0.)
    return cum, ccum

# %%----------------------------------------- devlpl
cdef inline double devlpl(double *a, int n, double x) noexcept nogil:
    cdef double temp = a[n-1]
    cdef int i

    for i in range(n-2,-1,-1):
        temp = a[i] + temp*x
    return temp

# %%-------------------------------------- dinvnr
cdef inline double dinvnr(double p, double q) noexcept nogil:
    cdef int maxit = 100
    cdef double eps = 1e-13
    cdef double r2pi = sqrt(1. / (2.*PI))
    cdef double strtx, xcur, cum, pp, dx
    cdef int i

    pp = q if (p > q) else p
    strtx = stvaln(pp)
    xcur = strtx

    for i in range(maxit):
        cum, _ = cumnor(xcur)
        dx = (cum - pp) / (r2pi * exp(0.5*xcur*xcur))
        xcur -= dx
        if abs(dx / xcur) < eps:
            return -xcur if (p > q) else xcur

    return -strtx if (p > q) else strtx

# %% ------------------------------------- dinvr
cdef struct DinvrState:
    double absstp
    double abstol
    double big
    double fbig
    double fx
    double fsmall
    double relstp
    double reltol
    double small
    int status
    double step
    double stpmul
    double x
    double xhi
    double xlb
    double xlo
    double xsave
    double xub
    double yy
    double zx
    double zy
    double zz
    int next_state
    bint qbdd
    bint qcond
    bint qdum1
    bint qdum2
    bint qhi
    bint qleft
    bint qincr
    bint qlim
    bint qok
    bint qup

cdef void dinvr(DinvrState *S, DzrorState *DZ) noexcept nogil:
    """Main zero-finding function. If not returned, cycles
    through different states and modifies S in place. Progress
    is achieved through reverse communication and tracking
    status variable.
    """
    while True:
        if S.next_state == 0:
            S.qcond = (S.small <= S.x <= S.big)
            if not S.qcond:
                S.status = -2
                return
            S.xsave = S.x
            S.x = S.small
            S.next_state = 10
            S.status = 1
            return

        elif S.next_state == 10:
            S.fsmall = S.fx
            S.x = S.big
            S.next_state = 20
            S.status = 1
            return

        elif S.next_state == 20:
            S.fbig = S.fx
            S.qincr = S.fbig > S.fsmall
            S.status = -1
            if not S.qincr:
                # 50
                if S.fsmall >= 0.:
                    # 60
                    if S.fbig <= 0.:
                        # 70
                        pass
                    else:
                        S.qleft = False
                        S.qhi = True
                        return
                else:
                    S.qleft = True
                    S.qhi = False
                    return
            else:
                if S.fsmall <= 0.:
                    # 30
                    if S.fbig >= 0.:
                        # 40
                        pass
                    else:
                        S.qleft = False
                        S.qhi = False
                        return
                else:
                    S.qleft = True
                    S.qhi = True
                    return

            S.status = 1
            S.x = S.xsave
            S.step = max(S.absstp, S.relstp*abs(S.x))
            S.next_state = 90
            return

        elif S.next_state == 90:
            S.yy = S.fx
            if S.yy == 0.:
                S.status = 0
                S.qok = True
                return
            S.next_state = 100

        elif S.next_state == 100:
            S.qup = (S.qincr and (S.yy < 0.)) or ((not S.qincr) and (S.yy > 0.))
            if S.qup:
                S.xlb = S.xsave
                S.xub = min(S.xlb + S.step, S.big)
                S.next_state = 120
            else:
                # 170
                S.xub = S.xsave
                S.xlb = max(S.xub - S.step, S.small)
                S.next_state = 190

        elif S.next_state == 120:
            S.x = S.xub
            S.status = 1
            S.next_state = 130
            return

        elif S.next_state == 130:
            S.yy = S.fx
            S.qbdd = (S.qincr and (S.yy >= 0.)) or ((not S.qincr) and (S.yy <= 0.))
            S.qlim = (S.xub >= S.big)
            S.qcond = S.qbdd or S.qlim
            if S.qcond:
                S.next_state = 150
            else:
                S.step *= S.stpmul
                S.xlb = S.xub
                S.xub = min(S.xlb + S.step, S.big)
                S.next_state = 120

        elif S.next_state == 150:
            if S.qlim and (not S.qbdd):
                S.status = -1
                S.qleft = False
                S.qhi = not S.qincr
                S.x = S.big
                return
            else:
                S.next_state = 240

        elif S.next_state == 190:
            S.x = S.xlb
            S.status = 1
            S.next_state = 200
            return

        elif S.next_state == 200:
            S.yy = S.fx
            S.qbdd = (S.qincr and (S.yy <= 0.)) or ((not S.qincr) and (S.yy >= 0.))
            S.qlim = (S.xlb <= S.small)
            S.qcond = S.qbdd or S.qlim
            if S.qcond:
                S.next_state = 220
            else:
                S.step *= S.stpmul
                S.xub = S.xlb
                S.xlb = max(S.xub - S.step, S.small)
                S.next_state = 190

        elif S.next_state == 220:
            if S.qlim and (not S.qbdd):
                S.status = -1
                S.qleft = True
                S.qhi = S.qincr
                S.x = S.small
                return
            else:
                S.next_state = 240

        elif S.next_state == 240:
            # Overwrite supplied DZ with the problem
            DZ.xhi = S.xub
            DZ.xlo = S.xlb
            DZ.atol = S.abstol
            DZ.rtol = S.reltol
            DZ.x = S.xlb
            DZ.b = S.xlb
            S.next_state = 250

        elif S.next_state == 250:
            dzror(DZ)
            if DZ.status == 1:
                S.next_state = 260
                S.status = 1
                S.x = DZ.x
                return
            else:
                S.x = DZ.xlo
                S.status = 0
                return

        elif S.next_state == 260:
            DZ.fx = S.fx
            S.next_state = 250

        else:
            S.status = -9999  # Bad state, should not be possible to get here
            return

# %%-------------------------------------- dt1
cdef inline double dt1(double p, double q, double df) noexcept nogil:
    cdef double ssum, term, x, xx
    cdef double denpow = 1.
    cdef double[4][5] coef = [[1., 1., 0., 0., 0.],
                              [3., 16., 5., 0., 0.],
                              [-15., 17., 19., 3., 0.],
                              [-945., -1920., 1482., 776., 79.]
                             ]
    cdef double[4] denom = [4., 96., 384.0, 92160.0]
    cdef int[4] ideg = [2, 3, 4, 5]

    x = abs(dinvnr(p, q))
    xx = x*x
    ssum = 1.
    for i in range(4):
        term = (devlpl(coef[i], ideg[i], xx))*x
        denpow *= df
        ssum += term / (denpow*denom[i])

    return -ssum if not (p >= 0.5) else ssum

# %% ------------------------------------- dzror
cdef struct DzrorState:
    double a
    double atol
    double b
    double c
    double d
    double fa
    double fb
    double fc
    double fd
    double fda
    double fdb
    double fx
    double m
    double mb
    double p
    double q
    double tol
    double rtol
    double w
    double xhi
    double xlo
    double x
    int ext
    int status
    int next_state
    bint first
    bint qrzero
    bint qleft
    bint qhi

cdef void dzror(DzrorState *S) noexcept nogil:
    """Main zero-finding function. If not returned, cycles
    through different states and modifies S in place. Progress
    is achieved through reverse communication and tracking
    status variable.
    """
    while True:
        if S.next_state == 0:
            S.next_state = 10
            S.status = 1
            return

        elif S.next_state == 10:
            S.fb = S.fx
            S.xlo = S.xhi
            S.a = S.xlo
            S.x = S.xlo
            S.next_state = 20
            S.status = 1
            return

        elif S.next_state == 20:
            if (not (S.fb > 0.)) or (not (S.fx > 0.)):
                if (not (S.fb < 0.)) or (not (S.fx < 0.)):
        # 60
                    S.fa = S.fx
                    S.first = True
        # 70
                    S.c = S.a
                    S.fc = S.fa
                    S.ext = 0
                    S.status = 1
                    S.next_state = 80
                else:
        # 40
                    S.status = -1
                    S.qleft = S.fx > S.fb
                    S.qhi = True
                    return
            else:
                S.status = -1
                S.qleft = S.fx < S.fb
                S.qhi = False
                return

        elif S.next_state == 80:
            if (abs(S.fc) < abs(S.fb)):
                if S.c != S.a:
                    S.d = S.a
                    S.fd = S.fa
                else:
                    S.a = S.b
                    S.fa = S.fb
                    S.xlo = S.c
                    S.b = S.xlo
                    S.fb = S.fc
                    S.c = S.a
                    S.fc = S.fa
            S.tol = 0.5 * max(S.atol, S.rtol * abs(S.xlo))
            S.m = (S.c + S.b) * 0.5
            S.mb = S.m - S.b

            if (not (abs(S.mb) > S.tol)):
                S.next_state = 240
                return
            else:
                if (not (S.ext > 3)):
                    S.tol = -S.tol if S.mb < 0 else S.tol
                    S.p = (S.b - S.a)*S.fb
                    if S.first:
                        S.q = S.fa - S.fb
                        S.first = False
                    else:
            # 120
                        S.fdb = (S.fd - S.fb)/(S.d - S.b)
                        S.fda = (S.fd - S.fa)/(S.d - S.a)
                        S.p = S.fda*S.p
                        S.q = S.fdb*S.fa - S.fda*S.fb
            # 130
                    if (S.p < 0.):
                        S.p = -S.p
                        S.q = -S.q
            # 140
                    if S.ext == 3:
                        S.p *= 2.

                    if not ((S.p*1. == 0.) or (S.p <= S.q*S.tol)):
            # 150
                        if not (S.p < (S.mb*S.q)):
            # 160
                            S.w = S.mb
                        else:
                            S.w = S.p / S.q
                    else:
                        S.w = S.tol
                else:
                    S.w = S.mb
            # 190
                S.d = S.a
                S.fd = S.fa
                S.a = S.b
                S.fa = S.fb
                S.b = S.b + S.w
                S.xlo = S.b
                S.x = S.xlo
                S.next_state = 200
                S.status = 1
                return

        elif S.next_state == 200:

            S.fb = S.fx
            if not (S.fc*S.fb >= 0.):
                if not (S.w == S.mb):
                    S.ext += 1
                else:
                    S.ext = 0
            else:
                S.c = S.a
                S.fc = S.fa
                S.ext = 0
            S.next_state = 80

        elif S.next_state == 240:
            S.xhi = S.c
            S.qrzero = (((S.fc >= 0.) and (S.fb <= 0.)) or
                        ((S.fc < 0.) and (S.fb >= 0.)))

            if S.qrzero:
                S.status = 0
            else:
                S.status = -1

            return

        else:
            S.status = -9999  # Bad state, should not be possible to get here
            return

# %%-------------------------------------- erf
cdef inline double erf(double x) noexcept nogil:
    cdef double ax, bot, t, top
    cdef double c = .564189583547756
    cdef double[5] a = [.771058495001320e-04, -.133733772997339e-02,
                        .323076579225834e-01, .479137145607681e-01,
                        .128379167095513e+00]
    cdef double[3] b = [.301048631703895e-02, .538971687740286e-01,
                        .375795757275549e+00]
    cdef double[8] p = [-1.36864857382717e-07, 5.64195517478974e-01,
                        7.21175825088309e+00, 4.31622272220567e+01,
                        1.52989285046940e+02, 3.39320816734344e+02,
                        4.51918953711873e+02, 3.00459261020162e+02]
    cdef double[8] q = [1.00000000000000e+00, 1.27827273196294e+01,
                        7.70001529352295e+01, 2.77585444743988e+02,
                        6.38980264465631e+02, 9.31354094850610e+02,
                        7.90950925327898e+02, 3.00459260956983e+02]
    cdef double[5] r = [2.10144126479064e+00, 2.62370141675169e+01,
                        2.13688200555087e+01, 4.65807828718470e+00,
                        2.82094791773523e-01]
    cdef double[4] s = [9.41537750555460e+01, 1.87114811799590e+02,
                        9.90191814623914e+01, 1.80124575948747e+01]

    ax = abs(x)
    if ax <= 0.5:
        t = x*x
        top = ((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]) + 1.
        bot = ((b[0]*t+b[1])*t+b[2])*t + 1.
        return x*(top/bot)

    if ax <= 4.:
        top = (((((((p[0]
                    )*ax+p[1]
                   )*ax+p[2]
                  )*ax+p[3]
                 )*ax+p[4]
                )*ax+p[5]
               )*ax+p[6]
              )*ax + p[7]
        bot = (((((((q[0]
                    )*ax+q[1]
                   )*ax+q[2]
                  )*ax+q[3]
                 )*ax+q[4]
                )*ax+q[5]
               )*ax+q[6])*ax + q[7]
        t = 0.5 + (0.5 - exp(-x*x)*(top/bot))
        return -t if (x < 0) else t

    if ax < 5.8:
        t = (1/x)**2
        top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t + r[4]
        bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t + 1.
        t = 0.5 + (0.5 - exp(-x*x) * (c - top/(x*x*bot))/ax)
        return -t if (x < 0) else t

    return -1. if (x < 0) else 1.

# %%-------------------------------------- erfc1
cdef inline double erfc1(int ind, double x) noexcept nogil:
    cdef double ax, bot, t, top, result
    cdef double c = 0.564189583547756
    cdef double[5] a = [.771058495001320e-04, -.133733772997339e-02,
                        .323076579225834e-01, .479137145607681e-01,
                        .128379167095513e+00]
    cdef double[3] b = [.301048631703895e-02, .538971687740286e-01,
                        .375795757275549e+00]
    cdef double[8] p = [-1.36864857382717e-07, 5.64195517478974e-01,
                        7.21175825088309e+00, 4.31622272220567e+01,
                        1.52989285046940e+02, 3.39320816734344e+02,
                        4.51918953711873e+02, 3.00459261020162e+02]
    cdef double[8] q = [1.00000000000000e+00, 1.27827273196294e+01,
                        7.70001529352295e+01, 2.77585444743988e+02,
                        6.38980264465631e+02, 9.31354094850610e+02,
                        7.90950925327898e+02, 3.00459260956983e+02]
    cdef double[5] r = [2.10144126479064e+00, 2.62370141675169e+01,
                        2.13688200555087e+01, 4.65807828718470e+00,
                        2.82094791773523e-01]
    cdef double[4] s = [9.41537750555460e+01, 1.87114811799590e+02,
                        9.90191814623914e+01, 1.80124575948747e+01]

    if x <= -5.6:
        return (2*exp(x*x)) if (ind != 0) else 2.

    # sqrt(log(np.finfo(np.float64).max)) ~= 26.64
    if (ind == 0) and (x > 26.64):
        return 0.

    ax = abs(x)

    if ax <= 0.5:
        t = x*x
        top = (((((a[0])*t+a[1])*t+a[2])*t+a[3])*t+a[4]) + 1.
        bot = (((b[0])*t+b[1])*t+b[2])*t + 1.
        result = 0.5 + (0.5 - x*(top/bot))
        return result if ind == 0 else result*exp(t)

    if 0.5 < ax <= 4.:
        top = (((((((p[0]
                    )*ax+p[1]
                   )*ax+p[2]
                  )*ax+p[3]
                 )*ax+p[4]
                )*ax+p[5]
               )*ax+p[6]
              )*ax + p[7]
        bot = (((((((q[0]
                  )*ax+q[1]
                 )*ax+q[2]
                )*ax+q[3]
               )*ax+q[4]
              )*ax+q[5]
             )*ax+q[6])*ax + q[7]
        result = top / bot
    else:
        t = (1 / x)**2
        top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t + r[4]
        bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t + 1.
        result = (c - t*(top/bot)) / ax

    if ind == 0:
        result *= exp(-(x*x))
        return (2. - result) if (x < 0) else result
    else:
        return (2.*exp(x*x) - result) if (x < 0) else result

# %%----------------------------------------- esum
cdef inline double esum(int mu, double x) noexcept nogil:
    if x > 0.:
        if (mu > 0.) or (mu + x < 0):
            return exp(mu)*exp(x)
        else:
            return exp(mu + x)
    else:
        if (mu < 0.) or (mu + x > 0.):
            return exp(mu)*exp(x)
        else:
            return exp(mu + x)

# %%----------------------------------------- fpser
cdef inline double fpser(double a, double b, double x, double eps) noexcept nogil:
    cdef double an, c, s, t, tol
    cdef double result = 1.

    result = 1.
    if (a > 1e-3*eps):
        result = 0.
        t = a*log(x)
        if t < -708.:
            return result
        result = exp(t)

    result *= (b / a)
    tol, an, t = eps / a, a + 1., x
    s = t / an
    while True:
        an += 1
        t *= x
        c = t / an
        s += c
        if abs(c) <= tol:
            break

    return result*(1. + a*s)

# %%----------------------------------------- gam1
cdef inline double gam1(double a) noexcept nogil:
    cdef double bot, d, t, top, w
    cdef double[7] p = [.577215664901533e+00, -.409078193005776e+00,
                        -.230975380857675e+00, .597275330452234e-01,
                        .766968181649490e-02, -.514889771323592e-02,
                        .589597428611429e-03]
    cdef double[5] q = [.100000000000000e+01, .427569613095214e+00,
                        .158451672430138e+00, .261132021441447e-01,
                        .423244297896961e-02]
    cdef double[9] r = [-.422784335098468e+00, -.771330383816272e+00,
                        -.244757765222226e+00, .118378989872749e+00,
                        .930357293360349e-03, -.118290993445146e-01,
                        .223047661158249e-02, .266505979058923e-03,
                        -.132674909766242e-03]
    cdef double[2] s = [.273076135303957e+00, .559398236957378e-01]

    d = a - 0.5
    t = d - 0.5 if (d > 0) else a

    if t == 0.:
        return 0.

    if t < 0:
        top = ((((((((r[8]
                     )*t+r[7]
                    )*t+r[6]
                   )*t+r[5]
                  )*t+r[4]
                 )*t+r[3]
                )*t+r[2]
               )*t+r[1]
              )*t + r[0]
        bot = (s[1]*t + s[0])*t + 1.
        w = top / bot
        if d > 0.:
            return t*w/a
        else:
            return a * ((w + 0.5) + 0.5)

    top = ((((((p[6]
               )*t+p[5]
              )*t+p[4]
             )*t+p[3]
            )*t+p[2]
           )*t+p[1]
          )*t + p[0]
    bot = ((((q[4]
             )*t+q[3]
            )*t+q[2]
           )*t+q[1]
          )*t + 1.
    w = top / bot
    if d > 0.:
        return (t/a) * ((w - 0.5) - 0.5)
    else:
        return a * w

# %%----------------------------------------- gaminv
cdef inline (double, int) gaminv(double a, double p,
                                 double q, double x0) noexcept nogil:
    cdef:
        double am1, ap1, ap2, ap3, apn, b, bot, d, g, h
        double pn, qg, qn, r, rta, s, s2, ssum, t, top, u, w, y, z
        double ln10 = log(10)
        double c = 0.57721566490153286060651209008
        double tol = 1e-5
        double e = spmpar[0]
        double e2 = 2*e
        double amax = 0.4e-10 / (e*e)
        double xmin = spmpar[1]
        double xmax = spmpar[2]
        double xn = x0
        double x = 0.
        int ierr = 0
        int iop = 1 if e > 1e-10 else 0
        bint use_p

    cdef double[4] arr = [3.31125922108741, 11.6616720288968,
                          4.28342155967104, 0.213623493715853]
    cdef double[4] barr = [6.61053765625462, 6.40691597760039,
                           1.27364489782223, 0.036117081018842]
    cdef double[2] eps0 = [1e-10, 1e-8]
    cdef double[2] amin = [500., 100.]
    cdef double[2] bmin = [1e-28, 1e-13]
    cdef double[2] dmin = [1e-6, 1e-4]
    cdef double[2] emin = [2e-3, 6e-3]
    cdef double eps = eps0[iop]

    if a <= 0.:
        return (x, -2)
    t = p + q - 1.
    if abs(t) > e:
        return (x, -4)

    ierr = 0
    if p == 0.:
        return (x, 0)
    if q == 0.:
        return (xmax, 0)
    if a == 1.:
        return (-log(q), 0)

    if x0 > 0.:
        use_p = False if p > 0.5 else True
        am1 = (a - 0.5) - 0.5
        if (p if use_p else q) <= 1.e10*xmin:
            return (xn, -8)

    elif a <= 1.:
        g = gamma(a + 1.)
        qg = q*g
        if qg == 0.:
            return (xmax, -8)
        b = qg / a
        if (qg > 0.6*a) or (((a >= 0.3) or (b < 0.35)) and (b > 0.45)):

            if (b*q > 1.e-8):
                if (p <= 0.9):
                    xn = exp(log(p*g)/a)
                else:
                    xn = exp((alnrel(-q) + gamln1(a))/a)
            else:
                xn = exp(-(q/a + c))

            if xn == 0.:
                return (x, -3)

            t = 0.5 + (0.5 - xn/(a + 1.))
            xn /= t
            use_p = True
            am1 = (a - 0.5) - 0.5
            if (p if use_p else q) <= 1.e10*xmin:
                return (xn, -8)

        elif not ((a >= 0.3) or (b < 0.35)):
            t = exp(-(b+c))
            u = t*exp(t)
            xn = t*exp(u)
            use_p = True
            am1 = (a - 0.5) - 0.5
            if p <= 1.e10*xmin:
                return (xn, -8)
        else:
            if b == 0.:
                return (xmax, -8)

            y = -log(b)
            s = 0.5 + (0.5 - a)
            z = log(y)
            t = y - s*z
            if b <= 0.01:
                xn = gaminv_helper_30(a, s, y, z)
                if (a <= 1.) or (b > bmin[iop]):
                    return (xn, 0)
            elif b < 0.15:
                u = ((t+2.*(3.-a))*t + (2.-a)* (3.-a))/((t+ (5.-a))*t+2.)
                xn = y - s*log(t) - log(u)
            else:
                xn = y - s*log(t)-log(1.+s/(t+1.))
            use_p = False
            am1 = (a - 0.5) - 0.5
            if q <= 1.e10*xmin:
                return (xn, -8)

    elif a > 1:
        w = log(p) if q > 0.5 else log(q)
        t = sqrt(-2.*w)
        top = ((arr[3]*t+arr[2])*t+arr[1])*t+arr[0]
        bot = (((barr[3]*t+barr[2])*t+barr[1])*t+barr[0])*t+1.
        s = t - top/bot
        if q >= 0.5:
            s = -s

        rta = sqrt(a)
        s2 = s*s
        xn = a + s*rta
        xn += (s2 - 1.)/3. + s*(s2-7.)/(36.*rta)
        xn -= ((3.*s2+7.)*s2-16.)/(810.*a)
        xn += s*((9.*s2+256.)*s2-433.)/(38880.*a*rta)
        xn = max(xn, 0.)
        if a >= amin[iop]:
            d = 0.5 + (0.5 - xn/a)
            if abs(d) <= dmin[iop]:
                return (xn, 0)

        if p > 0.5:
            # All branches that go to 220
            if xn < 3.*a:
                pass
            else:
                y = - (w + gamln(a))
                d = max(2., a*(a - 1.))
                if y < ln10*d:
                    t = a - 1.
                    # Note: recursion not duplicate
                    xn = y + t*log(xn) - alnrel(-t/(xn+1.))
                    xn = y + t*log(xn) - alnrel(-t/(xn+1.))
                else:
                    s = 1. - a
                    z = log(y)
                    xn = gaminv_helper_30(a, s, y, z)
                    if (a <= 1.):
                        return (xn, 0)

            use_p = False
            am1 = (a - 0.5) - 0.5
            if q <= 1.e10*xmin:
                return (xn, -8)

        else:
            # all branches that go to 170
            ap1 = a + 1.
            if xn > 0.7*ap1:
                pass
            else:
                # branches that go to 170 via 140
                w += gamln(ap1)
                ap2 = a + 2.
                ap3 = a + 3.
                # Note: recursion not duplicate
                x = exp((w+x)*a)
                x = exp((w+x-log(1.+(x/ap1)*(1.+x/ap2)))/a)
                x = exp((w+x-log(1.+(x/ap1)*(1.+x/ap2)))/a)
                x = exp((w+x-log(1.+(x/ap1)*(1.+(x/ap2)*(1.+x/ap3))))/a)

                # push xn = x inside the following if
                if (xn > 0.15*ap1) or (x > 0.01*ap1):
                    if not (xn > 0.15*ap1):
                        xn = x

                    # Go to 140
                    apn = ap1
                    t = xn/apn
                    ssum = 1. + t
                    while True:
                        t *= xn/apn
                        ssum += t
                        if t <= 1e-4:
                            break
                    t = w - log(ssum)
                    xn = exp((xn + t)/a)
                    xn *= (1. - (a*log(xn) - xn -t) / (a - xn))
                    # Go to 170
                elif x <= emin[iop]*ap1:
                    return (x, 0)

            # Go to 170
            use_p = True
            am1 = (a - 0.5) - 0.5
            if p <= 1.e10*xmin:
                return (xn, -8)

    # Schroder iteration using P or Q
    for ierr in range(20):
        if a > amax:
            d = 0.5 + (0.5 - xn / a)
            if abs(d) <= e2:
                return (xn, -8)

        pn, qn = gratio(a, xn, 0)
        if pn == 0. or qn == 0.:
            return (xn, -8)
        r = rcomp(a, xn)
        if r == 0.:
            return (xn, -8)
        t = (pn-p)/r if use_p else (q-qn)/r
        w = 0.5 * (am1 - xn)
        if (abs(t) <= 0.1) and (abs(w*t) <= 0.1):
            h = t * (1. + w*t)
            x = xn * (1. - h)
            if x <= 0.:
                return (x, -7)
            if (abs(w) >= 1.) and (abs(w)*t*t <= eps):
                return (x, ierr)
            d = abs(h)
        else:
            x = xn * (1. - t)
            if x <= 0.:
                return (x, -7)
            d = abs(t)

        xn = x

        if d <= tol:
            act_val = abs(pn-p) if use_p else abs(qn-q)
            tol_val = tol*p if use_p else tol*q
            if (d <= eps) or (act_val <= tol_val):
                return (x, ierr)
    else:
        return (x, -6)


cdef inline double gaminv_helper_30(double a, double s,
                                    double z, double y) noexcept nogil:
    cdef double c1, c2, c3, c4, c5
    c1 = -s*z
    c2 = -s*(1. + c1)
    c3 = s*((0.5*c1+ (2.-a))*c1 + (2.5-1.5*a))
    c4 = -s*(((c1/3. + (2.5-1.5*a))*c1 + ((a-6.)*a+7.))*c1 + ((11.*a-46.)*a+47.)/6.)
    c5 = -s*((((-c1/4.+ (11.*a-17.)/6.
               )*c1+ ((-3.*a+13.)*a-13.)
              )*c1 + 0.5*(((2.*a-25.)*a+72.)*a-61.)
             )*c1+ (((25.*a-195.)*a+477.)*a-379.)/12.)
    return ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y

# %%----------------------------------------- gamln
cdef inline double gamln(double a) noexcept nogil:
    cdef double t,w
    cdef double c[6]
    cdef double d = .418938533204673
    cdef int i,n

    c[0] = .833333333333333e-01
    c[1] = -.277777777760991e-02
    c[2] = .793650666825390e-03
    c[3] = -.595202931351870e-03
    c[4] = .837308034031215e-03
    c[5] = -.165322962780713e-02

    if a <= 0.8:
        return gamln1(a) - log(a)

    if a <= 2.25:
        t = (a-0.5) - 0.5
        return gamln1(t)

    if a < 10:
        n = <int>(a - 1.25)
        t = a
        w = 1.
        for i in range(n):
            t -= 1.
            w *= t
        return gamln1(t-1.) + log(w)

    t = (1/a)**2
    w = (((((c[5]*t+c[4])*t+c[3])*t+c[2])*t+c[1])*t+c[0])/a
    return (d + w) + (a-0.5)*(log(a) - 1.)

# %%----------------------------------------- gamln1
cdef inline double gamln1(double a) noexcept nogil:
    cdef double p[7]
    cdef double q[6]
    cdef double r[6]
    cdef double s[5]
    cdef double bot, top, w, x

    p[0] = .577215664901533e+00
    p[1] = .844203922187225e+00
    p[2] = -.168860593646662e+00
    p[3] = -.780427615533591e+00
    p[4] = -.402055799310489e+00
    p[5] = -.673562214325671e-01
    p[6] = -.271935708322958e-02
    q[1] = .288743195473681e+01
    q[2] = .312755088914843e+01
    q[3] = .156875193295039e+01
    q[4] = .361951990101499e+00
    q[5] = .325038868253937e-01
    q[6] = .667465618796164e-03
    r[0] = .422784335098467e+00
    r[1] = .848044614534529e+00
    r[2] = .565221050691933e+00
    r[3] = .156513060486551e+00
    r[4] = .170502484022650e-01
    r[5] = .497958207639485e-03
    s[1] = .124313399877507e+01
    s[2] = .548042109832463e+00
    s[3] = .101552187439830e+00
    s[4] = .713309612391000e-02
    s[5] = .116165475989616e-03

    if a < 0.6:
        top = ((((((p[6]
                   )*a+p[5]
                  )*a+p[4]
                 )*a+p[3]
                )*a+p[2]
               )*a+p[1]
              )*a+p[0]
        bot = ((((((q[6]
                   )*a+q[5]
                  )*a+q[4]
                 )*a+q[3]
                )*a+q[2]
               )*a+q[1]
              )*a+1.
        w = top/bot
        return -a*w
    else:
        x = (a - 0.5) - 0.5
        top = (((((r[5]
                  )*x+r[4]
                 )*x+r[3]
                )*x+r[2]
               )*x+r[1]
              )*x+r[0]
        bot = (((((s[5]
                  )*x+s[4]
                 )*x+s[3]
                )*x+s[2]
               )*x+s[1]
              )*x+1.
        w = top/bot
        return x*w

# %%-------------------------------------- gamma
cdef inline double gamma(double a) noexcept nogil:
    cdef double bot, g, lnx, s, t, top, w, z, result
    cdef int i, j, m, n

    cdef double d = 0.5*(log(2.*PI) - 1)
    cdef double x = a
    cdef double[7] p = [.539637273585445e-03, .261939260042690e-02,
                        .204493667594920e-01, .730981088720487e-01,
                        .279648642639792e+00, .553413866010467e+00,
                        1.0]
    cdef double[7] q = [-.832979206704073e-03, .470059485860584e-02,
                        .225211131035340e-01, -.170458969313360e+00,
                        -.567902761974940e-01, .113062953091122e+01,
                        1.0]

    cdef double[5] r = [.820756370353826e-03, -.595156336428591e-03,
                        .793650663183693e-03, -.277777777770481e-02,
                        .833333333333333e-01]

    result = 0.
    if abs(a) < 15:
        t = 1.
        m = <int>a - 1
        if m > 0:
            for j in range(m):
                x -= 1.
                t *= x
            x -= 1.
        elif m == 0:
            x -= 1.
        else:
            t = a
            if a <= 0.:
                m = -m - 1
                if m != 0.:
                    for j in range(m):
                        x += 1.
                        t *= x
                x += 0.5
                x += 0.5
                t *= x

                if t == 0.:
                    return result
            if abs(t) < 1e-30:
                if abs(t)*spmpar[2] <= 1.0001:
                    return result
                return 1./t

        top = p[0]
        bot = q[0]
        for i in range(1, 7):
            top *= x
            top += p[i]
            bot *= x
            bot += q[i]
        result = top / bot

        return result/t if (a < 1.) else result*t

    if abs(a) >= 1.e3:
        return result

    if a <= 0.:
        x = -a
        n = <int>x
        t = x - n
        if t > 0.9:
            t = 1. - t
        s = sin(PI*t) / PI
        if n % 2 == 0:
            s = -s
        if s == 0.:
            return result

    t = (1 / x)**2
    g = ((((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4]) / x
    lnx = log(x)
    z = x
    g = (d + g) + (z -0.5)*(lnx - 1.)
    w = g
    t = g - w
    if (w > 0.99999*709):
        return result
    result = exp(w)*(1. + t)
    return (1. / (result * s)) / x if (a < 0.) else result

# %%-------------------------------------- grat1
# Subroutine converted to function
cdef inline (double, double) grat1(double a, double x,
                                   double r, double eps) noexcept nogil:
    cdef double a2n, a2nm1, am0, an, an0, b2n, b2nm1, c, cma, g, h, j, l
    cdef double p, q, ssum, t, tol, w, z

    if a*x == 0.:
        return (1., 0.) if (x > a) else (0., 1.)

    if a == 0.5:
        if x < 0.25:
            p = erf(sqrt(x))
            return (p, 0.5 + (0.5 - p))
        else:
            q = erfc1(0, sqrt(x))
            return (0.5 + (0.5 - q), q)

    if x < 1.1:
        an = 3.
        c = x
        ssum = x / (a + 3.)
        tol = 0.1*eps / (a + 1.)
        while True:
            an += 1
            c *= -(x / an)
            t = c / (a + an)
            ssum += t
            if abs(t) < tol:
                break
        j = a*x*((ssum/6. - 0.5/(a+2.))*x + 1./(a+1.))

        z = a * log(x)
        h = gam1(a)
        g = 1. + h

        if ((x < 0.25) and (z > -.13394)) or (a < x / 2.59):
            l = rexp(z)
            w = 0.5 + (0.5 - l)
            q = (w*j - l)*g - h
            if q < 0.:
                return (1., 0.)
            p = 0.5 + (0.5 - q)
            return (p, q)
        else:
            w = exp(z)
            p = w*g*(0.5 + (0.5 - j))
            q = 0.5 + (0.5 - p)
            return (p, q)

    a2nm1 = 1.
    a2n = 1.
    b2nm1 = x
    b2n = x + (1. - a)
    c = 1.
    while True:
        a2nm1 = x*a2n + c*a2nm1
        b2nm1 = x*b2n + c*b2nm1
        am0 = a2nm1/b2nm1
        c = c + 1.
        cma = c - a
        a2n = a2nm1 + cma*a2n
        b2n = b2nm1 + cma*b2n
        an0 = a2n/b2n
        if (abs(an0-am0) < eps*an0):
            break
    q = r*an0
    p = 0.5 + (0.5 - q)
    return (p, q)

# %%-------------------------------------- gratio
cdef inline (double, double) gratio(double a, double x, int ind) noexcept nogil:
    cdef:
        double d10 = -.185185185185185e-02
        double d20 = .413359788359788e-02
        double d30 = .649434156378601e-03
        double d40 = -.861888290916712e-03
        double d50 = -.336798553366358e-03
        double d60 = .531307936463992e-03
        double d70 = .344367606892378e-03
        double alog10 = log(10)
        double rt2pi = sqrt(1. / (2.*PI))
        double rtpi = sqrt(PI)
        double eps = spmpar[0]
        double a2n, a2nm1, am0, amn, an, an0, apn, b2n, b2nm1
        double c, c0, c1, c2, c3, c4, c5, c6, cma, e0, g, h, j, l, r
        double rta, rtx, s, ssum, t, t1, tol, twoa, u, w, x0, y, z
        int i, m, n

    cdef double[20] wk
    cdef double[3] acc0 = [5.e-15, 5.e-7, 5.e-4]
    cdef double[3] big = [20., 14., 10.]
    cdef double[3] e00 = [.00025, .025, .14]
    cdef double[3] x00 = [31., 17., 9.7]
    cdef double[13] d0 = [.833333333333333e-01, -.148148148148148e-01,
                          .115740740740741e-02, .352733686067019e-03,
                          -.178755144032922e-03, .391926317852244e-04,
                          -.218544851067999e-05, -.185406221071516e-05,
                          .829671134095309e-06, -.176659527368261e-06,
                          .670785354340150e-08, .102618097842403e-07,
                          -.438203601845335e-08]
    cdef double[12] d1 = [-.347222222222222e-02, .264550264550265e-02,
                          -.990226337448560e-03, .205761316872428e-03,
                          -.401877572016461e-06, -.180985503344900e-04,
                          .764916091608111e-05, -.161209008945634e-05,
                          .464712780280743e-08, .137863344691572e-06,
                          -.575254560351770e-07, .119516285997781e-07]
    cdef double[10] d2 = [-.268132716049383e-02, .771604938271605e-03,
                          .200938786008230e-05, -.107366532263652e-03,
                          .529234488291201e-04, -.127606351886187e-04,
                          .342357873409614e-07, .137219573090629e-05,
                          -.629899213838006e-06, .142806142060642e-06]
    cdef double[8] d3 = [.229472093621399e-03, -.469189494395256e-03,
                         .267720632062839e-03, -.756180167188398e-04,
                         -.239650511386730e-06, .110826541153473e-04,
                         -.567495282699160e-05, .142309007324359e-05]
    cdef double[6] d4 = [.784039221720067e-03, -.299072480303190e-03,
                         -.146384525788434e-05, .664149821546512e-04,
                         -.396836504717943e-04, .113757269706784e-04]
    cdef double[4] d5 = [-.697281375836586e-04, .277275324495939e-03,
                         -.199325705161888e-03, .679778047793721e-04]
    cdef double[2] d6 = [-.592166437353694e-03, .270878209671804e-03]

    if a < 0. or x < 0.:
        return 2., 0.
    if a == 0. and x == 0.:
        return 2., 0.

    if a*x == 0.:
        return (1., 0.) if x > a else (0., 1.)

    if ind != 0 and ind != 1:
        ind = 2

    acc = max(acc0[ind], eps)
    e0, x0 = e00[ind], x00[ind]

    if a >= 1.:
        if a >= big[ind]:
            l = x / a
            if l == 0.:
                return (0., 1)
            s = 0.5 + (0.5 - l)
            z = rlog(l)
            if z >= (700. / a):
                if abs(s) <= 2.*eps:
                    return (2., 0)
                return (1., 0.) if x > a else (0., 1.)
            y = a*z
            rta = sqrt(a)
            if abs(s) <= (e0 / rta):
                if (a*eps*eps > 3.28e-3):
                    return (2., 0.)
                c = 0.5 + (0.5 - y)
                w = (0.5 - sqrt(y) * (0.5 + (0.5 - y/3.))/rtpi)/c
                u = 1. / a
                z = sqrt(z+z)
                if l < 1.:
                    z = -z
                if ind == 0:
                    c0 = (((((((d0[6]
                               )*z+d0[5]
                              )*z+d0[4]
                             )*z+d0[3]
                            )*z+d0[2]
                           )*z+d0[1]
                          )*z+d0[0]
                         )*z - (1./3.)
                    c1 = ((((((d1[5]
                              )*z+d1[4]
                             )*z+d1[3]
                            )*z+d1[2]
                           )*z+d1[1]
                          )*z+d1[0]
                         )*z + d10
                    c2 = (((((d2[4])*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z + d20
                    c3 = ((((d3[3])*z+d3[2])*z+d3[1])*z+d3[0])*z + d30
                    c4 = (d4[1]*z+d4[0])*z + d40
                    c5 = (d5[1]*z+d5[0])*z + d50
                    c6 = d6[0]*z + d60
                    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0

                elif ind == 1:
                    c0 = (d0[1]*z+d0[0])*z - (1. / 3.)
                    c1 = d1[0]*z + d10
                    t = (d20*u+c1)*u + c0

                else:
                    t = d0[0]*z - (1. / 3.)

                if l < 1.:
                    ans = c * (w - rt2pi*t/rta)
                    qans = 0.5 + (0.5 - ans)
                else:
                    qans = c * (w + rt2pi*t/rta)
                    ans = 0.5 + (0.5 - qans)

                return ans, qans
            if abs(s) <= 0.4:
                if (abs(s) <= 2.*eps) and (a*eps*eps > 3.28e-3):
                    return (2., 0.)
                c = exp(-y)
                w = 0.5*erfc1(1, sqrt(y))
                u = 1./a
                z = sqrt(z+z)
                if l < 1.:
                    z = -z
                if ind == 0:
                    if abs(s) <= 1e-3:
                        c0 = (((((((d0[6]
                                   )*z+d0[5]
                                  )*z+d0[4]
                                 )*z+d0[3]
                                )*z+d0[2]
                               )*z+d0[1]
                              )*z+d0[0]
                             )*z - (1./3.)
                        c1 = ((((((d1[5]
                                  )*z+d1[4]
                                 )*z+d1[3]
                                )*z+d1[2]
                               )*z+d1[1]
                              )*z+d1[0]
                             )*z + d10
                        c2 = (((((d2[4])*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z + d20
                        c3 = ((((d3[3])*z+d3[2])*z+d3[1])*z+d3[0])*z + d30
                        c4 = (d4[1]*z+d4[0])*z + d40
                        c5 = (d5[1]*z+d5[0])*z + d50
                        c6 = d6[0]*z + d60
                        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
                    else:
                        c0 = (((((((((((((d0[12]
                                         )*z+d0[11]
                                        )*z+d0[10]
                                       )*z+d0[9]
                                      )*z+d0[8]
                                     )*z+d0[7]
                                    )*z+d0[6]
                                   )*z+d0[5]
                                  )*z+d0[4]
                                 )*z+d0[3]
                                )*z+d0[2]
                               )*z+d0[1]
                              )*z+d0[0]
                             )*z - (1./3.)
                        c1 = ((((((((((((d1[11]
                                        )*z+d1[10]
                                       )*z+d1[9]
                                      )*z+d1[8]
                                     )*z+d1[7]
                                    )*z+d1[6]
                                   )*z+d1[5]
                                  )*z+d1[4]
                                 )*z+d1[3]
                                )*z+d1[2]
                               )*z+d1[1]
                              )*z+d1[0]
                             )*z + d10
                        c2 = ((((((((((d2[9]
                                      )*z+d2[8]
                                     )*z+d2[7]
                                    )*z+d2[6]
                                   )*z+d2[5]
                                  )*z+d2[4]
                                 )*z+d2[3]
                                )*z+d2[2]
                               )*z+d2[1]
                              )*z+d2[0]
                             )*z + d20
                        c3 = ((((((((d3[7]
                                    )*z+d3[6]
                                   )*z+d3[5]
                                  )*z+d3[4]
                                 )*z+d3[3]
                                )*z+d3[2]
                               )*z+d3[1]
                              )*z+d3[0]
                             )*z + d30
                        c4 = ((((((d4[5]
                                  )*z+d4[4]
                                 )*z+d4[3]
                                )*z+d4[2]
                               )*z+d4[1]
                              )*z+d4[0]
                             )*z + d40
                        c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z + d50
                        c6 = (d6[1]*z+d6[0])*z + d60
                        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
                elif ind == 1:
                    c0 = (((((d0[5]*z+d0[4])*z+d0[3])*z+d0[2])*z+d0[1])*z+d0[0])*z - (1./3.)
                    c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z + d10
                    c2 = d2[0]*z + d20
                    t = (c2*u+c1)*u + c0
                else:
                    t = ((d0[2]*z+d0[1])*z+d0[0])*z - (1./3.)

                if l < 1.:
                    ans = c * (w - rt2pi*t/rta)
                    qans = 0.5 + (0.5 - ans)
                else:
                    qans = c * (w + rt2pi*t/rta)
                    ans = 0.5 + (0.5 - qans)

                return ans, qans

            t = (1 / a)**2
            t1 = (((0.75*t-1.)*t+3.5)*t-105.0)/ (a*1260.0)
            t1 -= y
            r = rt2pi*rta*exp(t1)
            # 40
            if r == 0.:
                return (1., 0.) if x > a else (0., 1.)
            # 50
            if x <= max(a, alog10):
                apn = a + 1.
                t = x / apn
                wk[0] = t

                for n in range(1, 20):
                    apn += 1.
                    t *= (x / apn)
                    if t <= 1e-3:
                        break
                    wk[n] = t

                ssum = t
                tol = 0.5 * acc
                while True:
                    apn += 1.
                    t *= x / apn
                    ssum += t
                    if t <= tol:
                        break

                for m in range(n - 1, -1, -1):
                    ssum += wk[m]
                ans = (r/a) * (1. + ssum)
                qans = 0.5 + (0.5 - ans)
                return ans, qans

            # 250
            if x <= max(a, alog10):
                apn = a + 1.
                t = x / apn
                wk[0] = t

                for n in range(1, 20):
                    apn += 1.
                    t *= (x / apn)
                    if t <= 1e-3:
                        break
                    wk[n] = t

                ssum = t
                tol = 0.5 * acc
                while True:
                    apn += 1.
                    t *= x / apn
                    ssum += t
                    if t <= tol:
                        break

                for m in range(n - 1, -1, -1):
                    ssum += wk[m]
                ans = (r/a) * (1. + ssum)
                qans = 0.5 + (0.5 - ans)
                return ans, qans

            # 100
            if x < x0:
                tol = max(5.*eps, acc)
                a2nm1, a2n, c = 1., 1., 1.
                b2nm1, b2n = x, x + (1. - a)

                while True:
                    a2nm1 = x*a2n + c*a2nm1
                    b2nm1 = x*b2n + c*b2nm1
                    am0 = a2nm1 / b2nm1
                    c += 1.
                    cma = c - a
                    a2n = a2nm1 + cma*a2n
                    b2n = b2nm1 + cma*b2n
                    an0 = a2n/b2n
                    if (abs(an0-am0) < tol*an0):
                        break
                qans = r*an0
                ans = 0.5 + (0.5 - qans)
                return ans, qans

            amn = a - 1.
            t = amn / x
            wk[0] = t
            for n in range(2, 20):
                amn -= 1.
                t *= amn / x
                if abs(t) <= 1e-3:
                    break

            ssum = t
            while not (abs(t) <= acc):
                amn -= 1.
                t *= amn / x
                ssum += t

            for m in range(n - 1, -1, -1):
                ssum += wk[n]

            qans = (r/x) * (1. + ssum)
            ans = 0.5 + (0.5 - qans)
            return ans, qans



        twoa = a + a
        m = <int>(twoa)

        if ((a > x) or (x >= x0)) or (twoa != m):
            t1 = a*log(x) - x
            r = exp(t1)/gamma(a)

            if r == 0.:
                return (1., 0.) if x > a else (0., 1.)

            if x <= max(a, alog10):
                apn = a + 1.
                t = x / apn
                wk[0] = t

                for n in range(1, 20):
                    apn += 1.
                    t *= (x / apn)
                    if t <= 1e-3:
                        break
                    wk[n] = t

                ssum = t
                tol = 0.5 * acc
                while True:
                    apn += 1.
                    t *= x / apn
                    ssum += t
                    if t <= tol:
                        break

                for m in range(n - 1, -1, -1):
                    ssum += wk[m]
                ans = (r/a) * (1. + ssum)
                qans = 0.5 + (0.5 - ans)
                return ans, qans

            if x < x0:
                tol = max(5.*eps, acc)
                a2nm1, a2n, c = 1., 1., 1.
                b2nm1, b2n = x, x + (1. - a)

                while True:
                    a2nm1 = x*a2n + c*a2nm1
                    b2nm1 = x*b2n + c*b2nm1
                    am0 = a2nm1 / b2nm1
                    c += 1.
                    cma = c - a
                    a2n = a2nm1 + cma*a2n
                    b2n = b2nm1 + cma*b2n
                    an0 = a2n/b2n
                    if (abs(an0-am0) < tol*an0):
                        break
                qans = r*an0
                ans = 0.5 + (0.5 - qans)
                return ans, qans

            amn = a - 1.
            t = amn / x
            wk[0] = t
            for n in range(2, 20):
                amn -= 1.
                t *= amn / x
                if abs(t) <= 1e-3:
                    break

            ssum = t
            while not (abs(t) <= acc):
                amn -= 1.
                t *= amn / x
                ssum += t

            for m in range(n - 1, -1, -1):
                ssum += wk[n]

            qans = (r/x) * (1. + ssum)
            ans = 0.5 + (0.5 - qans)
            return ans, qans

        i = m // 2
        if a == i:
            ssum = exp(-x)
            t = ssum
            n = 1
            c = 0.
        else:
            rtx = sqrt(x)
            ssum = erfc1(0, rtx)
            t = exp(-x) / (rtpi*rtx)
            n = 0
            c = -0.5

        while n < i:
            n += 1
            c += 1.
            t *= x/c
            ssum += t

        qans = ssum
        ans = 0.5 + (0.5 - qans)
        return ans, qans

    if a == 0.5:
        if x >= 0.5:
            qans = erfc1(0, sqrt(x))
            ans = 0.5 + (0.5 - qans)
        else:
            ans = erf(sqrt(x))
            qans = 0.5 + (0.5 - ans)
        return ans, qans

    if x < 1.1:
        an = 3.
        c = x
        ssum = x / (a + 3.)
        tol = 3.*acc / (a + 1.)
        while True:
            an += 1.
            c *= -(x / an)
            t = c / (a + an)
            ssum += t
            if abs(t) <= tol:
                break
        j = a*x*((ssum / 6. - 0.5 / (a + 2.))*x + 1./(a + 1.))
        z = a*log(x)
        h = gam1(a)
        g = 1. + h

        if ((x < 0.25) and (z > -0.13394)) or (a < x/2.59):
            l = rexp(z)
            w = 0.5 + (0.5 + l)
            qans = (w*j - l)*g - h
            if qans < 0.:
                return (1., 0.)
            ans = 0.5 + (0.5 - qans)
        else:
            w = exp(z)
            ans = w*g*(0.5 + (0.5 - j))
            qans = 0.5 + (0.5 - ans)
        return ans, qans

    t1 = a*log(x) - x
    u = a*exp(t1)
    if u == 0.:
        return (1., 0.)
    r = u * (1. + gam1(a))
    tol = max(5.*eps, acc)
    a2nm1, a2n, c = 1., 1., 1.
    b2nm1, b2n = x, x + (1. - a)

    while True:
        a2nm1 = x*a2n + c*a2nm1
        b2nm1 = x*b2n + c*b2nm1
        am0 = a2nm1 / b2nm1
        c += 1.
        cma = c - a
        a2n = a2nm1 + cma*a2n
        b2n = b2nm1 + cma*b2n
        an0 = a2n/b2n
        if (abs(an0-am0) < tol*an0):
            break
    qans = r*an0
    ans = 0.5 + (0.5 - qans)
    return ans, qans

# %%-------------------------------------- gsumln
cdef inline double gsumln(double a, double b) noexcept nogil:
    cdef double x

    x = a + b - 2
    if x <= 0.25:
        return gamln1(1. + x)

    if x <= 1.25:
        return gamln1(x) + alnrel(x)

    return gamln1(x - 1.) + log(x*(1. + x))


# %%-------------------------------------- psi_fort
cdef inline double psi(double xx) noexcept nogil:
    cdef double aug, den, dx0, sgn, upper, w, x, xmax1, xmx0, xsmall, z
    cdef double p1[7]
    cdef double q1[6]
    cdef double p2[4]
    cdef double q2[4]
    cdef int nq, i

    dx0 = 1.461632144968362341262659542325721325
    p1[0] = .895385022981970e-02
    p1[1] = .477762828042627e+01
    p1[2] = .142441585084029e+03
    p1[3] = .118645200713425e+04
    p1[4] = .363351846806499e+04
    p1[5] = .413810161269013e+04
    p1[6] = .130560269827897e+04
    q1[0] = .448452573429826e+02
    q1[1] = .520752771467162e+03
    q1[2] = .221000799247830e+04
    q1[3] = .364127349079381e+04
    q1[4] = .190831076596300e+04
    q1[5] = .691091682714533e-05
    p2[0] = -.212940445131011e+01
    p2[1] = -.701677227766759e+01
    p2[2] = -.448616543918019e+01
    p2[3] = -.648157123766197e+00
    q2[0] = .322703493791143e+02
    q2[1] = .892920700481861e+02
    q2[2] = .546117738103215e+02
    q2[3] = .777788548522962e+01

    xmax1 = 4503599627370496.0
    xsmall = 1e-9

    x = xx
    aug = 0.
    if x < 0.5:
        if abs(x) <= xsmall:
            if x == 0.:
                return 0.
            aug = -1./x
        else:
            w = -x
            sgn = PI / 4
            if w <= 0.:
                w = -w
                sgn = -sgn

            if w >= xmax1:
                return 0.

            w = w - <int>w
            nq = <int>(w*4.)
            w = 4.*(w - 0.25*nq)

            if nq % 2 == 1:
                w = 1. - w
            z = (PI / 4.)*w

            if (nq // 2) % 2 == 1:
                sgn = -sgn

            if ((nq + 1) // 2) % 2 == 1:
                aug = sgn * (tan(x)*4.)
            else:
                if z == 0.:
                    return 0.
                aug = sgn * (4./tan(x))

        x = 1 - x

    if x <= 3.:
        den = x
        upper = p1[0]*x

        for i in range(5):
            den = (den + q1[i])*x
            upper = (upper + p1[i+1])*x

        den = (upper + p1[6]) / (den + q1[5])
        xmx0 = x - dx0
        return (den * xmx0) + aug

    else:
        if x < xmax1:
            w = 1. / (x*x)
            den = w
            upper = p2[0]*w

            for i in range(3):
                den = (den + q2[i])*w
                upper = (upper + p2[i+1])*w

            aug += upper / (den + q2[3]) - 0.5/x

        return aug + log(x)

# %%-------------------------------------- rcomp
cdef inline double rcomp(double a, double x) noexcept nogil:
    cdef double t, t1, u
    cdef double r2pi = sqrt(1. / (2.*PI))

    if a < 20:
        t = a*log(x) - x
        return a*exp(t)*(1. + gam1(a)) if (a < 1) else exp(t) / gamma(a)
    else:
        u = x / a
        if u == 0.:
            return 0.
        t = (1 / a)**2
        t1 = (((0.75*t-1.)*t+3.5)*t - 105.) / (a*1260)
        t1 -= a*rlog(u)
        return r2pi*sqrt(a)*exp(t1)

# %%-------------------------------------- rexp
cdef inline double rexp(double x) noexcept nogil:
    cdef double[2] p = [.914041914819518e-09, .238082361044469e-01]
    cdef double[4] q = [-.499999999085958e+00, .107141568980644e+00,
                        -.119041179760821e-01, .595130811860248e-03]
    cdef double w

    if abs(x) <= 0.15:
        return x*(((p[1]*x+p[0])*x+1.)/((((q[3]*x+q[2])*x+q[1])*x+q[0])*x + 1.))
    else:
        w = exp(x)
        if x > 0.:
            return w* (0.5+ (0.5 - 1./w))
        else:
            return (w - 0.5) - 0.5

# %%-------------------------------------- rlog
cdef inline double rlog(double x) noexcept nogil:
    cdef double r, t, u, w, w1
    cdef double a = .566749439387324e-01
    cdef double b = .456512608815524e-01
    cdef double[3] p = [.333333333333333, -.224696413112536,
                        .620886815375787e-02]
    cdef double[2] q = [-.127408923933623e+01, .354508718369557]

    if (x < 0.61) or (x > 1.57):
        return ((x - 0.5) - 0.5) - log(x)

    if x < 0.82:
        u = (x - 0.7) / 0.7
        w1 = a - u*0.3
    elif x > 1.18:
        u = 0.75*x - 1.
        w1 = b + u/3.
    else:
        u = (x - 0.5) - 0.5
        w1 = 0.

    r = u / (u + 2.)
    t = r*r
    w = ((p[2]*t+p[1])*t+p[0])/ ((q[1]*t+q[0])*t+1.)
    return 2.*t*(1. / (1. - r) - r*w) + w1

# %%-------------------------------------- rlog1
cdef inline double rlog1(double x) noexcept nogil:
    cdef:
        double a = .566749439387324e-01
        double b = .456512608815524e-01
        double p0 = .333333333333333e+00
        double p1 = -.224696413112536e+00
        double p2 = .620886815375787e-02
        double q1 = -.127408923933623e+01
        double q2 = .354508718369557e+00
        double r
        double t
        double w
        double w1

    if -0.39 <= x <= 0.57:
        if -0.18 <= x <= 0.18:
            h = x
            w1 = 0.

        elif x < -0.18:
            h = (x + 0.3)/0.7
            w1 = a - h*0.3

        else:  # 0.57 >= x > 0.18
            h = 0.75*x - 0.25
            w1 = b + h/3.

        r = h / (h + 2)
        t = r*r
        w = ((p2*t + p1)*t + p0) / ((q2*t + q1)*t + 1.)
        return 2*t*(1./(1.-r) - r*w) + w1
    else:
        return x - log((x + 0.5) + 0.5)

# %%-------------------------------------- stvaln
cdef inline double stvaln(double p) noexcept nogil:
    cdef double y, z
    cdef double[5] xnum = [-0.322232431088, -1.000000000000,
                           -0.342242088547, -0.204231210245e-1,
                           -0.453642210148e-4]
    cdef double[5] xden = [0.993484626060e-1, 0.588581570495,
                           0.531103462366, 0.103537752850,
                           0.38560700634e-2]

    z = 1. - p if (p > 0.5) else p
    y = sqrt(-2.*log(z))
    z = y + (devlpl(xnum, 5, y) / devlpl(xden, 5, y))
    return z if (p > 0.5) else -z
