import math
import numpy as np


#pythran export _Aij(float[:,:], int, int)
#pythran export _Aij(int[:,:], int, int)
def _Aij(A, i, j):
    """Sum of upper-left and lower right blocks of contingency table."""
    # See `somersd` References [2] bottom of page 309
    return A[:i, :j].sum() + A[i+1:, j+1:].sum()


#pythran export _Dij(float[:,:], int, int)
#pythran export _Dij(int[:,:], int, int)
def _Dij(A, i, j):
    """Sum of lower-left and upper-right blocks of contingency table."""
    # See `somersd` References [2] bottom of page 309
    return A[i+1:, :j].sum() + A[:i, j+1:].sum()


# pythran export _concordant_pairs(float[:,:])
# pythran export _concordant_pairs(int[:,:])
def _concordant_pairs(A):
    """Twice the number of concordant pairs, excluding ties."""
    # See `somersd` References [2] bottom of page 309
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*_Aij(A, i, j)
    return count


# pythran export _discordant_pairs(float[:,:])
# pythran export _discordant_pairs(int[:,:])
def _discordant_pairs(A):
    """Twice the number of discordant pairs, excluding ties."""
    # See `somersd` References [2] bottom of page 309
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*_Dij(A, i, j)
    return count


#pythran export _a_ij_Aij_Dij2(float[:,:])
#pythran export _a_ij_Aij_Dij2(int[:,:])
def _a_ij_Aij_Dij2(A):
    """A term that appears in the ASE of Kendall's tau and Somers' D."""
    # See `somersd` References [2] section 4:
    # Modified ASEs to test the null hypothesis...
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*(_Aij(A, i, j) - _Dij(A, i, j))**2
    return count


#pythran export _compute_outer_prob_inside_method(int64, int64, int64, int64)
def _compute_outer_prob_inside_method(m, n, g, h):
    """
    Count the proportion of paths that do not stay strictly inside two
    diagonal lines.

    Parameters
    ----------
    m : integer
        m > 0
    n : integer
        n > 0
    g : integer
        g is greatest common divisor of m and n
    h : integer
        0 <= h <= lcm(m,n)

    Returns
    -------
    p : float
        The proportion of paths that do not stay inside the two lines.

    The classical algorithm counts the integer lattice paths from (0, 0)
    to (m, n) which satisfy |x/m - y/n| < h / lcm(m, n).
    The paths make steps of size +1 in either positive x or positive y
    directions.
    We are, however, interested in 1 - proportion to computes p-values,
    so we change the recursion to compute 1 - p directly while staying
    within the "inside method" a described by Hodges.

    We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.
    Hodges, J.L. Jr.,
    "The Significance Probability of the Smirnov Two-Sample Test,"
    Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.

    For the recursion for 1-p see
    Viehmann, T.: "Numerically more stable computation of the p-values
    for the two-sample Kolmogorov-Smirnov test," arXiv: 2102.08037

    """
    # Probability is symmetrical in m, n.  Computation below uses m >= n.
    if m < n:
        m, n = n, m
    mg = m // g
    ng = n // g

    # Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    # |nx/g - my/g| < h.
    # Compute matrix A such that:
    #  A(x, 0) = A(0, y) = 1
    #  A(x, y) = A(x, y-1) + A(x-1, y), for x,y>=1, except that
    #  A(x, y) = 0 if |x/m - y/n|>= h
    # Probability is A(m, n)/binom(m+n, n)
    # Optimizations exist for m==n, m==n*p.
    # Only need to preserve a single column of A, and only a
    # sliding window of it.
    # minj keeps track of the slide.
    minj, maxj = 0, min(int(np.ceil(h / mg)), n + 1)
    curlen = maxj - minj
    # Make a vector long enough to hold maximum window needed.
    lenA = min(2 * maxj + 2, n + 1)
    # This is an integer calculation, but the entries are essentially
    # binomial coefficients, hence grow quickly.
    # Scaling after each column is computed avoids dividing by a
    # large binomial coefficient at the end, but is not sufficient to avoid
    # the large dynamic range which appears during the calculation.
    # Instead we rescale based on the magnitude of the right most term in
    # the column and keep track of an exponent separately and apply
    # it at the end of the calculation.  Similarly when multiplying by
    # the binomial coefficient
    dtype = np.float64
    A = np.ones(lenA, dtype=dtype)
    # Initialize the first column
    A[minj:maxj] = 0.0
    for i in range(1, m + 1):
        # Generate the next column.
        # First calculate the sliding window
        lastminj, lastlen = minj, curlen
        minj = max(int(np.floor((ng * i - h) / mg)) + 1, 0)
        minj = min(minj, n)
        maxj = min(int(np.ceil((ng * i + h) / mg)), n + 1)
        if maxj <= minj:
            return 1.0
        # Now fill in the values. We cannot use cumsum, unfortunately.
        val = 0.0 if minj == 0 else 1.0
        for jj in range(maxj - minj):
            j = jj + minj
            val = (A[jj + minj - lastminj] * i + val * j) / (i + j)
            A[jj] = val
        curlen = maxj - minj
        if lastlen > curlen:
            # Set some carried-over elements to 1
            A[maxj - minj:maxj - minj + (lastlen - curlen)] = 1

    return A[maxj - minj - 1]


# pythran export siegelslopes(float32[:], float32[:], str)
# pythran export siegelslopes(float64[:], float64[:], str)
def siegelslopes(y, x, method):
    deltax = np.expand_dims(x, 1) - x
    deltay = np.expand_dims(y, 1) - y
    slopes, intercepts = [], []

    for j in range(len(x)):
        id_nonzero, = np.nonzero(deltax[j, :])
        slopes_j = deltay[j, id_nonzero] / deltax[j, id_nonzero]
        medslope_j = np.median(slopes_j)
        slopes.append(medslope_j)
        if method == 'separate':
            z = y*x[j] - y[j]*x
            medintercept_j = np.median(z[id_nonzero] / deltax[j, id_nonzero])
            intercepts.append(medintercept_j)

    medslope = np.median(np.asarray(slopes))
    if method == "separate":
        medinter = np.median(np.asarray(intercepts))
    else:
        medinter = np.median(y - medslope*x)

    return medslope, medinter


# pythran export _poisson_binom_pmf(float64[:])
def _poisson_binom_pmf(p):
    # implemented from poisson_binom [2] Equation 2
    n = p.shape[0]
    pmf = np.zeros(n + 1, dtype=np.float64)
    pmf[:2] = 1 - p[0], p[0]
    for i in range(1, n):
        tmp = pmf[:i+1] * p[i]
        pmf[:i+1] *= (1 - p[i])
        pmf[1:i+2] += tmp
    return pmf


# pythran export _poisson_binom(int64[:], float64[:, :], str)
def _poisson_binom(k, args, tp):
    # PDF/CDF of Poisson binomial distribution
    # k - arguments, shape (m,)
    # args - shape parameters, shape (n, m)
    # kind - {'pdf', 'cdf'}
    n, m = args.shape  # number of shapes, batch size
    cache = {}
    out = np.zeros(m, dtype=np.float64)
    for i in range(m):
        p = tuple(args[:, i])
        if p not in cache:
            pmf = _poisson_binom_pmf(args[:, i])
            cache[p] = np.cumsum(pmf) if tp=='cdf' else pmf
        out[i] = cache[p][k[i]]
    return out


# function p = phid(z), p = erfc( -z/sqrt(2) )/2; % Normal cdf
def phid(z):
    return math.erfc(-z / math.sqrt(2)) / 2


def np_dot(x, y):
    return np.sum(x * y)


# function p = bvnu( dh, dk, r )
#pythran export _bvnu(float64, float64, float64)
def _bvnu(dh, dk, r):
    math_inf, math_pi = np.inf, np.pi
    # if dh ==  inf | dk ==  inf:p = 0;
    if (dh == math_inf) or (dk == math_inf):
        p = 0.
    # elseif dh == -inf, if dk == -inf, p = 1; else p = phid(-dk); end
    elif dh == -math_inf:
        if dk == -math_inf:
            p = 1.
        else:
            p = phid(-dk)
    # elseif dk == -inf, p = phid(-dh);
    elif dk == -math_inf:
        p = phid(-dh)
    # elseif r == 0, p = phid(-dh)*phid(-dk);
    elif r == 0:
        p = phid(-dh) * phid(-dk)
    # else, tp = 2*pi; h = dh; k = dk; hk = h*k; bvn = 0;
    else:
        tp = 2*math_pi
        h = dh
        k = dk
        hk = h*k
        bvn = 0.
        # if abs(r) < 0.3      % Gauss Legendre points and weights, n =  6
        #     w(1:3) = [0.1713244923791705 0.3607615730481384 0.4679139345726904];
        #     x(1:3) = [0.9324695142031522 0.6612093864662647 0.2386191860831970];
        if abs(r) < 0.3:
            w = [0.1713244923791705, 0.3607615730481384, 0.4679139345726904]
            x = [0.9324695142031522, 0.6612093864662647, 0.2386191860831970]
        # elseif abs(r) < 0.75 % Gauss Legendre points and weights, n = 12
        #     w(1:3) = [.04717533638651177 0.1069393259953183 0.1600783285433464];
        #     w(4:6) = [0.2031674267230659 0.2334925365383547 0.2491470458134029];
        #     x(1:3) = [0.9815606342467191 0.9041172563704750 0.7699026741943050];
        #     x(4:6) = [0.5873179542866171 0.3678314989981802 0.1252334085114692];
        elif abs(r) < 0.75:
            w = [.04717533638651177, 0.1069393259953183, 0.1600783285433464,
                 0.2031674267230659, 0.2334925365383547, 0.2491470458134029]
            x = [0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
                 0.5873179542866171, 0.3678314989981802, 0.1252334085114692]
        # else,                % Gauss Legendre points and weights, n = 20
        #     w(1:3) = [.01761400713915212 .04060142980038694 .06267204833410906];
        #     w(4:6) = [.08327674157670475 0.1019301198172404 0.1181945319615184];
        #     w(7:9) = [0.1316886384491766 0.1420961093183821 0.1491729864726037];
        #     w(10) =   0.1527533871307259;
        #     x(1:3) = [0.9931285991850949 0.9639719272779138 0.9122344282513259];
        #     x(4:6) = [0.8391169718222188 0.7463319064601508 0.6360536807265150];
        #     x(7:9) = [0.5108670019508271 0.3737060887154196 0.2277858511416451];
        #     x(10) =   0.07652652113349733;
        else:
            w = [.01761400713915212, .04060142980038694, .06267204833410906,
                 .08327674157670475, 0.1019301198172404, 0.1181945319615184,
                 0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
                 0.1527533871307259]
            x = [0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
                 0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
                 0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
                 0.07652652113349733]
        # end, w = [w  w]; x = [1-x 1+x];
        w, x = np.asarray(w), np.asarray(x)
        w = np.concatenate((w, w))
        x = np.concatenate((1-x, 1+x))
        # if abs(r) < 0.925, hs = ( h*h + k*k )/2; asr = asin(r)/2;
        if abs(r) < 0.925:
            hs = ( h*h + k*k )/2
            asr = math.asin(r)/2
        #     sn = sin(asr*x); bvn = exp((sn*hk-hs)./(1-sn.^2))*w';
            sn = np.sin(asr*x)
            bvn = np_dot(np.exp((sn*hk-hs) / (1-sn**2)), w)
        #     bvn = bvn*asr/tp + phid(-h)*phid(-k);
            bvn = bvn*asr/tp + phid(-h)*phid(-k)
        # else, if r < 0, k = -k; hk = -hk; end
        else:
            if r < 0:
                k = -k
                hk = -hk
            # if abs(r) < 1, as = 1-r^2; a = sqrt(as); bs = (h-k)^2;
            if abs(r) < 1:
                as_ = 1-r**2
                a = math.sqrt(as_)
                bs = (h-k)**2
                # asr = -( bs/as + hk )/2; c = (4-hk)/8 ; d = (12-hk)/80;
                asr = -( bs/as_ + hk )/2
                c = (4-hk)/8
                d = (12-hk)/80
                # if asr > -100, bvn = a*exp(asr)*(1-c*(bs-as)*(1-d*bs)/3+c*d*as^2); end
                if asr > -100:
                    bvn = a*math.exp(asr)*(1-c*(bs-as_)*(1-d*bs)/3+c*d*as_**2)
                # if hk  > -100, b = sqrt(bs); sp = sqrt(tp)*phid(-b/a);
                if hk  > -100:
                    b = math.sqrt(bs)
                    sp = math.sqrt(tp)*phid(-b/a)
                    # bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1-d*bs)/3 );
                    bvn = bvn - math.exp(-hk/2)*sp*b*( 1 - c*bs*(1-d*bs)/3 )

                # end, a = a/2; xs = (a*x).^2; asr = -( bs./xs + hk )/2;
                a = a/2
                xs = (a*x)**2
                asr = -( bs / xs + hk )/2
                # ix = find( asr > -100 ); xs = xs(ix); sp = ( 1 + c*xs.*(1+5*d*xs) );
                ix = asr > -100
                xs = xs[ix]
                sp = 1 + c*xs * (1+5*d*xs)
                # rs = sqrt(1-xs); ep = exp( -(hk/2)*xs./(1+rs).^2 )./rs;
                rs = np.sqrt(1-xs)
                ep = np.exp( -(hk/2)*xs / (1+rs)**2 )/rs
                # bvn = ( a*( (exp(asr(ix)).*(sp-ep))*w(ix)' ) - bvn )/tp;
                bvn = ( a*np_dot( (np.exp(asr[ix]) * (sp-ep)), w[ix] ) - bvn )/tp
            # end
            # if r > 0, bvn =  bvn + phid( -max( h, k ) );
            if r > 0:
                bvn =  bvn + phid( -max( h, k ) )
            # elseif h >= k, bvn = -bvn;
            elif h >= k:
                bvn = -bvn
            # else, if h < 0, L = phid(k)-phid(h); else, L = phid(-h)-phid(-k); end
            else:
                if h < 0:
                    L = phid(k)-phid(h)
                else:
                    L = phid(-h)-phid(-k)
                # bvn =  L - bvn;
                bvn =  L - bvn
            # end
        # end, p = max( 0, min( 1, bvn ) );
        p = max( 0, min( 1, bvn ) )
    # end
    return p
