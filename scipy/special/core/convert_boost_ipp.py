"""Script to convert .ipp test data from boost to csv format. The data are
found in boost-trunk/libs/math/test."""

import re

from numpy import arccosh, arcsinh, arctanh, log1p, expm1
from scipy.special import gamma, gammaln, digamma, erf, erfc, expi, expn

CXX_COMMENT = re.compile(r'^\s+//')
DATA_REGEX = re.compile(r'^\s+/*\{*\s*SC_')
ITEM_REGEX = re.compile(r'SC_\(([-\d\.e]+)\)')
HEADER_REGEX = re.compile(
r'const boost::array\<boost::array\<T, (\d+)\>, (\d+)\> ([a-zA-Z_\d]+)')

def id(x):
    return x

TEST_DATA = {
    'acosh_data.ipp': (id, arccosh),
    'asinh_data.ipp': (id, arcsinh),
    'assoc_legendre_p.ipp': (id, ),
    'atanh_data.ipp': (id, arctanh),
    'bessel_i_data.ipp': (id, ),
    'bessel_i_int_data.ipp': (id, ),
    'bessel_j_data.ipp': (id, ),
    'bessel_j_int_data.ipp': (id, ),
    'bessel_j_large_data.ipp': (id, ),
    'bessel_k_data.ipp': (id, ),
    'bessel_k_int_data.ipp': (id, ),
    'bessel_y01_data.ipp': (id, ),
    'bessel_yn_data.ipp': (id, ),
    'bessel_yv_data.ipp': (id, ),
    'beta_exp_data.ipp': (id, ),
    'beta_med_data.ipp': (id, ),
    'beta_small_data.ipp': (id, ),
    'binomial_data.ipp': (id, ),
    'binomial_large_data.ipp': (id, ),
    'binomial_quantile.ipp': (id, ),
    'cbrt_data.ipp': (id, ),
    'digamma_data.ipp': (id, digamma),
    'digamma_neg_data.ipp': (id, digamma),
    'digamma_root_data.ipp': (id, digamma),
    'digamma_small_data.ipp': (id, digamma),
    'ellint_e2_data.ipp': (id, ),
    'ellint_e_data.ipp': (id, ),
    'ellint_f_data.ipp': (id, ),
    'ellint_k_data.ipp': (id, ),
    'ellint_pi2_data.ipp': (id, ),
    'ellint_pi3_data.ipp': (id, ),
    'ellint_pi3_large_data.ipp': (id, ),
    'ellint_rc_data.ipp': (id, ),
    'ellint_rd_data.ipp': (id, ),
    'ellint_rf_data.ipp': (id, ),
    'ellint_rj_data.ipp': (id, ),
    'erfc_inv_big_data.ipp': (id, erfc),
    'erfc_inv_data.ipp': (id, erfc),
    'erf_data.ipp': (id, erf, erfc),
    'erf_inv_data.ipp': (id, erf),
    'erf_large_data.ipp': (id, erf, erfc),
    'erf_small_data.ipp': (id, erf, erfc),
    'expint_1_data.ipp': (id, expi),
    'expint_data.ipp': (id, expi),
    'expinti_data_double.ipp': (id, expi),
    'expinti_data.ipp': (id, expi),
    'expinti_data_long.ipp': (id, expi),
    'expint_small_data.ipp': (id, expi, expn),
    'gamma_inv_big_data.ipp': (id, ),
    'gamma_inv_data.ipp': (id, ),
    'gamma_inv_small_data.ipp': (id, ),
    'hermite.ipp': (id, ),
    'ibeta_data.ipp': (id, ),
    'ibeta_int_data.ipp': (id, ),
    'ibeta_inva_data.ipp': (id, ),
    'ibeta_inv_data.ipp': (id, ),
    'ibeta_large_data.ipp': (id, ),
    'ibeta_small_data.ipp': (id, ),
    'igamma_big_data.ipp': (id, ),
    'igamma_int_data.ipp': (id, ),
    'igamma_inva_data.ipp': (id, ),
    'igamma_med_data.ipp': (id, ),
    'igamma_small_data.ipp': (id, ),
    'laguerre2.ipp': (id, ),
    'laguerre3.ipp': (id, ),
    'legendre_p.ipp': (id, ),
    'legendre_p_large.ipp': (id, ),
    'log1p_expm1_data.ipp': (id, log1p, expm1),
    'ncbeta_big.ipp': (id, ),
    'ncbeta.ipp': (id, ),
    #'nccs_big.ipp': (id, ),
    'nccs.ipp': (id, ),
    'nct.ipp': (id, ),
    'negative_binomial_quantile.ipp': (id, ),
    'poisson_quantile.ipp': (id, ),
    'sph_bessel_data.ipp': (id, ),
    'spherical_harmonic.ipp': (id, ),
    'sph_neumann_data.ipp': (id, ),
    'test_gamma_data.ipp': (id, gamma, gammaln),
    'tgamma_delta_ratio_data.ipp': (id, ),
    'tgamma_delta_ratio_int2.ipp': (id, ),
    'tgamma_delta_ratio_int.ipp': (id, ),
    'tgamma_ratio_data.ipp': (id, ),
    'zeta_1_below_data.ipp': (id, ),
    'zeta_1_up_data.ipp': (id, ),
    'zeta_data.ipp': (id, ),
    'zeta_neg_data.ipp': (id, ),
}

def _raw_data(line):
    items = line.split(',')
    l = []
    for item in items:
        m = ITEM_REGEX.search(item)
        if m:
            l.append(m.group(1))
    return l

def get_data(filename):
    a = open(filename, 'r')
    lines = a.readlines()
    data = {}
    i = 0
    while (i < len(lines)):
        line = lines[i]
        m = HEADER_REGEX.search(line)
        if m:
            d = int(m.group(1))
            n = int(m.group(2))
            cdata = []
            i += 1
            line = lines[i]
            # Skip comments
            while CXX_COMMENT.match(line):
                i += 1
                line = lines[i]
            while DATA_REGEX.match(line):
                cdata.append(_raw_data(line))
                i += 1
                line = lines[i]
                # Skip comments
                while CXX_COMMENT.match(line):
                    i += 1
                    line = lines[i]
            if not len(cdata) == n:
                raise ValueError("parsed data: %d, expected %d" % (len(cdata), n))
            data[m.group(3)] = cdata
        else:
            i += 1

    return data

if __name__ == '__main__':
    for file, funcs in TEST_DATA.items():
        print "================= %s ===============" % file
        data = get_data(file)
        print "%d data sets" % len(data)
        for k, d in data.items():
            fid = open('data/%s.txt' % k, 'w')
            for line in d:
                fid.write("%s,\n" % ",".join(line))

    #    for items in data:
    #        assert len(items) == len(funcs)
    #        for i in range(1, len(items)):
    #            print funcs[i](float(items[0])) - float(items[i])
    #file = 'poisson_quantile.ipp'
    #data = get_data(file)
