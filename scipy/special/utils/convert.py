# This script is used to parse BOOST special function test data into something
# we can easily import in numpy. It is ugly as hell, but it works.
from __future__ import division, print_function, absolute_import

import re
import os

# Where to put the data (directory will be created)
DATA_DIR = 'data'
# where to pull out boost data: assume a SVN checkout of boost (here in
# /Users/david/src/dev/boost/trunk)
BOOST_SRC = "boostmath/test"

CXX_COMMENT = re.compile(r'^\s+//')
DATA_REGEX = re.compile(r'^\s+/*\{*\s*SC_')
ITEM_REGEX = re.compile(r'[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?')
HEADER_REGEX = re.compile(
r'const boost::array\<boost::array\<.*, (\d+)\>, (\d+)\> ([a-zA-Z_\d]+)')

# List of boost test data files to parse
DATA_FILES = [
    'acosh_data.ipp',
    'asinh_data.ipp',
    'assoc_legendre_p.ipp',
    'atanh_data.ipp',
    'bessel_i_data.ipp',
    'bessel_i_int_data.ipp',
    'bessel_j_data.ipp',
    'bessel_j_int_data.ipp',
    'bessel_j_large_data.ipp',
    'bessel_k_data.ipp',
    'bessel_k_int_data.ipp',
    'bessel_y01_data.ipp',
    'bessel_yn_data.ipp',
    'bessel_yv_data.ipp',
    'beta_exp_data.ipp',
    'beta_med_data.ipp',
    'beta_small_data.ipp',
    'binomial_data.ipp',
    'binomial_large_data.ipp',
    'binomial_quantile.ipp',
    'cbrt_data.ipp',
    'digamma_data.ipp',
    'digamma_neg_data.ipp',
    'digamma_root_data.ipp',
    'digamma_small_data.ipp',
    'ellint_e2_data.ipp',
    'ellint_e_data.ipp',
    'ellint_f_data.ipp',
    'ellint_k_data.ipp',
    'ellint_pi2_data.ipp',
    'ellint_pi3_data.ipp',
    'ellint_pi3_large_data.ipp',
    'ellint_rc_data.ipp',
    'ellint_rd_data.ipp',
    'ellint_rf_data.ipp',
    'ellint_rj_data.ipp',
    'erfc_inv_big_data.ipp',
    'erfc_inv_data.ipp',
    'erf_data.ipp',
    'erf_inv_data.ipp',
    'erf_large_data.ipp',
    'erf_small_data.ipp',
    'expint_1_data.ipp',
    'expint_data.ipp',
    'expinti_data_double.ipp',
    'expinti_data.ipp',
    'expinti_data_long.ipp',
    'expint_small_data.ipp',
    'gamma_inv_big_data.ipp',
    'gamma_inv_data.ipp',
    'gamma_inv_small_data.ipp',
    'hermite.ipp',
    'hypergeometric_dist_data2.ipp',
    'hypergeometric_test_data.ipp',
    'ibeta_data.ipp',
    'ibeta_int_data.ipp',
    'ibeta_inva_data.ipp',
    'ibeta_inv_data.ipp',
    'ibeta_large_data.ipp',
    'ibeta_small_data.ipp',
    'igamma_big_data.ipp',
    'igamma_int_data.ipp',
    'igamma_inva_data.ipp',
    'igamma_med_data.ipp',
    'igamma_small_data.ipp',
    'laguerre2.ipp',
    'laguerre3.ipp',
    'legendre_p.ipp',
    'legendre_p_large.ipp',
    'log1p_expm1_data.ipp',
    'ncbeta_big.ipp',
    'ncbeta.ipp',
    'nccs_big.ipp',
    'nccs.ipp',
    'nct.ipp',
    'negative_binomial_quantile.ipp',
    'poisson_quantile.ipp',
    'powm1_sqrtp1m1_test.cpp',
    'sph_bessel_data.ipp',
    'spherical_harmonic.ipp',
    'sph_neumann_data.ipp',
    #'test_bessel_i.cpp',
    #'test_bessel_j.cpp',
    #'test_bessel_k.cpp',
    #'test_bessel_y.cpp',
    # Those 3 files use arithmetic operations whithin the data, so we can't parse
    # them naively
    #'test_ellint_1.cpp',
    #'test_ellint_2.cpp',
    #'test_ellint_3.cpp',
    'test_gamma_data.ipp',
    'tgamma_delta_ratio_data.ipp',
    'tgamma_delta_ratio_int2.ipp',
    'tgamma_delta_ratio_int.ipp',
    'tgamma_ratio_data.ipp',
    'zeta_1_below_data.ipp',
    'zeta_1_up_data.ipp',
    'zeta_data.ipp',
    'zeta_neg_data.ipp',
]


def _raw_data(line):
    items = line.split(',')
    l = []
    for item in items:
        m = ITEM_REGEX.search(item)
        if m:
            q = m.group(0)
            l.append(q)
    return l


def parse_ipp_file(filename):
    print(filename)
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
            print("d = {0}, n = {1}".format(d, n))
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


def dump_dataset(filename, data):
    fid = open(filename, 'w')
    try:
        for line in data:
            fid.write("%s\n" % " ".join(line))
    finally:
        fid.close()


def dump_datasets(filename):
    base, ext = os.path.splitext(os.path.basename(filename))
    base += '_%s' % ext[1:]
    datadir = os.path.join(DATA_DIR, base)
    os.makedirs(datadir)
    datasets = parse_ipp_file(filename)
    for k, d in datasets.items():
        print(k, len(d))
        dfilename = os.path.join(datadir, k) + '.txt'
        dump_dataset(dfilename, d)

if __name__ == '__main__':
    for filename in DATA_FILES:
        filename = os.path.join(BOOST_SRC, filename)
        print("================= %s ===============" % filename)
        dump_datasets(filename)
