# This script is used to parse BOOST special function test data into something
# we can easily import in numpy.
import re
import os

# Where to put the data (directory will be created)
DATA_DIR = 'scipy/special/tests/data/boost'
# Where to pull out boost data
BOOST_SRC = "boostmath/test"

CXX_COMMENT = re.compile(r'^\s+//')
DATA_REGEX = re.compile(r'^\s*/*\{*\s*SC_')
ITEM_REGEX = re.compile(r'[+-]?\d*\.?\d+(?:[eE][+-]?\d+)?')
HEADER_REGEX = re.compile(
r'const boost::array\<boost::array\<.*, (\d+)\>, (\d+)\> ([a-zA-Z_\d]+)')

IGNORE_PATTERNS = [
    # Makes use of ldexp and casts
    "hypergeometric_1F1_big_double_limited.ipp",
    "hypergeometric_1F1_big_unsolved.ipp",

    # Makes use of numeric_limits and ternary operator
    "beta_small_data.ipp",

    # Doesn't contain any data
    "almost_equal.ipp",

    # Derivatives functions don't exist
    "bessel_y01_prime_data.ipp",
    "bessel_yn_prime_data.ipp",
    "sph_bessel_prime_data.ipp",
    "sph_neumann_prime_data.ipp",

    # Data files not needed by scipy special tests.
    "ibeta_derivative_",
    r"ellint_r[cdfjg]_[^d]",
    r"ellint_d2?_",
    "jacobi_",
    "heuman_lambda_",
    "hypergeometric_",
    "nct_",
    r".*gammap1m1_",
    "trig_",
    "powm1_data.ipp",
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
    with open(filename, 'r') as a:
        lines = a.readlines()
    data = {}
    i = 0
    while (i < len(lines)):
        line = lines[i]
        m = HEADER_REGEX.search(line)
        if m:
            d = int(m.group(1))
            n = int(m.group(2))
            print(f"d = {d}, n = {n}")
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
                raise ValueError(f"parsed data: {len(cdata)}, expected {n}")
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
    for filename in sorted(os.listdir(BOOST_SRC)):
        # Note: Misses data in hpp files (e.x. powm1_sqrtp1m1_test.hpp)
        if filename.endswith(".ipp"):
            if any(re.match(pattern, filename) for pattern in IGNORE_PATTERNS):
                continue

            path = os.path.join(BOOST_SRC, filename)
            print(f"================= {path} ===============")
            dump_datasets(path)
