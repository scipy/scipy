"""This script tests scipy's implementation of hyp2f1 against mpmath's.

Author: Albert Steppi, with credit to Adam Kullberg (FormerPhycisist) for
the implementation of mp_hyp2f1 below, which modifies mpmath's hyp2f1 to
return the same branch as scipy's on the standard branch cut.

Produces a tab separated values file with 11 columns. The first four columns
contain the parameters a, b, c and the argument z. The next two contain |z| and
a region code for which region of the complex plane belongs to. The regions are

    1) |z| < 0.9 and real(z) >= 0
    2) |z| < 1 and real(z) < 0
    3) 0.9 <= |z| <= 1 and |1 - z| < 1.0:
    4) 0.9 <= |z| <= 1 and |1 - z| >= 1.0:
    5) |z| > 1

Parameters a, b, c are taken from a 10 * 10 * 10 grid with values at

    -16, -8, -4, -2, -1, 1, 2, 4, 8, 16

with random perturbations applied. The following cases are handled.

    1) A, B, C, B - A, C - A, C - B, C - A - B all non-integral.
    2) B - A integral
    3) C - A integral
    4) C - B integral
    5) C - A - B integral
    6) A integral
    7) B integral
    8) C integral

The seventh column of the output file is an integer between 1 and 8 specifying
the parameter group as above.

The argument z is taken from a 20 * 20 grid in the box
    -2 <= real(z) <= 2, -2 <= imag(z) <= 2.

The final four columns have the expected value of hyp2f1 for the given
parameters and argument as calculated with mpmath, the observed value
calculated with scipy's hyp2f1, the relative error, and the absolute error.

As special cases of hyp2f1 are moved from the original Fortran implementation
into Cython, this script can be used to ensure that no regressions occur and
to point out where improvements are needed.

The generated file is roughly 700MB in size. The script has two arguments;
a positional argument for specifying the path to the location where the output
file is to be placed, and an optional argument specifying the number of
processes to use for parallel execution.

Takes around 40 minutes using an Intel(R) Core(TM) i5-8250U CPU with n_jobs
set to 8 (full utilization).
"""


import os
import csv
import argparse
import numpy as np
from itertools import product
from multiprocessing import Pool


from scipy.special import hyp2f1


try:
    import mpmath
except ImportError:
    pass


def mp_hyp2f1(a, b, c, z):
    """Return mpmath hyp2f1 calculated on same branch as scipy hyp2f1.

    For most values of a,b,c mpmath returns the x - 0j branch of hyp2f1 on the
    branch cut x=(1,inf) whereas scipy's hyp2f1 calculates the x + 0j branch.
    Thus, to generate the right comparison values on the branch cut, we
    evaluate mpmath.hyp2f1 at x + 1e-15*j.

    The exception to this occurs when c-a=-m in which case both mpmath and
    scipy calculate the x + 0j branch on the branch cut. When this happens
    mpmath.hyp2f1 will be evaluated at the original z point.
    """
    on_branch_cut = z.real > 1.0 and abs(z.imag) < 1.0e-15
    cond1 = abs(c - a - round(c - a)) < 1.0e-15 and round(c - a) <= 0
    cond2 = abs(c - b - round(c - b)) < 1.0e-15 and round(c - b) <= 0
    # Make sure imaginary part is *exactly* zero
    if on_branch_cut:
        z = z.real + 0.0j
    if on_branch_cut and not (cond1 or cond2):
        z_mpmath = z.real + 1.0e-15j
    else:
        z_mpmath = z
    return complex(mpmath.hyp2f1(a, b, c, z_mpmath))


def get_region(z):
    """Assign numbers for regions where hyp2f1 is calculated differently."""
    if abs(z) < 0.9 and z.real >= 0:
        return 1
    elif abs(z) < 1 and z.real < 0:
        return 2
    elif 0.9 <= abs(z) <= 1 and abs(1 - z) < 1.0:
        return 3
    elif 0.9 <= abs(z) <= 1 and abs(1 - z) >= 1.0:
        return 4
    else:
        return 5


def get_result(a, b, c, z, group):
    """Get results for given parameter and value combination."""
    expected, observed = mp_hyp2f1(a, b, c, z), hyp2f1(a, b, c, z)
    if expected != 0:
        relative_error = abs(expected - observed) / abs(expected)
    else:
        relative_error = float('nan')
    return (
        a,
        b,
        c,
        z,
        abs(z),
        get_region(z),
        group,
        expected,
        observed,
        relative_error,
        abs(expected - observed),
    )


def get_results(params, Z, n_jobs=1):
    """Batch compute results for multiple parameter and argument values.

    Parameters
    ----------
    params : iterable
        iterable of tuples of floats (a, b, c) specificying parameter values
        a, b, c for hyp2f1
    Z : iterable of complex
        Arguments at which to evaluate hyp2f1
    n_jobs : Optional[int]
        Number of jobs for parallel execution.

    Returns
    -------
    list
        List of tuples of results values. See return value in source code
        of `get_result`.
    """
    input_ = (
        (a, b, c, z, group) for (a, b, c, group), z in product(params, Z)
    )

    with Pool(n_jobs) as pool:
        rows = pool.starmap(get_result, input_)
    return rows


def _make_hyp2f1_test_case(a, b, c, z, rtol):
    """Generate string for single test case as used in test_hyp2f1.py."""
    expected = mp_hyp2f1(a, b, c, z)
    return (
        "    pytest.param(\n"
        "        Hyp2f1TestCase(\n"
        f"            a={a},\n"
        f"            b={b},\n"
        f"            c={c},\n"
        f"            z={z},\n"
        f"            expected={expected},\n"
        f"            rtol={rtol},\n"
        "        ),\n"
        "    ),"
    )


def make_hyp2f1_test_cases(rows):
    """Generate string for a list of test cases for test_hyp2f1.py.

    Parameters
    ----------
    rows : list
        List of lists of the form [a, b, c, z, rtol] where a, b, c, z are
        parameters and the argument for hyp2f1 and rtol is an expected
        relative error for the associated test case.

    Returns
    -------
    str
        String for a list of test cases. The output string can be printed
        or saved to a file and then copied into an argument for
        `pytest.mark.parameterize` within `scipy.special.tests.test_hyp2f1.py`.
    """
    result = "[\n"
    result += '\n'.join(
        _make_hyp2f1_test_case(a, b, c, z, rtol)
        for a, b, c, z, rtol in rows
    )
    result += "\n]"
    return result


def main(outpath, n_jobs=1):
    outpath = os.path.realpath(os.path.expanduser(outpath))

    random_state = np.random.RandomState(1234)
    root_params = np.array(
        [-16, -8, -4, -2, -1, 1, 2, 4, 8, 16]
    )

    perturbations = 0.1 * random_state.random_sample(
        size=(3, len(root_params))
    )

    params = []
    # No integer differences
    A = root_params + perturbations[0, :]
    B = root_params + perturbations[1, :]
    C = root_params + perturbations[2, :]
    params.extend(
        sorted(
            ((a, b, c, 1) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # B - A an integer
    A = root_params + 0.5
    B = root_params + 0.5
    C = root_params + perturbations[1, :]
    params.extend(
        sorted(
            ((a, b, c, 2) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # C - A an integer
    A = root_params + 0.5
    B = root_params + perturbations[1, :]
    C = root_params + 0.5
    params.extend(
        sorted(
            ((a, b, c, 3) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # C - B an integer
    A = root_params + perturbations[0, :]
    B = root_params + 0.5
    C = root_params + 0.5
    params.extend(
        sorted(
            ((a, b, c, 4) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # C - A - B an integer
    A = root_params + 0.25
    B = root_params + 0.25
    C = root_params + 0.5
    params.extend(
        sorted(
            ((a, b, c, 5) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # A an integer
    A = root_params
    B = root_params + perturbations[0, :]
    C = root_params + perturbations[1, :]
    params.extend(
        sorted(
            ((a, b, c, 6) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # B an integer
    A = root_params + perturbations[0, :]
    B = root_params
    C = root_params + perturbations[1, :]
    params.extend(
        sorted(
            ((a, b, c, 7) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    # C an integer
    A = root_params + perturbations[0, :]
    B = root_params + perturbations[1, :]
    C = root_params
    params.extend(
        sorted(
            ((a, b, c, 8) for a, b, c in product(A, B, C)),
            key=lambda x: max(abs(x[0]), abs(x[1])),
        )
    )

    X, Y = np.meshgrid(np.linspace(-2, 2, 20), np.linspace(-2, 2, 20))
    Z = X + Y * 1j
    Z = Z.flatten()

    rows = get_results(params, Z, n_jobs=n_jobs)

    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "a",
                "b",
                "c",
                "z",
                "|z|",
                "region",
                "parameter_group",
                "expected",
                "observed",
                "relative_error",
                "absolute_error",
            ]
        )
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test hyp2f1 against mpmath on a grid in the complex plane"
        " over a grid of parameter values."
    )
    parser.add_argument(
        "outpath", type=str, help="Path to output tsv file."
    )
    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="Number of jobs for multiprocessing.",
    )
    args = parser.parse_args()
    main(args.outpath, n_jobs=args.n_jobs)
