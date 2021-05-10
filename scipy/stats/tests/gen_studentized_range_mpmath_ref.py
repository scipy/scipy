# To run this script, run
# `python gen_studentized_range_mpmath_ref.py`
# in the "scipy/stats/tests/" directory

# This script generates a JSON file "./data/studentized_range_mpmath_ref.json"
# that is used to compare the accuracy of `studentized_range` functions against
# precise (20 DOP) results generated using `mpmath`.

# Note: I would have prefered to use pickle rather than JSON, but -
# due to security concerns - decided against it.

from collections import namedtuple
import json
import time

from mpmath import gamma, pi, sqrt, quad, inf, mpf, mp, exp
from mpmath import npdf as phi
from mpmath import ncdf as Phi

results_filepath = "data/studentized_range_mpmath_ref.json"

MPResult = namedtuple("MPResult", ["src_case", "mp_result"])

CdfCase = namedtuple("CdfCase",
                     ["q", "k", "v", "expected_atol", "expected_rtol"])

MomentCase = namedtuple("MomentCase",
                        ["m", "k", "v", "expected_atol", "expected_rtol"])

general_atol = 1e-11
general_rtol = 1e-11
mp.dps = 24

cdf_cases = [
    CdfCase(0.1, 3, 10, general_atol, general_rtol),
    CdfCase(1, 3, 10, general_atol, general_rtol),
    CdfCase(10, 3, 10, general_atol, general_rtol),

    CdfCase(0.1, 10, 10, general_atol, general_rtol),
    CdfCase(1, 10, 10, general_atol, general_rtol),
    CdfCase(10, 10, 10, general_atol, general_rtol),

    CdfCase(0.1, 3, 100, general_atol, general_rtol),
    CdfCase(1, 3, 100, general_atol, general_rtol),
    CdfCase(10, 3, 100, general_atol, general_rtol),

    CdfCase(0.1, 10, 100, general_atol, general_rtol),
    CdfCase(1, 10, 100, general_atol, general_rtol),
    CdfCase(10, 10, 100, general_atol, general_rtol),

    CdfCase(1, 3, 1000, general_atol, general_rtol),
    CdfCase(10, 3, 1000, 1e-10, general_rtol),
]

# PDF can use same cases as CDF
pdf_cases = cdf_cases

mom_atol, mom_rtol = 1e-9, 1e-9
# These are EXTREMELY slow - Multiple days each in worst case.
moment_cases = [
    MomentCase(1, 3, 10, mom_atol, mom_rtol),
    MomentCase(2, 3, 10, mom_atol, mom_rtol),
    MomentCase(3, 3, 10, mom_atol, mom_rtol),
    MomentCase(4, 3, 10, mom_atol, mom_rtol)
]


def to_dict(named_tuple):
    return dict(named_tuple._asdict())


def mp_res_to_dict(mp_result):
    """Formats an MPResult namedtuple into a dict for JSON dumping"""
    return {
        "src_case": to_dict(mp_result.src_case),

        # np assert can't handle mpf, so take the accuracy hit here.
        "mp_result": float(mp_result.mp_result)
    }


def cdf_mp(q, k, nu):
    """Calculates the CDF integrand for the given q, k, nu vals"""
    q, k, nu = mpf(q), mpf(k), mpf(nu)

    def inner(s, z):
        return phi(z) * (Phi(z + q * s) - Phi(z)) ** (k - 1)

    def outer(s, z):
        return s ** (nu - 1) * phi(sqrt(nu) * s) * inner(s, z)

    def whole(s, z):
        return (sqrt(2 * pi) * k * nu ** (nu / 2) /
                (gamma(nu / 2) * 2 ** (nu / 2 - 1)) * outer(s, z))

    res = quad(whole, [0, inf], [-inf, inf],
               method="gauss-legendre", maxdegree=10)
    return res


def pdf_mp(q, k, nu):
    """Calculates the PDF integrand for the given q, k, nu vals"""
    q, k, nu = mpf(q), mpf(k), mpf(nu)

    def integral(s, z):
        return (k * (k - 1) * s * phi(z) * phi(s * q + z)
                * (Phi(s * q + z) - Phi(z)) ** (k - 2))

    def whole(s, z):
        return (nu ** (nu / 2) / (gamma(nu / 2) * 2 ** (nu / 2 - 1)) *
                s ** (nu - 1) * exp(-nu * s ** 2 / 2) * integral(s, z))

    res = quad(whole, [0, inf], [-inf, inf],
               method="gauss-legendre", maxdegree=10)
    return res


def moment_mp(m, k, nu):
    """Calculates the moment integrand for the given m, k, nu vals"""
    m, k, nu = mpf(m), mpf(k), mpf(nu)

    def integral(q, s, z):
        return (k * (k - 1) * s * phi(z) * phi(s * q + z)
                * (Phi(s * q + z) - Phi(z)) ** (k - 2))

    def whole(q, s, z):
        return (q ** m * nu ** (nu / 2) / (gamma(nu / 2) * 2 ** (nu / 2 - 1))
                * s ** (nu - 1) * exp(-nu * s ** 2 / 2) * integral(q, s, z))

    res = quad(whole, [0, inf], [0, inf], [-inf, inf],
               method="gauss-legendre", maxdegree=10)
    return res


total_cases = len(cdf_cases) + len(pdf_cases)
print(f"Processing {total_cases} test cases")


def run(case, run_lambda, index=0, total_cases=0):
    """Runs the single passed case. Returns a result dict"""
    t_start = time.perf_counter()
    res = run_lambda(case)
    print(f"Finished {index + 1}/{total_cases} in batch. "
          f"(Took {time.perf_counter() - t_start}s)")

    return mp_res_to_dict(MPResult(case, res))


def run_cases(cases, run_lambda):
    """Runs an array of cases, returns an array of dicts"""
    return [run(case, run_lambda, index, len(cases))
            for index, case in enumerate(cases)]


def run_pdf(case):
    return pdf_mp(case.q, case.k, case.v)


def run_cdf(case):
    return cdf_mp(case.q, case.k, case.v)


def run_moment(case):
    return moment_mp(case.m, case.k, case.v)


t_start = time.perf_counter()

print(f"Running 1st batch ({len(pdf_cases)} PDF cases). "
      f"These take about 30s each.")
pdf_data = run_cases(pdf_cases, run_pdf)

print(f"Running 2nd batch ({len(pdf_cases)} CDF cases). "
      f"These take about 30s each.")
cdf_data = run_cases(cdf_cases, run_cdf)

print(f"Running 3rd batch ({len(moment_cases)} moment cases). "
      f"These take about anywhere from a few hours to days each.")
moment_data = run_cases(moment_cases, run_moment)

print(f"Test data generated in {time.perf_counter() - t_start}s")

# store data with the function type as a top level key to allow future
# expansion
result_dict = {
    "COMMENT": "!!!!!!!!!! THIS FILE WAS AUTOGENERATED BY RUNNING "
               "`python gen_studentized_range_mpmath_ref.py` !!!!!!!!!!",
    "cdf_data": cdf_data,
    "pdf_data": pdf_data,
    "moment_data": moment_data
}

json.dump(result_dict, open(results_filepath, mode="w"), indent=2)
