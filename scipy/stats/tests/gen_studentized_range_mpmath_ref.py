# To run this script, use
# `python gen_studentized_range_mpmath_ref.py`

# This script generates pickled file "./data/studentized_range_mpmath_ref.p"
# that is used to test the accuracy of `studentized_range` functions against
# arbitrarily precise results generated using `mpmath`.

# Note: I would have prefered to use pickle rather than JSON, but -
# due to security concerns - decided against it.

from collections import namedtuple
import json
import time

from mpmath import gamma, pi, erf, exp, sqrt, quad, inf, mpf, mp
from mpmath import npdf as phi
from mpmath import ncdf as Phi

results_filepath = "data/studentized_range_mpmath_ref.json"

MPResult = namedtuple("MPResult", ["src_case", "mp_result"])

CdfCase = namedtuple("CdfCase",
                     ["q", "k", "v", "expected_atol", "expected_rtol"])

general_atol = 1e-11
general_rtol = 1e-11

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
    CdfCase(10, 3, 1000, general_atol, general_rtol),
]


def to_dict(named_tuple):
    return dict(named_tuple._asdict())


def mp_res_to_dict(mp_result):
    """Formats and MPResults tuple into a dict for JSON dumping"""
    return {
        "src_case": to_dict(mp_result.src_case),

        # np assert can't handle mpf, so take the accuracy hit here.
        "mp_result": float(mp_result.mp_result)
    }


def cdf_mp(q, k, nu):
    """Calculates an integrand for the given q, k, nu vals"""
    mp.dps = 24
    q, k, nu = mpf(q), mpf(k), mpf(nu)

    def inner(s, z):
        return phi(z) * (Phi(z + q * s) - Phi(z)) ** (k - 1)

    def outer(s, z):
        return s ** (nu - 1) * phi(sqrt(nu) * s) * inner(s, z)

    def whole(s, z):
        return sqrt(2 * pi) * k * nu ** (nu / 2) / (
                gamma(nu / 2) * 2 ** (nu / 2 - 1)) * outer(s, z)

    res = quad(whole, [0, inf], [-inf, inf],
               method="gauss-legendre", maxdegree=10)
    return res


print(f"Processing {len(cdf_cases)} CDF cases. This should take about "
      f"{0.5 * len(cdf_cases)} minutes at 30s per case on a "
      f"midrange CPU.")

t_start = time.perf_counter()
# [mp_res_to_dict(MPResult(cdf_cases[0], mpf(2)))]
cdf_data = [mp_res_to_dict(MPResult(case, cdf_mp(case.q, case.k, case.v)))
            for case in cdf_cases]

print(f"CDF data generated in {time.perf_counter() - t_start}s")

result_dict = {
    "cdf_data": cdf_data
}
print(result_dict)
json.dump(result_dict, open(results_filepath, mode="w"), indent=2)
