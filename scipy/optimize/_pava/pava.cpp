#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>

#include <cstdint>
#include <string>
#include <tuple>

namespace nb = nanobind;

namespace {
using Array = nb::ndarray<double, nb::device::cpu, nb::c_contig>;
using IndexArray = nb::ndarray<intptr_t, nb::device::cpu, nb::c_contig>;

template <typename T>
static void check_ndim(const nb::ndarray<T, nb::device::cpu, nb::c_contig> &a) {
    if (a.ndim() != 1) {
        const std::string msg =
            "array has incorrect number of dimensions: " +
            std::to_string(a.ndim()) +
            "; expected 1";
        throw nb::value_error(msg.c_str());
    }
}

auto pava(Array xa, Array wa, IndexArray ra) {
    // x is the response variable (often written as y). Its ordering is crucial.
    // Usually, it is sorted according to some other data (feature or covariate), e.g.
    //   indices = np.argsort(z)
    //   x = x[indices]
    // Note that x is modified inplace and, on return, will contain the solution.
    // w is an array of case weights, modified inplace.
    // r is an array of indices such that x[r[i]:r[i+1]] contains the i-th block,
    // modified inplace.
    check_ndim(xa);
    check_ndim(wa);
    check_ndim(ra);
    auto x = xa.view<nb::ndim<1>>();
    auto w = wa.view<nb::ndim<1>>();
    auto r = ra.view<nb::ndim<1>>();

    intptr_t n = static_cast<intptr_t>(x.shape(0));

    // Algorithm 1 of
    // Busing, F. M. T. A. (2022).
    // Monotone Regression: A Simple and Fast O(n) PAVA Implementation.
    // Journal of Statistical Software, Code Snippets, 102(1), 1-25.
    // https://doi.org/10.18637/jss.v102.c01
    // Notes:
    //  - We translated it to 0-based indices.
    //  - xb, wb, sb instead of x, w and S to avoid name collisions
    //  - xb_prev and wb_prev instead of x_hat and w_hat
    //  - ERROR CORRECTED: Lines 9 and 10 have index i instead of b.
    //  - MODIFICATIONS: Lines 11 and 22 both have >= instead of >
    //    to get correct block indices in r. Otherwise, same values can get in
    //    different blocks, e.g. x = [2, 2] would produce
    //    r = [0, 1, 2] instead of r = [0, 2].
    //
    // procedure monotone(n, x, w)      // 1: x in expected order and w nonnegative
    r(0) = 0;  // 2: initialize index 0
    r(1) = 1;  // 3: initialize index 1
    intptr_t b = 0;  // 4: initialize block counter
    double xb_prev = x(0);  // 5: set previous block value
    double wb_prev = w(0);  // 6: set previous block weight

    for (intptr_t i = 1; i < n; ++i) {  // 7: loop over elements
        b++;  // 8: increase number of blocks
        double xb = x(i);  // 9: set current block value xb (i, not b)
        double wb = w(i);  // 10: set current block weight wb (i, not b)
        double sb = 0;

        if (xb_prev >= xb) {  // 11: check for down violation of x (>= instead of >)
            b--;  // 12: decrease number of blocks
            sb = wb_prev * xb_prev + wb * xb;  // 13: set current weighted block sum
            wb += wb_prev;  // 14: set new current block weight
            xb = sb / wb;  // 15: set new current block value

            while (i < n - 1 && xb >= x(i + 1)) {  // 16: repair up violations
                i++;
                sb += w(i) * x(i);  // 18: set new current weighted block sum
                wb += w(i);
                xb = sb / wb;
            }

            while (b > 0 && x(b - 1) >= xb) {  // 22: repair down violations (>= instead of >)
                b--;
                sb += w(b) * x(b);
                wb += w(b);
                xb = sb / wb;  // 26: set new current block value
            }
        }

        x(b) = xb_prev = xb;  // 29: save block value
        w(b) = wb_prev = wb;  // 30: save block weight
        r(b + 1) = i + 1;  // 31: save block index
    }

    intptr_t f = n - 1;  // 33: initialize "from" index
    for (intptr_t k = b; k >= 0; --k) {  // 34: loop over blocks
        intptr_t t = r(k);  // 35: set "to" index
        double xk = x(k);
        for (intptr_t i = f; i >= t; --i) {  // 37: loop "from" downto "to"
            x(i) = xk;  // 38: set all elements equal to block value
        }
        f = t - 1;  // 40: set new "from" equal to old "to" minus one
    }

    return std::make_tuple(xa, wa, ra, b + 1);  // b + 1 is number of blocks
}

NB_MODULE(_pava, m) {
    m.def(
        "pava",
        &pava,
        R"(Pool adjacent violators algorithm (PAVA) for isotonic regression

The routine might modify the input arguments x, w and r inplace.

Parameters
----------
xa : contiguous ndarray of shape (n,) and dtype np.float64
wa : contiguous ndarray of shape (n,) and dtype np.float64
ra : contiguous ndarray of shape (n+1,) and dtype np.intp

Returns
-------
x : ndarray
    The isotonic solution.
w : ndarray
    The array of weights for each block.
r : ndarray
    The array of indices for each block, such that xa[ra[i]:ra[i+1]]
    is the i-th block with all elements having the same value.
b : np.intp
    Number of blocks.
)",
        nb::arg("x"), nb::arg("w"), nb::arg("indices")
    );
}

}  // namespace
