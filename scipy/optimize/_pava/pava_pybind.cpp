#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <numpy/arrayobject.h>
#include <tuple>

namespace py = pybind11;

namespace {

auto pava(
    py::array_t<double, py::array::c_style | py::array::forcecast> xa,
    py::array_t<double, py::array::c_style | py::array::forcecast> wa,
    py::array_t<intptr_t, py::array::c_style | py::array::forcecast> ra
) {
    // x is the response variable (often written as y). Its ordering is crucial.
    // Usually, it is sorted according to some other data (feature or covariate), e.g.
    //   indices = np.argsort(z)
    //   x = x[indices]
    // Note that x is modified inplace and, on return, will contain the solution.
    // w is an array of case weights, modified inplace.
    // r is an array of indices such that x[r[i]:r[i+1]] contains the i-th block,
    // modified inplace.

    auto x = xa.mutable_unchecked<1>();
    intptr_t n = x.shape(0);
    auto w = wa.mutable_unchecked<1>();
    auto r = ra.mutable_unchecked<1>();

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
    r[0] = 0;  // 2: initialize index 0
    r[1] = 1;  // 3: initialize index 1
    intptr_t b = 0;  // 4: initialize block counter
    double xb_prev = x[b];  // 5: set previous block value
    double wb_prev = w[b];  // 6: set previous block weight
    for (intptr_t i = 1; i < n; ++i) {  // 7: loop over elements
        b++;  // 8: increase number of blocks
        double xb = x[i];  // 9: set current block value xb (i, not b)
        double wb = w[i];  // 10: set current block weight wb (i, not b)
        double sb = 0;
        if (xb_prev >= xb) {  // 11: check for down violation of x (>= instead of >)
            b--;  // 12: decrease number of blocks
            sb = wb_prev * xb_prev + wb * xb;  // 13: set current weighted block sum
            wb += wb_prev;  // 14: set new current block weight
            xb = sb / wb;  // 15: set new current block value
            while (i < n - 1 && xb >= x[i + 1]) {  // 16: repair up violations
                i++;
                sb += w[i] * x[i];  // 18: set new current weighted block sum
                wb += w[i];
                xb = sb / wb;
            }
            while (b > 0 && x[b - 1] >= xb) {  // 22: repair down violations (>= instead of >)
                b--;
                sb += w[b] * x[b];
                wb += w[b];
                xb = sb / wb;  // 26: set new current block value
            }
        }
        x[b] = xb_prev = xb;  // 29: save block value
        w[b] = wb_prev = wb;  // 30: save block weight
        r[b + 1] = i + 1;  // 31: save block index
    }

    intptr_t f = n - 1;  // 33: initialize "from" index
    for (intptr_t k = b; k >= 0; --k) {  // 34: loop over blocks
        intptr_t t = r[k];  // 35: set "to" index
        double xk = x[k];
        for (intptr_t i = f; i >= t; --i) {  // 37: loop "from" downto "to"
            x[i] = xk;  // 38: set all elements equal to block value
        }
        f = t - 1;  // 40: set new "from" equal to old "to" minus one
    }
    return std::make_tuple(xa, wa, ra, b + 1);  // b + 1 is number of blocks
}

PYBIND11_MODULE(_pava_pybind, m) {
    if (_import_array() != 0) {
        throw py::error_already_set();
    }
    m.def(
        "pava",
        &pava,
        "Pool adjacent violators algorithm (PAVA) for isotonic regression\n"
        "\n"
        "The routine might modify the input arguments x, w and r inplace.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "xa : contiguous ndarray of shape (n,) and dtype np.float64\n"
        "wa : contiguous ndarray of shape (n,) and dtype np.float64\n"
        "ra : contiguous ndarray of shape (n+1,) and dtype np.intp\n"
        "\n"
        "Returns\n"
        "-------\n"
        "x : ndarray\n"
        "    The isotonic solution.\n"
        "w : ndarray\n"
        "    The array of weights for each block.\n"
        "r : ndarray\n"
        "    The array of indices for each block, such that xa[ra[i]:ra[i+1]]\n"
        "    is the i-th block with all elements having the same value.\n"
        "b : np.intp\n"
        "    Number of blocks.\n",
        py::arg("x"), py::arg("w"), py::arg("indices")
    );
}

}  // namespace (anonymous)
