/*
 * C++ implementation of column grouping for finite difference Jacobian
 * estimation.
 *
 * Greedy sequential algorithm from:
 *   A. Curtis, M. J. D. Powell, and J. Reid,
 *   "On the estimation of sparse Jacobian matrices",
 *   Journal of the Institute of Mathematics and its Applications,
 *   13 (1974), pp. 117-120.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <numpy/arrayobject.h>
#include <vector>

namespace py = pybind11;

namespace {

template <typename index_t>
py::array_t<npy_intp> group_dense_impl(
    int m, int n,
    py::array_t<index_t, py::array::c_style> A_arr
) {
    const index_t *A = A_arr.data();

    auto groups_arr = py::array_t<npy_intp>(n);
    npy_intp *groups = groups_arr.mutable_data();
    for (int i = 0; i < n; ++i) {
        groups[i] = -1;
    }

    std::vector<int> u(m);

    int current_group = 0;

    for (int i = 0; i < n; ++i) {
        if (groups[i] >= 0) {
            continue;
        }

        groups[i] = current_group;
        bool all_grouped = true;

        for (int k = 0; k < m; ++k) {
            u[k] = A[static_cast<npy_intp>(k) * n + i];
        }

        for (int j = 0; j < n; ++j) {
            if (groups[j] < 0) {
                all_grouped = false;
            } else {
                continue;
            }

            bool intersect = false;
            for (int k = 0; k < m; ++k) {
                if (u[k] > 0 && A[static_cast<npy_intp>(k) * n + j] > 0) {
                    intersect = true;
                    break;
                }
            }

            if (!intersect) {
                for (int k = 0; k < m; ++k) {
                    u[k] += A[static_cast<npy_intp>(k) * n + j];
                }
                groups[j] = current_group;
            }
        }

        if (all_grouped) {
            break;
        }

        current_group += 1;
    }

    return groups_arr;
}

py::array_t<npy_intp> group_dense(
    int m, int n,
    py::array A_arr
) {
    auto dtype = A_arr.dtype();

    if (dtype.is(py::dtype::of<int32_t>())) {
        return group_dense_impl<int32_t>(m, n, A_arr);
    } else if (dtype.is(py::dtype::of<int64_t>())) {
        return group_dense_impl<int64_t>(m, n, A_arr);
    } else {
        throw py::type_error("A must have dtype int32 or int64");
    }
}

template <typename index_t>
py::array_t<npy_intp> group_sparse_impl(
    int m, int n,
    py::array_t<index_t> indices_arr,
    py::array_t<index_t> indptr_arr
) {
    const index_t *indices = indices_arr.data();
    const index_t *indptr = indptr_arr.data();

    auto groups_arr = py::array_t<npy_intp>(n);
    npy_intp *groups = groups_arr.mutable_data();
    for (int i = 0; i < n; ++i) {
        groups[i] = -1;
    }

    std::vector<int> u(m, 0);

    int current_group = 0;

    for (int i = 0; i < n; ++i) {
        if (groups[i] >= 0) {
            continue;
        }

        groups[i] = current_group;
        bool all_grouped = true;

        for (int k = 0; k < m; ++k) {
            u[k] = 0;
        }
        for (index_t k = indptr[i]; k < indptr[i + 1]; ++k) {
            u[indices[k]] = 1;
        }

        for (int j = 0; j < n; ++j) {
            if (groups[j] < 0) {
                all_grouped = false;
            } else {
                continue;
            }

            bool intersect = false;
            for (index_t k = indptr[j]; k < indptr[j + 1]; ++k) {
                if (u[indices[k]] == 1) {
                    intersect = true;
                    break;
                }
            }
            if (!intersect) {
                for (index_t k = indptr[j]; k < indptr[j + 1]; ++k) {
                    u[indices[k]] = 1;
                }
                groups[j] = current_group;
            }
        }

        if (all_grouped) {
            break;
        }

        current_group += 1;
    }

    return groups_arr;
}

py::array_t<npy_intp> group_sparse(
    int m, int n,
    py::array indices_arr,
    py::array indptr_arr
) {
    auto dtype = indices_arr.dtype();

    if (dtype.is(py::dtype::of<int32_t>())) {
        return group_sparse_impl<int32_t>(
            m, n,
            py::cast<py::array_t<int32_t>>(indices_arr),
            py::cast<py::array_t<int32_t>>(indptr_arr)
        );
    } else if (dtype.is(py::dtype::of<int64_t>())) {
        return group_sparse_impl<int64_t>(
            m, n,
            py::cast<py::array_t<int64_t>>(indices_arr),
            py::cast<py::array_t<int64_t>>(indptr_arr)
        );
    } else {
        throw py::type_error(
            "indices and indptr must have dtype int32 or int64"
        );
    }
}

PYBIND11_MODULE(_group_columns, m, py::mod_gil_not_used()) {
    if (_import_array() != 0) {
        throw py::error_already_set();
    }
    m.def("group_dense", &group_dense,
        py::arg("m"), py::arg("n"), py::arg("A"));
    m.def("group_sparse", &group_sparse,
        py::arg("m"), py::arg("n"), py::arg("indices"), py::arg("indptr"));
}

}  // namespace (anonymous)
