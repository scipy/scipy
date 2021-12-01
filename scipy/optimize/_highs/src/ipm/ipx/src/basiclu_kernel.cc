// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "basiclu_kernel.h"
#include <cassert>
#include <new>                  // bad_alloc
#include <stdexcept>
#include "basiclu.h"

namespace ipx {

// Wraps the BASICLU object into a struct with contructor/destructor for RAII.
struct BasicLuHelper {
    static_assert(sizeof(Int) == sizeof(lu_int),
                  "IPX integer type does not match BASICLU integer type");

    BasicLuHelper(Int dim) {
        Int err = basiclu_obj_initialize(&obj, dim);
        if (err == BASICLU_ERROR_out_of_memory)
            throw std::bad_alloc();
        if (err != BASICLU_OK)
            throw std::logic_error("basiclu_obj_initialize failed");
    }

    ~BasicLuHelper() {
        basiclu_obj_free(&obj);
    }

    struct basiclu_object obj;
};

void BasicLuKernel::_Factorize(Int dim, const Int* Bbegin, const Int* Bend,
                               const Int* Bi, const double* Bx, double pivottol,
                               bool strict_abs_pivottol,
                               SparseMatrix* L, SparseMatrix* U,
                               std::vector<Int>* rowperm,
                               std::vector<Int>* colperm,
                               std::vector<Int>* dependent_cols) {
    BasicLuHelper lu(dim);
    lu.obj.xstore[BASICLU_REL_PIVOT_TOLERANCE] = pivottol;
    if (strict_abs_pivottol) {
        lu.obj.xstore[BASICLU_ABS_PIVOT_TOLERANCE] = kLuDependencyTol;
        lu.obj.xstore[BASICLU_REMOVE_COLUMNS] = 1;
    }
    Int err = basiclu_obj_factorize(&lu.obj, Bbegin, Bend, Bi, Bx);
    if (err == BASICLU_ERROR_out_of_memory)
        throw std::bad_alloc();
    if (err != BASICLU_OK && err != BASICLU_WARNING_singular_matrix)
        throw std::logic_error("basiclu_obj_factorize failed");

    // Extract factors. Dependent columns are at the end of the BASICLU pivot
    // sequence.
    Int rank = lu.obj.xstore[BASICLU_RANK];
    dependent_cols->clear();
    for (Int k = rank; k < dim; k++)
        dependent_cols->push_back(k);
    L->resize(dim, dim, dim + lu.obj.xstore[BASICLU_LNZ]);
    U->resize(dim, dim, dim + lu.obj.xstore[BASICLU_UNZ]);
    rowperm->resize(dim);
    colperm->resize(dim);
    err = basiclu_obj_get_factors(&lu.obj, rowperm->data(), colperm->data(),
                                  L->colptr(), L->rowidx(), L->values(),
                                  U->colptr(), U->rowidx(), U->values());
    if (err != BASICLU_OK)
        throw std::logic_error("basiclu_obj_get_factors failed");

    // Remove unit diagonal from L.
    Int num_dropped = RemoveDiagonal(*L, nullptr);
    assert(num_dropped == dim);
    assert(L->entries() == lu.obj.xstore[BASICLU_LNZ]);
    (void)(num_dropped);
}

}  // namespace ipx
