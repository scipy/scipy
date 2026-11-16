/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSort.h
 * @brief Sorting routines for HiGHS
 */
#ifndef UTIL_HIGHSSORT_H_
#define UTIL_HIGHSSORT_H_

#include <vector>

#include "lp_data/HConst.h"

using std::vector;

void addToDecreasingHeap(HighsInt& n, HighsInt mx_n, vector<double>& heap_v,
                         vector<HighsInt>& heap_ix, double v, HighsInt ix);
void sortDecreasingHeap(const HighsInt n, vector<double>& heap_v,
                        vector<HighsInt>& heap_ix);
/**
 * @brief Sort values[1..n] of an array by increasing value
 */
void maxheapsort(HighsInt* heap_v,  //!< HighsInt values to be sorted
                 HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Sort values[1..n] of an array by increasing value with corresponding
 * indices
 */
void maxheapsort(
    HighsInt* heap_v,  //!< Values to be sorted
    HighsInt* heap_i,  //!< Indices corresponding to (sorted) values
    HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Sort values[1..n] of an array by increasing value with corresponding
 * indices
 */
void maxheapsort(
    double* heap_v,    //!< Values to be sorted
    HighsInt* heap_i,  //!< Indices corresponding to (sorted) values
    HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Build a value heap for sorting values[1..n] of an array by increasing
 * value
 */
void buildMaxheap(HighsInt* heap_v,  //!< HighsInt values to be sorted
                  HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Build a value-index heap for sorting values[1..n] of an array by
 * increasing value
 */
void buildMaxheap(
    HighsInt* heap_v,  //!< Values to be sorted
    HighsInt* heap_i,  //!< Indices corresponding to (sorted) values
    HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Build a value-index heap for sorting values[1..n] of an array by
 * increasing value
 */
void buildMaxheap(
    double* heap_v,    //!< Values to be sorted
    HighsInt* heap_i,  //!< Indices corresponding to (sorted) values
    HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Sort by increasing value a heap built with buildMaxheap
 */
void maxHeapsort(HighsInt* heap_v,  //!< HighsInt values to be sorted
                 HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Sort by increasing value a heap built with buildMaxheap
 */
void maxHeapsort(
    HighsInt* heap_v,  //!< Values to be sorted
    HighsInt* heap_i,  //!< Indices corresponding to (sorted) values
    HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Sort by increasing value a heap built with buildMaxheap
 */
void maxHeapsort(
    double* heap_v,    //!< Values to be sorted
    HighsInt* heap_i,  //!< Indices corresponding to (sorted) values
    HighsInt n         //!< Number of values to be sorted
);
/**
 * @brief Heapify function for sorting by increasing value
 */
void maxHeapify(HighsInt* heap_v, HighsInt i, HighsInt n);

/**
 * @brief Heapify function for sorting by increasing value
 */
void maxHeapify(HighsInt* heap_v, HighsInt* heap_i, HighsInt i, HighsInt n);

/**%
 * @brief Heapify function for sorting by increasing value
 */
void maxHeapify(double* heap_v, HighsInt* heap_i, HighsInt i, HighsInt n);

/**
 * @brief Check that a set of integers is in increasing order and in bounds
 */
bool increasingSetOk(const vector<HighsInt>& set,
                     const HighsInt set_entry_lower,
                     const HighsInt set_entry_upper, bool strict);

/**
 * @brief Check that a set of doubles is in increasing order and in bounds
 */
bool increasingSetOk(const vector<double>& set, const double set_entry_lower,
                     const double set_entry_upper, bool strict);

void sortSetData(const HighsInt num_entries, vector<HighsInt>& set,
                 const double* data0, const double* data1, const double* data2,
                 double* sorted_data0, double* sorted_data1,
                 double* sorted_data2);

void sortSetData(const HighsInt num_entries, vector<HighsInt>& set,
                 const HighsVarType* data0, HighsVarType* sorted_data0);

#endif /* UTIL_HIGHSSORT_H_ */
