/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsCallbackStruct.h
 * @brief
 */
#ifndef LP_DATA_HIGHSCALLBACKSTRUCT_H_
#define LP_DATA_HIGHSCALLBACKSTRUCT_H_

#include "util/HighsType.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Struct to handle callback output data
 *
 */
typedef struct HighsCallbackDataOut {
  void* cbdata;  // cast of HighsCallbackOutput
  int log_type;  // cast of HighsLogType
  double running_time;
  HighsInt simplex_iteration_count;
  HighsInt ipm_iteration_count;
  HighsInt pdlp_iteration_count;
  HighsInt qpasm_iteration_count;
  double objective_function_value;
  int64_t mip_node_count;
  int64_t mip_total_lp_iterations;
  double mip_primal_bound;
  double mip_dual_bound;
  double mip_gap;
  double* mip_solution;
  HighsInt mip_solution_size;
  HighsInt cutpool_num_col;
  HighsInt cutpool_num_cut;
  HighsInt cutpool_num_nz;
  HighsInt* cutpool_start;
  HighsInt* cutpool_index;
  double* cutpool_value;
  double* cutpool_lower;
  double* cutpool_upper;
  double* qp_solution;
  HighsInt qp_solution_size;
  HighsInt external_solution_query_origin;
} HighsCallbackDataOut;

// Some external packages (e.g., jump) currently assume that the first 2 fields
// of this struct are interrupt and solution. Rearranging the struct may be a
// breaking change.
typedef struct HighsCallbackDataIn {
  int user_interrupt;
  double* user_solution;
  void* cbdata;  // cast of HighsCallbackInput (for internal use)
  int user_has_solution;
  HighsInt user_solution_size;
} HighsCallbackDataIn;

// Additional callback handling
typedef void (*HighsCCallbackType)(int, const char*,
                                   const HighsCallbackDataOut*,
                                   HighsCallbackDataIn*, void*);

typedef HighsInt (*HighsCHessianFunctionType)(const HighsInt, const HighsInt*,
                                              const HighsInt*, const double*,
                                              HighsInt*, HighsInt*, double*,
                                              void*);

#ifdef __cplusplus
}
#endif

#endif /* LP_DATA_HIGHSCALLBACKSTRUCT_H_ */
