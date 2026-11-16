#ifndef IPX_STATUS_H_
#define IPX_STATUS_H_

/* for status */
#define IPX_STATUS_not_run                  0
#define IPX_STATUS_solved                   1000
#define IPX_STATUS_stopped                  1005
#define IPX_STATUS_no_model                 1006
#define IPX_STATUS_out_of_memory            1003
#define IPX_STATUS_internal_error           1004

/* for status_ipm and status_crossover */
#define IPX_STATUS_optimal                  1
#define IPX_STATUS_imprecise                2
#define IPX_STATUS_primal_infeas            3
#define IPX_STATUS_dual_infeas              4
#define IPX_STATUS_user_interrupt           5
#define IPX_STATUS_time_limit               6
#define IPX_STATUS_iter_limit               7
#define IPX_STATUS_no_progress              8
#define IPX_STATUS_failed                   9
#define IPX_STATUS_debug                    10

/* error flags for invalid input */
#define IPX_ERROR_argument_null             102
#define IPX_ERROR_invalid_dimension         103
#define IPX_ERROR_invalid_matrix            104
#define IPX_ERROR_invalid_vector            105
#define IPX_ERROR_invalid_basis             107

/* error flags CR method */
#define IPX_ERROR_cr_iter_limit             201
#define IPX_ERROR_cr_matrix_not_posdef      202
#define IPX_ERROR_cr_precond_not_posdef     203
#define IPX_ERROR_cr_no_progress            204
#define IPX_ERROR_cr_inf_or_nan             205

/* error flags basis factorization */
#define IPX_ERROR_basis_singular            301
#define IPX_ERROR_basis_almost_singular     302
#define IPX_ERROR_basis_update_singular     303
#define IPX_ERROR_basis_repair_overflow     304
#define IPX_ERROR_basis_repair_search       305
#define IPX_ERROR_basis_too_ill_conditioned 306

#define IPX_ERROR_lapack_chol               401
#define IPX_ERROR_not_implemented           901
#define IPX_ERROR_user_interrupt            998
#define IPX_ERROR_time_interrupt            999

#define IPX_basic                            0
#define IPX_nonbasic                        -1
#define IPX_nonbasic_lb                     -1
#define IPX_nonbasic_ub                     -2
#define IPX_superbasic                      -3

#endif  /* IPX_STATUS_H_ */
