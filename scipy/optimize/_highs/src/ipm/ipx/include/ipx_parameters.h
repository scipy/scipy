#ifndef IPX_PARAMETERS_H_
#define IPX_PARAMETERS_H_

#include "ipx_config.h"

struct ipx_parameters {
    /* Solver control */
    ipxint display;
    const char* logfile;
    double print_interval;
    double time_limit;

    /* Preprocessing */
    ipxint dualize;
    ipxint scale;

    /* Interior point method */
    ipxint ipm_maxiter;
    double ipm_feasibility_tol;
    double ipm_optimality_tol;
    double ipm_drop_primal;
    double ipm_drop_dual;

    /* Linear solver */
    double kkt_tol;

    /* Basis construction in IPM */
    ipxint crash_basis;
    double dependency_tol;
    double volume_tol;
    ipxint rows_per_slice;
    ipxint maxskip_updates;

    /* LU factorization */
    ipxint lu_kernel;
    double lu_pivottol;

    /* Crossover */
    ipxint crossover;
    double crossover_start;
    double pfeasibility_tol;
    double dfeasibility_tol;

    /* Debugging */
    ipxint debug;
    ipxint switchiter;
    ipxint stop_at_switch;
    ipxint update_heuristic;
    ipxint maxpasses;

    #ifdef __cplusplus
    ipx_parameters() {
        display = 1;
        logfile = nullptr;
        print_interval = 5.0;
        time_limit = -1.0;
        dualize = -1;
        scale = 1;
        ipm_maxiter = 300;
        ipm_feasibility_tol = 1e-6;
        ipm_optimality_tol = 1e-8;
        ipm_drop_primal = 1e-9;
        ipm_drop_dual = 1e-9;
        kkt_tol = 0.3;
        crash_basis = 1;
        dependency_tol = 1e-6;
        volume_tol = 2.0;
        rows_per_slice = 10000;
        maxskip_updates = 10;
        lu_kernel = 0;
        lu_pivottol = 0.0625;
        crossover = 1;
        crossover_start = 1e-8;
        pfeasibility_tol = 1e-7;
        dfeasibility_tol = 1e-7;
        debug = 0;
        switchiter = -1;
        stop_at_switch = 0;
        update_heuristic = 1;
        maxpasses = -1;
    }
    #endif
};

#endif  /* IPX_PARAMETERS_H_ */
