#ifndef IPX_INFO_H_
#define IPX_INFO_H_

#include "ipx_config.h"

struct ipx_info {
    ipxint status;
    ipxint status_ipm;
    ipxint status_crossover;
    ipxint errflag;

    /* dimension of LP model as given by user */
    ipxint num_var;
    ipxint num_constr;
    ipxint num_entries;

    /* dimension of constraint matrix in solver (including slack columns) */
    ipxint num_rows_solver;
    ipxint num_cols_solver;
    ipxint num_entries_solver;

    ipxint dualized;            /* dualized model? */
    ipxint dense_cols;          /* # columns classified "dense" */

    /* reductions in IPM */
    ipxint dependent_rows;      /* # dependent rows (to eq constr) removed */
    ipxint dependent_cols;      /* # dependent cols (to free vars) removed */
    ipxint rows_inconsistent;   /* dependent rows inconsistent? */
    ipxint cols_inconsistent;   /* dependent cols inconsistent? */
    ipxint primal_dropped;      /* # primal variables dropped to bound */
    ipxint dual_dropped;        /* # dual variables dropped to zero */

    /* interior solution */
    double abs_presidual;       /* max violation A*x+s=rhs, x-xl=lb, x+xu=ub */
    double abs_dresidual;       /* max violation A'y+zl-zu=obj */
    double rel_presidual;       /* relative primal residual */
    double rel_dresidual;       /* relative dual residual */
    double pobjval;             /* primal objective */
    double dobjval;             /* dual objective */
    double rel_objgap;          /* relative gap between primal and dual obj */
    double complementarity;     /* zl'xl + zu'xu + y'slack */
    double normx;               /* infnorm(x) */
    double normy;               /* infnorm(y) */
    double normz;               /* infnorm(zl,zu) */

    /* basic solution */
    double objval;              /* (primal) objective */
    double primal_infeas;       /* max violation of lb <= x <= ub */
    double dual_infeas;         /* max violation of sign condition on (y,z) */

    /* operation counts */
    ipxint iter;                /* # interior point iterations */
    ipxint kktiter1;            /* # linear solver iterations before switch */
    ipxint kktiter2;            /* # linear solver iterations after switch */
    ipxint basis_repairs;       /* # basis repairs after crash, < 0 discarded */
    ipxint updates_start;       /* # basis updates for starting basis */
    ipxint updates_ipm;         /* # basis updates in IPM */
    ipxint updates_crossover;   /* # basis updates in crossover */
    ipxint pushes_crossover;    /* # Primal/Dual pushes in crossover */

    /* major computation times */
    double time_total;          /* total runtime (wallclock) */
    double time_ipm1;           /* IPM before switch */
    double time_ipm2;           /* IPM after switch (without starting basis) */
    double time_starting_basis; /* constructing starting basis */
    double time_crossover;      /* crossover */

    /* profiling linear solver */
    double time_kkt_factorize;  /* factorize/build precond for KKT matrix */
    double time_kkt_solve;      /* linear system solves with KKT matrix */
    double time_maxvol;         /* algorithm "maxvolume" */
    double time_cr1;            /* CR method with diag precond */
    double time_cr1_AAt;        /* ... matrix-vector products with AA' */
    double time_cr1_pre;        /* ... preconditioner (diag+dense col) */
    double time_cr2;            /* CR method with basis precond */
    double time_cr2_NNt;        /* ... matrix-vector products with NN' */
    double time_cr2_B;          /* ... solves with B */
    double time_cr2_Bt;         /* ... solves with B' */

    /* profiling basis factorization */
    double ftran_sparse;        /* fraction of FTRAN solutions sparse */
    double btran_sparse;        /* fraction of BTRAN solutions sparse */
    double time_ftran;          /* all FTRAN operations */
    double time_btran;          /* all BTRAN operations */
    double time_lu_invert;      /* LU factorizations */
    double time_lu_update;      /* LU updates */
    double mean_fill;           /* geom. mean of LU fill factors */
    double max_fill;            /* max LU fill factor */
    double time_symb_invert;    /* computing row/column counts of inverse(B) */

    /* analysis of algorithm maxvolume */
    ipxint maxvol_updates;      /* # basis updates */
    ipxint maxvol_skipped;      /* # columns computed but basis not updated */
    ipxint maxvol_passes;       /* # passes over tableau matrix */
    ipxint tbl_nnz;             /* nnz in tbl matrix computed */
    double tbl_max;             /* max entry in tbl matrix computed */
    double frobnorm_squared;    /* Frobnorm^2 of tbl matrix computed */
    double lambdamax;           /* max eigenval of transformed normal matrix */
    double volume_increase;     /* base-2 log of volume(new)/volume(old) */

    #ifdef __cplusplus
    ipx_info();                 /* initializes all members to zero */
    #endif
};

#endif  /* IPX_INFO_H_ */
