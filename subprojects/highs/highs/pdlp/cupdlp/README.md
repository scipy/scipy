# cuPDLP-C observations

This directory contains files from [cuPDLP-C v0.3.0](https://github.com/COPT-Public/cuPDLP-C/tree/v0.3.0). Below are some issues experienced when integrating them into HiGHS. In March 2025 it was updated to https://github.com/COPT-Public/cuPDLP-C/commit/7b94c41efca42fa176a129a935b59718e1455fad. 

## Termination of cuPDLP-C

cuPDLP-C terminates when either the current or averaged iterates satisfy primal/dual feasibility, using a 2-norm measure relative to the size of the RHS/costs, and after scaling the LP. 

HiGHS assesses primal/dual feasibility using a infinity-norm absolute measure for the unscaled LP. Thus the cuPDLP-C result frequently fails to satisfy HiGHS primal/dual feasibility. To get around this partially, `iInfNormAbsLocalTermination` has been introduced into cuPDLP-C. 

By default, `iInfNormAbsLocalTermination` is false, so that the original cuPDLP-C termination criteria are used.

When `iInfNormAbsLocalTermination` is true, cuPDLP-C terminates only when primal/dual feasibility is satisfied for the infinity-norm absolute measure of the current iterate, so that HiGHS primal/dual feasibility is satisfied. 

However, the cuPDLP-C scaling may still result in the HiGHS tolerances not being satisfied. Users can inspect `HighsInfo` values for the maximum and sum of infeasibilities, and the new `HighsInfo` values measuring the maximum and sum of complementarity violations.

## Preprocessing issue

The following line is not recognised by g++, 

> #if !defined(CUPDLP_CPU)

so I've had to replace all ocurrences by

> #ifndef CUPDLP_CPU

This yields a compiler warning about "extra tokens at end of #ifndef
directive" in the case of the following, but it's not a problem for
now, as CUPDLP_CPU is set

> #ifndef CUPDLP_CPU & USE_KERNELS

## Use of macro definitions within C++

When definitions in `glbopts.h` are used in `CupdlpWrapper.cpp` there is a g++ compiler error, because `typeof` isn't recognised

> #define CUPDLP_INIT(var, size)                                  \
  {                                                             \
    (var) = (typeof(var))malloc((size) * sizeof(typeof(*var))); \
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }

Hence there is a set of type-specific definitions in `CupdlpWrapper.h`, such as 

>#define cupdlp_init_double(var, size)\
   {\
     (var) = (double*)malloc((size) * sizeof(double));\
   }

## C methods not picked up by g++

Three methods
* `double infNorm(double *x, cupdlp_int n);`
* `void cupdlp_haslb(cupdlp_float *haslb, const cupdlp_float *lb, const cupdlp_float bound, const cupdlp_int len);`
* `void cupdlp_hasub(cupdlp_float *hasub, const cupdlp_float *ub, const cupdlp_float bound, const cupdlp_int len);`

are declared in `cupdlp_linalg.h` and defined in `cupdlp_linalg.c` but not picked up by g++. Hence duplicate methods are declared and defined in `CupdlpWrapper.h` and `CupdlpWrapper.cpp`.

## Use of macro definitions within C

Although the macro definitions in `glbopts.h` are fine when used in C under Linux, they cause the following compiler errors on Windows.

> error C2146: syntax error: missing ';' before identifier 'calloc' (or 'malloc')

In HiGHS, all the macros using `typeof` have been replaced by multiple type-specific macros

## Problem with sys/time.h

The HiGHS branch add-pdlp compiles and runs fine on @jajhall's Linux machine, but CI tests on GitHub fail utterly due to `sys/time.h` not being found. Until this is fixed, or HiGHS passes its own timer for use within `cuPDLP-c`, timing within `cuPDLP-c` can be disabled using the compiler directive `CUPDLP_TIMER`. By default this is defined, so the `cuPDLP-c` is retained.

## Controlling the `cuPDLP-c` logging

As a research code, `cuPDLP-c` naturally produces a lot of logging output. HiGHS must be able to run with less logging output, or completely silently. This is achieved using the `nLogLevel` parameter in `cuPDLP-c`. 

By default, `nLogLevel` is 2, so all the original `cuPDLP-c` logging is produced.

* If `nLogLevel` is 1, then the `cuPDLP-c` logging is less verbose 
* If `nLogLevel` is 0, then there is no `cuPDLP-c` logging

A related issue is the use of `fp` and `fp_sol`. HiGHS won't be using these, so sets them to null pointers. `cuPDLP-c` already doesn't print the solution if `fp_sol` is a null pointer, so the call to `writeJson(fp, pdhg);` is now conditional on `if (fp)`. 

## Handling infeasible or unbounded problems

`cuPDLP-c` now terminates with status `INFEASIBLE_OR_UNBOUNDED` for the infeasible and unbounded LPs in unit tests `pdlp-infeasible-lp` and `pdlp-unbounded-lp` in `highs/check/TestPdlp.cpp`. In the case of the unbounded LP, PDLP identifies a primal feasible point, so unboundedness can be deduced. This is done in `HighsSolve.cpp:131.

## Returning the iteration count

The `cuPDLP-c` iteration count is held in `pdhg->timers->nIter`, but `pdhg` is destroyed in `LP_SolvePDHG`, so `cupdlp_int* num_iter` has been added to the parameter list of this method.

## To be done

- Make CupldlpWrapper.cpp look more like C++ than C
