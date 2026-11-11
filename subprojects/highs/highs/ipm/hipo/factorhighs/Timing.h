#ifndef FACTORHIGHS_TIMING_H
#define FACTORHIGHS_TIMING_H

#include "FactorHighsSettings.h"

namespace hipo {

enum TimeItems {
  kTimeAnalyse,                           // TIMING_LEVEL 1
  kTimeAnalyseTree,                       // TIMING_LEVEL 2
  kTimeAnalyseCount,                      // TIMING_LEVEL 2
  kTimeAnalysePattern,                    // TIMING_LEVEL 2
  kTimeAnalyseSn,                         // TIMING_LEVEL 2
  kTimeAnalyseReorder,                    // TIMING_LEVEL 2
  kTimeAnalyseRelInd,                     // TIMING_LEVEL 2
  kTimeFactorise,                         // TIMING_LEVEL 1
  kTimeFactorisePrepare,                  // TIMING_LEVEL 2
  kTimeFactoriseAssembleOriginal,         // TIMING_LEVEL 2
  kTimeFactoriseAssembleChildrenFrontal,  // TIMING_LEVEL 2
  kTimeFactoriseAssembleChildrenClique,   // TIMING_LEVEL 2
  kTimeFactoriseDenseFact,                // TIMING_LEVEL 2
  kTimeDenseFact_main,                    // TIMING_LEVEL 2
  kTimeDenseFact_schur,                   // TIMING_LEVEL 2
  kTimeDenseFact_kernel,                  // TIMING_LEVEL 2
  kTimeDenseFact_convert,                 // TIMING_LEVEL 2
  kTimeDenseFact_pivoting,                // TIMING_LEVEL 2
  kTimeFactoriseTerminate,                // TIMING_LEVEL 2
  kTimeSolve,                             // TIMING_LEVEL 1
  kTimeSolvePrepare,                      // TIMING_LEVEL 2
  kTimeSolveSolve,                        // TIMING_LEVEL 2
  kTimeSolveSolve_dense,                  // TIMING_LEVEL 2
  kTimeSolveSolve_sparse,                 // TIMING_LEVEL 2
  kTimeSolveSolve_swap,                   // TIMING_LEVEL 2
  kTimeBlasStart,                         //
  kTimeBlas_copy = kTimeBlasStart,        // TIMING_LEVEL 3
  kTimeBlas_axpy,                         // TIMING_LEVEL 3
  kTimeBlas_scal,                         // TIMING_LEVEL 3
  kTimeBlas_swap,                         // TIMING_LEVEL 3
  kTimeBlas_gemv,                         // TIMING_LEVEL 3
  kTimeBlas_trsv,                         // TIMING_LEVEL 3
  kTimeBlas_tpsv,                         // TIMING_LEVEL 3
  kTimeBlas_ger,                          // TIMING_LEVEL 3
  kTimeBlas_trsm,                         // TIMING_LEVEL 3
  kTimeBlas_syrk,                         // TIMING_LEVEL 3
  kTimeBlas_gemm,                         // TIMING_LEVEL 3
  kTimeBlasEnd = kTimeBlas_gemm,          //
  kTimeSize                               //
};

// Macros for timing.
// - Use HIPO_CLOCK_CREATE to create the Clock objects.
// - Use HIPO_CLOCK_START(id) to start the corresponding clock.
// - Use HIPO_CLOCK_STOP(id, data, item) to read the corresponding clock and
//   save the time in the corresponding time item of the data collector.
// If HIPO_TIMING_LEVEL is not >= id (at compile time), then nothing happens.
//
#define HIPO_CLOCK_CREATE \
  HIPO_CLOCK_CREATE_1;    \
  HIPO_CLOCK_CREATE_2;    \
  HIPO_CLOCK_CREATE_3;
#define HIPO_CLOCK_START(id) HIPO_CLOCK_START_##id
#define HIPO_CLOCK_STOP(id, data, item) HIPO_CLOCK_STOP_##id(data, item)

#if HIPO_TIMING_LEVEL >= 1
#define HIPO_CLOCK_CREATE_1 Clock clock_1
#else
#define HIPO_CLOCK_CREATE_1
#endif
#if HIPO_TIMING_LEVEL >= 2
#define HIPO_CLOCK_CREATE_2 Clock clock_2
#else
#define HIPO_CLOCK_CREATE_2
#endif
#if HIPO_TIMING_LEVEL >= 3
#define HIPO_CLOCK_CREATE_3 Clock clock_3
#else
#define HIPO_CLOCK_CREATE_3
#endif

#if HIPO_TIMING_LEVEL >= 1
#define HIPO_CLOCK_START_1 clock_1.start()
#else
#define HIPO_CLOCK_START_1
#endif
#if HIPO_TIMING_LEVEL >= 2
#define HIPO_CLOCK_START_2 clock_2.start()
#else
#define HIPO_CLOCK_START_2
#endif
#if HIPO_TIMING_LEVEL >= 3
#define HIPO_CLOCK_START_3 clock_3.start()
#else
#define HIPO_CLOCK_START_3
#endif

#if HIPO_TIMING_LEVEL >= 1
#define HIPO_CLOCK_STOP_1(data, item) (data).sumTime(item, clock_1.stop())
#else
#define HIPO_CLOCK_STOP_1(data, item) (void)data  // suppress warning
#endif
#if HIPO_TIMING_LEVEL >= 2
#define HIPO_CLOCK_STOP_2(data, item) (data).sumTime(item, clock_2.stop())
#else
#define HIPO_CLOCK_STOP_2(data, item)
#endif
#if HIPO_TIMING_LEVEL >= 3
#define HIPO_CLOCK_STOP_3(data, item) (data).sumTime(item, clock_3.stop())
#else
#define HIPO_CLOCK_STOP_3(data, item)
#endif

}  // namespace hipo

#endif