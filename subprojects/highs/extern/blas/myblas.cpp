#include "mycblas.h"

#if defined(HIPO_USES_OPENBLAS)
#include <mutex>
#endif

// Returns the number of threads used by OpenBLAS after setting them, or -1 if
// OpenBLAS is not used. If num_threads <= 0, the number of threads is not
// changed and is just returned.
int highs_openblas_set_num_threads(int num_threads) {
#if defined(HIPO_USES_OPENBLAS)
  // openblas_set_num_threads is not thread safe: serialise calls from
  // concurrent HiGHS instances, and only call it if the number of threads
  // needs to change, so that OpenBLAS is not reconfigured while another
  // instance may be running BLAS routines
  static std::mutex openblas_threads_mutex;
  std::lock_guard<std::mutex> lock(openblas_threads_mutex);
  if (num_threads > 0 && openblas_get_num_threads() != num_threads)
    openblas_set_num_threads(num_threads);
  return openblas_get_num_threads();
#else
  return -1;
#endif
}
